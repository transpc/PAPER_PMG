# Paper Draft: An Unstructured-Grid PMG-Preconditioned BiCGSTAB Solver with Application to MSFR and iSMR

> Date: 2026-06-10 (draft v0.1)
> Status: Skeleton + draft of Introduction / Related Work / Methodology. Results section to be filled after experiments.
> Related note: [[Literature Survey - PMG-Preconditioned BiCGSTAB]]

---

# 0. Meta Information

- **Working title (English)**: *A Parallel Semi-Coarsening Geometric Multigrid Preconditioner for BiCGSTAB on Unstructured Grids: Application to MSFR Steady State and iSMR Natural Circulation*
- **Target journals**: *Annals of Nuclear Energy* / *Nuclear Engineering and Design* / *Journal of Computational Physics* (computational nuclear engineering)
- **Key contributions**:
  1. Integration of a semi-coarsening geometric MG applicable to unstructured grids as a **preconditioner for BiCGSTAB**. Whereas the preceding meshless GMG work [@do2024meshless] was a standalone solver, this work extends it to a preconditioner, simultaneously securing robustness and acceleration.
  2. Parallel scalability secured through **communication-volume reduction** and a **serial-computation-below-a-fixed-level** strategy.
  3. Demonstration of speedup and robustness on the **MSFR steady state**, **iSMR ECT natural circulation**, and **iSMR core daily load-following transient** problems.

---

# 1. Introduction

## 1-1. Background

In reactor thermal-hydraulic analysis, solving the pressure correction equation (pressure correction / Poisson type) for incompressible flow dominates the overall computation time. This elliptic equation has a large condition number, and the convergence of iterative solvers degrades sharply with mesh refinement, making it the bottleneck of large-scale 3D analysis.

KAERI's CUPID code analyzes reactor cores and systems with an unstructured-grid-based two-fluid, three-field model, and has traditionally used BiConjugate Gradient (BiCG)-family solvers for the pressure correction equation [@do2024meshless]. However, BiCG-family methods have an intrinsic limitation: the iteration count grows with mesh size (non-scalable).

Multigrid (MG) is almost the only family of methods that theoretically provides a convergence rate independent of mesh size. In prior work [@do2024meshless; @do2022highaspect; @ha2024aiaa], the authors' group developed a **meshless geometric multigrid (GMG)** directly applicable to unstructured grids and demonstrated superior scalability compared with BiCG. This work **reconstructs it as a preconditioner** and extends it with **parallel communication reduction** and application to **MSFR and iSMR**.

## 1-2. Motivation: Why MG-Preconditioned BiCGSTAB?

- **Limitations of a standalone MG solver**: Geometric MG is sensitive to the smoothing and grid-transfer operators with respect to the problem, so its convergence may stagnate when used alone on strongly anisotropic, irregular unstructured grids.
- **Effect of Krylov acceleration**: When MG is used as a preconditioner for BiCGSTAB, MG efficiently removes low-frequency error while BiCGSTAB handles the residual spectral components that MG misses, greatly improving robustness.
- **The variable-preconditioner issue**: In large-scale parallel runs, if the MG hierarchy is restricted and the coarsest grid is approximated by an iterative method to reduce communication, the preconditioner changes at every iteration, becoming a **variable preconditioner**. In this case standard BiCGSTAB may break down in convergence, and a **flexible variant (FBiCGSTAB)** is required [@chen2016fbicgstab].

## 1-3. Objectives

This work aims to (1) integrate an unstructured-grid semi-coarsening geometric MG preconditioner into BiCGSTAB (FBiCGSTAB when necessary), (2) implement communication-volume reduction and serial computation below a fixed level, and (3) quantitatively demonstrate the acceleration performance on the MSFR steady state, iSMR ECT natural circulation, and iSMR core daily load-following transient problems.

---

# 2. Related Work

## 2-1. Unstructured-Grid Geometric Multigrid (Authors' Group Prior Work)

The **scalable solver** lineage that this work directly extends is as follows.

- **Meshless GMG (ANE 2024)** [@do2024meshless]: Introduced a node-coarsening-based meshless geometric MG into the CUPID code. It converges stably on both structured and unstructured grids and substantially reduces computation time compared with BiCG. **This is the starting point of the present paper and the basis for the unstructured-grid application strategy.**
- **High aspect ratio GMG (JMST 2022)** [@do2022highaspect]: Meshless GMG on high-aspect-ratio (stretched), anisotropic grids. Corresponds to boundary layers and thin-film flows, and connects directly to the semi-coarsening motivation.
- **Complex geometry cell-coarsening (AIAA J.)** [@ha2024aiaa]: An improved cell-coarsening algorithm for complex geometries. The basis for application to complex shapes such as MSFR and iSMR.
- **Node-coarsening for linear FE (2021)** [@ha2021nodecoarsening]: The foundation of the node-coarsening algorithm for linear finite element discretization.
- **Optimal aggregation level for parallel meshless MG (DDM)** [@ha2025aggregation]: A study of the optimal aggregation level for a domain-decomposition (DDM)-based parallel meshless MG — **directly connected to this paper's "serial computation below a fixed level" and communication-reduction strategy.**

## 2-2. Semi-Coarsening Geometric Multigrid

Semi-coarsening coarsens only one coordinate direction and is robust to the anisotropy of stretched grids [@osti440722]. Since Mulder's proposal [@mulder1989], MSG (Multiple Semi-Coarsened Multigrid) [@dendy2019msg] added extra coarse grids coarsened along each coordinate direction, improving the treatment of anisotropic diffusion. The fact that coordinate-alignment directions are ambiguous on unstructured grids is the central challenge and contribution point of this work.

## 2-3. Parallel MG Communication Reduction and Coarse-Level Handling

The central challenge of large-scale parallel MG is coarse-level communication. Key strategies:

- **Coarse-grid redistribution (incremental agglomeration)** [@reisner2018]: Gradually reduces processes using a predictive performance model, scaling structured-grid MG to 500K+ cores.
- **Node-aware multi-stage communication** [@bienz2020]: Simultaneously reduces the number and size of inter-node messages, prioritizing intra-node channels.
- **Post-hierarchy sparsification** [@manteuffel2016]: Reduces communication by removing coarse-matrix entries after hierarchy formation, with a convergence-degradation trade-off when overdone.

These have each been validated on structured grids / AMG, and **their combined application to unstructured reactor CFD grids remains an open area.**

## 2-4. BiCGSTAB and Variable Preconditioning

BiCGSTAB [@vandervorst1992] is a standard Krylov solver more efficient than CGS for nonsymmetric systems. When MG is used as a variable preconditioner, restricting the hierarchy to about three levels and approximating the coarsest grid iteratively, then applying **FBiCGSTAB**, has been reported to yield a 2–3× speedup over Chebyshev-preconditioned IBiCGSTAB at large scale [@chen2016fbicgstab]. There is also a theoretical basis that Krylov convergence accelerates when the spectrum of the preconditioned operator clusters at a few points [@lucerolorca2026].

## 2-5. Numerical Analysis of MSFR / SMR Natural Circulation

In coupled neutronics/thermal-hydraulics MSFR analysis [@msfr_coupled], the fuel salt circulates the primary loop in roughly 3–4 seconds, which is directly tied to delayed neutron precursor (DNP) transport and is a reference design parameter. SMR natural circulation is a buoyancy-driven, low-Reynolds, thermally stratified condition that is qualitatively different from MSFR forced turbulence, requiring a separate examination of preconditioner robustness.

---

# 3. Methodology

## 3-1. Governing Equations and Discretization

The pressure correction equation derived from SIMPLE-family pressure-velocity coupling in CUPID's incompressible (or low-Mach) two-fluid model:

$$
\nabla \cdot \left( \frac{1}{a_P} \nabla p' \right) = \nabla \cdot \mathbf{u}^*
$$

Unstructured finite-volume discretization yields a nonsymmetric (or weakly symmetric) sparse linear system $A p' = b$. The condition number of $A$ deteriorates with grid aspect ratio and non-orthogonality.

## 3-2. Semi-Coarsening Geometric MG Preconditioner (Unstructured Grids)

Based on the node-coarsening of the preceding meshless GMG [@do2024meshless; @ha2021nodecoarsening], we construct a semi-coarsening variant that detects the anisotropy direction and coarsens only along that direction.

- **Direction detection**: Since there are no coordinate axes on an unstructured grid, the "strong-coupling direction" is estimated from the local grid aspect ratio and connection strength (matrix coefficient magnitude).
- **Grid-transfer operators**: Meshless interpolation (restriction/prolongation) operating on mixed structured/unstructured grids.
- **Smoother**: Point/line relaxation, or multi-color Gauss-Seidel.

> ⚠️ Open issue: Ensuring semi-coarsening direction consistency on polyhedral grids where coordinate-alignment directions are undefined. See [[Literature Survey - PMG-Preconditioned BiCGSTAB]] §4, Open Question 1.

## 3-3. Parallelization: Communication Reduction and Switch to Serial Computation

- **Communication reduction**: Applying the node-aware multi-stage messaging [@bienz2020] concept to reduce the number of inter-node messages in halo exchange. Eliminating redundant communication during grid transfer.
- **Serial computation below a fixed level**: From the point where the aggregation level becomes sufficiently small (unknowns per process < threshold) that communication dominates computation, the coarse problem is gathered onto a single (or a few) processes and handled by a **serial direct solver**. The optimal switching level is determined by the performance model of the prior DDM aggregation study [@ha2025aggregation].
- **Alternative**: Comparative evaluation against coarse-grid redistribution (incremental agglomeration) [@reisner2018].

## 3-4. BiCGSTAB / FBiCGSTAB Coupling

- When the MG preconditioner is **fixed** (full V-cycle): apply standard BiCGSTAB.
- When the MG preconditioner is **variable** (iterative approximation of the coarsest grid): apply **FBiCGSTAB** [@chen2016fbicgstab]. The trade-off between inner MG accuracy and outer convergence is determined experimentally.

### Algorithm Overview (FBiCGSTAB + PMG Preconditioner)

```
Given A, b, initial guess x0
r0 = b - A x0,  r̂0 arbitrary (r̂0·r0 ≠ 0)
for k = 0, 1, 2, ...
    y_k = M_k^{-1} p_k        # M_k: PMG preconditioner (may be variable)
    v_k = A y_k
    α = (r̂0·r_k) / (r̂0·v_k)
    s_k = r_k - α v_k
    z_k = M_k^{-1} s_k        # PMG preconditioner reapplied
    t_k = A z_k
    ω = (t_k·s_k) / (t_k·t_k)
    x_{k+1} = x_k + α y_k + ω z_k
    r_{k+1} = s_k - ω t_k
    convergence test: ||r_{k+1}|| / ||b|| < tol
end
```

## 3-5. Verification and Performance Metrics

- **Verification**: Confirm accuracy with the method of manufactured solutions (MMS) and standard benchmarks (lid-driven cavity, channel flow).
- **Performance metrics**: Mesh-independence of iteration count, strong/weak scaling, wall-clock speedup ratio over standalone BiCG, and parallel efficiency.

---

# 4. Case Studies

All three case studies were selected to demonstrate the value of a **multidimensional high-speed pressure solver**. The key common motivations are as follows.

- **Time advancement over long physical time scales**: In all three problems, the phenomena of interest are governed by slow transport/circulation/feedback time scales ranging from seconds to tens of thousands of seconds (a day). Whether reaching a steady state (4-1, 4-2) or a long-duration operational transient (4-3), one must time-advance many steps with small time intervals, solving the pressure correction equation at every step. Therefore, **the acceleration of a single linear solve accumulates multiplicatively over the whole analysis time**, and accelerating the pressure solver directly shortens the total wall-clock.
- **Inevitability of multidimensionality (3D)**: The core physics of the three problems (advection-diffusion of DNP in the fuel salt, thermal stratification and recirculation of buoyancy-driven natural circulation, core flow redistribution during power variations) is inherently coupled to a 3D flow field and cannot be reduced to a 1D system code. Because local velocity/temperature distributions are directly tied to safety margins, multidimensional CFD analysis is essential.
- **Iterative analysis for design and safety evaluation**: Steady-state and operational-transient solutions are computed tens to hundreds of times in design-variable sweeps, sensitivity analysis, and safety licensing evaluation. A high-speed solver lowers this iterative cost to a practical level and shortens the design optimization cycle.

## 4-1. MSFR Steady-State Analysis

**Why high-speed multidimensional steady-state analysis matters**: The MSFR has a unique structure in which the nuclear fuel circulates the primary loop in a liquid state, and as the fuel salt passes through the core in roughly 3–4 seconds, delayed neutron precursors (DNP) are transported outside the active region [@msfr_coupled]. The spatial distribution of this DNP depends strongly on the core geometry and velocity field, so it cannot be captured by point kinetics or a 1D model; only by solving a **coupled 3D steady-state distribution of neutron flux, temperature, and velocity** can the effective reactivity and power distribution be accurately predicted. The steady-state operating point is the reference for fuel-salt composition, flow rate, and core geometry design, and is recomputed at every design iteration, so the acceleration of the pressure correction solver governs the efficiency of design evaluation.

- **Target**: 3000 MWth-class MSFR primary loop, steady turbulence [@msfr_coupled].
- **Physics**: Fuel-salt circulation in 3–4 seconds, asymmetric advection term due to DNP transport → assessment of impact on preconditioner quality.
- **Measurement**: Pressure-correction convergence and total computation time of PMG-BiCGSTAB vs standalone BiCG.

## 4-2. iSMR ECT Natural Circulation

**Why high-speed multidimensional steady-state analysis matters**: The passive decay-heat removal of the iSMR relies on natural circulation driven by buoyancy alone without pumps, so the core of safety demonstration is to show that natural circulation stably converges to a **steady circulation state** with sufficient heat-removal flow. The natural circulation flow rate is a global, multidimensional equilibrium determined by the loop-wide density difference, thermal stratification, and recirculation structure; due to the weak driving force, reaching steady state is slow, and transient analysis requires a very large number of time steps. Moreover, buoyancy-driven low-Reynolds flow converges slowly and unstably, so **without a robust yet fast pressure solver, the steady-state analysis itself becomes impractically slow**. This problem therefore best reveals the practical value of a high-speed, robust solver.

- **Target**: Emergency core cooling / heat-removal (ECT) natural-circulation scenario of an innovative SMR (iSMR).
- **Physics**: Buoyancy-driven low-Reynolds, thermally stratified → strong anisotropy + weak convection. An ideal case for testing semi-coarsening robustness.
- **Measurement**: Speedup ratio to reach natural-circulation steady state, robustness (divergence-free convergence).

## 4-3. iSMR Core Daily Load-Following Flow Distribution Transient Analysis

**The value of a "long-duration transient" beyond a single steady-state solution**: Whereas the previous two cases (4-1, 4-2) address the acceleration of a single steady-state solution, this case analyzes a **long-duration transient on the scale of 24 hours of physical time** by time advancement. This is the configuration in which the multiplicative accumulation effect of solver acceleration is revealed most dramatically: a tiny acceleration of a single pressure solve accumulates over tens of thousands to hundreds of thousands of time steps and determines the total analysis time. Thus, in addition to 4-1 (preconditioner quality) and 4-2 (robustness), it provides a third evaluation axis: **cumulative acceleration and sustained robustness over long-duration time advancement**.

### 4-3-1. Why the "Day (24-hour)" Scale — Physical Justification

1. **Intrinsic long-period nature of load-following operation**: The iSMR is designed as a flexible power source to compensate for the variability of renewable energy, targeting **daily load-following** operation with a soluble-boron-free core and control-rod-based reactivity control. A representative load-following pattern (e.g., a 100%–50%–100% daily cycle) has 24 hours as one cycle, so directly analyzing the core behavior of this operational scenario requires time-advancing the full 24 hours of physical time.
2. **Slow feedback time constants**: The core's key feedback mechanisms have time constants on the scale of hours to a day.
   - **Xenon-135 kinetics**: Reactivity feedback from the I-135 (half-life ≈ 6.6 h) → Xe-135 (half-life ≈ 9.2 h) decay chain develops on a scale of hours, and the period of spatial xenon oscillations is about 15–30 hours. The re-equilibration of the xenon distribution after a power change is captured only on the scale of a day.
   - **Thermal inertia**: Temperature-distribution changes due to the heat capacity of the core, structures, and primary-system coolant proceed on a scale of minutes to hours.
   - These feedbacks cannot be reduced to a single steady-state solution, and the **entire time history** is the subject of safety/operability evaluation.
3. **3D flow redistribution following power variation**: During a power ramp, the core power distribution changes, and accordingly the flow rate and temperature of each subchannel are redistributed. Since the relative contributions of natural-circulation and mixed-convection components vary with time, it cannot be reduced to a 1D system code, and evaluating the time history of local DNBR and temperature margins requires a **transient analysis of the 3D core flow field**.

### 4-3-2. Meaningful Demonstration Plan

| Metric | What it shows | Expected result |
|---|---|---|
| **Cumulative wall-clock curve** | Compare cumulative pressure-solver time versus physical time (0→24 h) for PMG-BiCGSTAB vs standalone BiCG | The gap between the two curves widens linearly over time → an intuitive presentation of multiplicative accumulation ("an analysis that takes N days with BiCG done in M hours with PMG") |
| **Per-time-step iteration-count stability** | Whether the iteration count stays flat, independent of grid and time, even through intervals where the flow field and condition number change due to power ramps and xenon oscillations | PMG stays flat; standalone BiCG shows a sharp rise in iteration count and risk of divergence in the rapid-ramp intervals → robustness demonstration |
| **Time-step-size sensitivity** | Preconditioner robustness when the condition number worsens due to weakened diagonal dominance at large time steps | PMG allows larger time steps → double acceleration of fewer steps × per-step speedup |
| **(Optional) Computational-resource perspective** | Core-hours and power consumption of long-duration HPC analysis | Quantify the real cost savings accumulated over design iterations (sensitivity / licensing) |

- **Scenario definition (example)**: 100%–50%–100% daily load-following (example ramp rate: a few % of rated per minute), or a partial-power step change that induces spatial xenon oscillations. The scenario is designed to sufficiently excite changes in the core flow distribution.
- **Key message**: This case most directly shows the value of a high-speed, robust pressure solver not as a one-off benchmark but in an **actual SMR operational analysis workload**. It is the capstone case that integrally validates steady-state acceleration (4-1) and robustness (4-2) at a real operational time scale.

- **Target**: Daily load-following operational transient of the iSMR core (or a representative partial core), weakly or strongly coupled neutronics/thermal-hydraulics.
- **Physics**: Power ramp + xenon kinetics + thermal inertia → time-varying flow field. Long-duration time advancement.
- **Measurement**: 24-hour cumulative speedup ratio, iteration-count stability over the entire transient, allowable time-step enlargement.

---

# 5. Results (TODO — to be written after experiments)

- [ ] Mesh-independent convergence curve (iteration count vs unknowns)
- [ ] Speedup ratio table over standalone BiCG (MSFR / iSMR)
- [ ] Strong/weak scaling (core count vs efficiency)
- [ ] Communication-reduction effect (communication-time fraction vs level)
- [ ] Serial-switch level sensitivity analysis
- [ ] Fixed MG vs variable MG (FBiCGSTAB) comparison
- [ ] (4-3) Cumulative wall-clock curve (physical time 0→24 h vs cumulative solver time)
- [ ] (4-3) Per-time-step iteration-count stability over the entire transient (ramp / xenon-oscillation intervals)
- [ ] (4-3) Time-step-size sensitivity (large steps allowed by PMG robustness)

---

# 6. Conclusion (TODO)

We integrate a semi-coarsening geometric MG applicable to unstructured grids as a BiCGSTAB preconditioner, secure parallel scalability through communication reduction and serial switching, and demonstrate acceleration on the MSFR and iSMR problems (planned).

---

# 7. References (BibTeX)

> Note: Because publisher pages (ScienceDirect/Springer/SSRN) were blocked, the **author order, pages, article number, and DOI** of some entries need cross-checking against the published versions (marked ⚠️ below). arXiv/DOI links are hyperlinked.

## 7-1. Authors' Group Prior Work (Scalable Solver Lineage)

```bibtex
@article{do2024meshless,
  title   = {Highly scalable meshless multigrid solver for 3D thermal-hydraulic analysis of nuclear reactors},
  author  = {Do, Seong Ju and Ha, Sang Truong and Choi, Hyoung Gwon and Yoon, Han Young},
  journal = {Annals of Nuclear Energy},
  volume  = {207},
  pages   = {110713},
  year    = {2024},
  doi     = {10.1016/j.anucene.2024.110713},
  note    = {\url{https://www.sciencedirect.com/science/article/pii/S0306454924003761}}
}
% ⚠️ Verify author order and article number (110713) against published version

@article{do2022highaspect,
  title   = {A meshless geometric multigrid method for a grid with a high aspect ratio},
  author  = {Do, Seong Ju and Ha, Sang Truong and Choi, Hyoung Gwon and Yoon, Han Young},
  journal = {Journal of Mechanical Science and Technology},
  year    = {2022},
  doi     = {10.1007/s12206-022-1019-4},
  note    = {\url{https://link.springer.com/article/10.1007/s12206-022-1019-4}}
}
% ⚠️ Verify authors, volume, pages

@article{ha2024aiaa,
  title   = {Meshless Geometric Multigrid Method for Complex Geometries with Improved Cell Coarsening Algorithm},
  author  = {Ha, Sang Truong and Do, Seong Ju and Choi, Hyoung Gwon and Yoon, Han Young},
  journal = {AIAA Journal},
  year    = {2024},
  doi     = {10.2514/1.J063127},
  note    = {\url{https://arc.aiaa.org/doi/10.2514/1.J063127}}
}
% ⚠️ Verify author order, volume, year

@article{ha2021nodecoarsening,
  title   = {A meshless geometric multigrid method based on a node-coarsening algorithm for the linear finite element discretization},
  author  = {Ha, Sang Truong and Choi, Hyoung Gwon},
  year    = {2021},
  note    = {\url{https://www.researchgate.net/publication/351771284}}
}
% ⚠️ Verify authors, journal, DOI (ResearchGate-registered version)

@misc{ha2025aggregation,
  title        = {Investigation on an Optimal Aggregation Level for a Parallel Meshless Multigrid Method Based on Domain Decomposition Method},
  author       = {Ha, Sang Truong and Yoon, Han Young and Choi, Hyoung Gwon},
  howpublished = {SSRN preprint 5081886},
  year         = {2025},
  note         = {\url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5081886}}
}
% ⚠️ Verify year and publication status
```

## 7-2. Semi-Coarsening Multigrid

```bibtex
@article{mulder1989,
  title   = {A new multigrid approach to convection problems},
  author  = {Mulder, W. A.},
  journal = {Journal of Computational Physics},
  volume  = {83},
  number  = {2},
  pages   = {303--323},
  year    = {1989},
  doi     = {10.1016/0021-9991(89)90121-6}
}
% ⚠️ Foundational semi-coarsening reference. Verify exact source (JCP 1989 vs SIAM J. Numer. Anal.)

@techreport{osti440722,
  title       = {A parallel semicoarsening multigrid algorithm for solving the Reynolds-averaged Navier-Stokes equations},
  institution = {OSTI},
  number      = {440722},
  note        = {\url{https://www.osti.gov/biblio/440722}}
}
% ⚠️ Verify authors and year

@article{dendy2019msg,
  title   = {Multiple Semicoarsened Multigrid for anisotropic diffusion problems},
  author  = {Anonymous},
  journal = {arXiv preprint arXiv:1907.12334},
  year    = {2019},
  note    = {\url{https://arxiv.org/abs/1907.12334}}
}
% ⚠️ Verify authors and formal publication venue
```

## 7-3. Parallel MG Communication Reduction

```bibtex
@article{reisner2018,
  title   = {Scaling Structured Multigrid to 500K+ Cores through Coarse-Grid Redistribution},
  author  = {Reisner, Andrew and Olson, Luke N. and Moulton, J. David},
  journal = {SIAM Journal on Scientific Computing},
  year    = {2018},
  note    = {\url{https://arxiv.org/abs/1803.02481}}
}
% ⚠️ Verify volume, pages, DOI

@article{bienz2020,
  title   = {Node-Aware Improvements to Allreduce / Reducing Communication in Algebraic Multigrid},
  author  = {Bienz, Amanda and Gropp, William D. and Olson, Luke N.},
  journal = {International Journal of High Performance Computing Applications},
  volume  = {34},
  number  = {5},
  pages   = {547--561},
  year    = {2020},
  note    = {\url{https://arxiv.org/abs/1904.05838}}
}
% ⚠️ Verify exact title

@article{manteuffel2016,
  title   = {Nonsymmetric Algebraic Multigrid Sparsification / Reducing Communication Costs in AMG},
  author  = {Manteuffel, Thomas A. and others},
  journal = {SIAM Journal on Scientific Computing},
  year    = {2016},
  note    = {\url{https://arxiv.org/abs/1512.04629}}
}
% ⚠️ Verify exact title, full author list, volume, pages
```

## 7-4. BiCGSTAB / Preconditioned Convergence Theory

```bibtex
@article{vandervorst1992,
  title   = {Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems},
  author  = {van der Vorst, Henk A.},
  journal = {SIAM Journal on Scientific and Statistical Computing},
  volume  = {13},
  number  = {2},
  pages   = {631--644},
  year    = {1992},
  doi     = {10.1137/0913035},
  note    = {\url{https://epubs.siam.org/doi/10.1137/0913035}}
}

@article{chen2016fbicgstab,
  title   = {A Flexible Variant of Bi-CGSTAB with a Multigrid Preconditioner (FBiCGStab)},
  author  = {Chen, Jie and McInnes, Lois Curfman and Zhang, Hong},
  journal = {Journal of Scientific Computing},
  volume  = {68},
  number  = {2},
  pages   = {803--825},
  year    = {2016},
  note    = {\url{https://jiechenjiechen.github.io/pub/fbcgs.pdf}}
}
% ⚠️ Verify exact title (PFLOTRAN application)

@article{lucerolorca2026,
  title   = {Two-point spectral clustering and Krylov convergence of multigrid-preconditioned solvers},
  author  = {Lucero Lorca, Jose Pablo and McCoid, Conor and Outrata, Michal},
  journal = {arXiv preprint arXiv:2511.12298},
  year    = {2026},
  note    = {\url{https://arxiv.org/abs/2511.12298}; preprint, pre-peer-review}
}
% ⚠️ Verify title and authors (2026 preprint)
```

## 7-5. MSFR / Reactor CFD

```bibtex
@article{msfr_coupled,
  title   = {Coupled neutronics and thermal-hydraulics numerical simulations of a Molten Salt Fast Reactor (MSFR)},
  author  = {Anonymous},
  note    = {\url{https://www.academia.edu/70734179/}}
}
% ⚠️ Verify authors, journal, year

@article{fnuen2025,
  title   = {Multigrid preconditioning for reactor thermal-hydraulics CFD},
  journal = {Frontiers in Nuclear Engineering},
  year    = {2025},
  note    = {\url{https://www.frontiersin.org/journals/nuclear-engineering/articles/10.3389/fnuen.2025.1597165/full}}
}
% ⚠️ Verify title and authors

@inproceedings{nuclearcfd_springer,
  title     = {Parallel multigrid in nuclear CFD},
  booktitle = {Springer LNCS},
  doi       = {10.1007/978-3-030-39647-3_20},
  note      = {\url{https://link.springer.com/chapter/10.1007/978-3-030-39647-3_20}}
}
% ⚠️ Verify title, authors, year
```

---

# Appendix: Writing Notes

- This draft is based on the verified findings of [[Literature Survey - PMG-Preconditioned BiCGSTAB]]. The 13 rejected claims (exaggeration/misattribution) are not cited.
- The ⚠️ BibTeX entries are unconfirmed due to blocked publisher-page access — cross-check author order, volume, pages, and DOI directly against the published versions before finalizing.
- The Results (§5) and Conclusion (§6) are to be written after the actual numerical experiments.
