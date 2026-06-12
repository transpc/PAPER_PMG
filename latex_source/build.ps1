# LaTeX Workshop 의 기본 latexmk 레시피와 동일하게 main.tex 를 빌드한다.
# pdflatex -> bibtex -> pdflatex -> pdflatex 를 latexmk 가 알아서 반복 실행한다.

# 스크립트가 있는 디렉터리로 이동
Set-Location -Path $PSScriptRoot

latexmk -pdf -synctex=1 -interaction=nonstopmode -file-line-error main.tex
