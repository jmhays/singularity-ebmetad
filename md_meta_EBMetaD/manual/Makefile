manual.pdf : manual.ind manual.bbl
	pdflatex manual
manual.ind : manual.idx manual.toc
	latex manual
	makeindex manual
manual.idx : *.tex
	latex manual
manual.bbl : *.tex
	bibtex manual
	latex manual
manual.toc : *.tex
	latex manual
clean :
	rm  manual.dvi manual.aux  manual.idx  manual.ilg  manual.ind  manual.log  manual.out  manual.pdf manual.toc  manual.bbl  manual.blg
