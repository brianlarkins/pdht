TEX = paper.tex

CLS = acmsmall.cls

STY = acmcopyright.sty

PACKAGE = acmsmall

PDF = ${TEX:%.tex=%.pdf}

PWD= ${shell pwd}

all: $(PDF)

%.pdf: %.tex $(CLS) $(STY)
	pdflatex $*
	bibtex $*
	pdflatex $*
	while ( grep -q '^LaTeX Warning: Label(s) may have changed' $*.log)\
	do pdflatex $*; done

clean:
	$(RM) *.aux *.log  *.blg

distclean: clean
	$(RM) $(PDF) acm-update.pdf *.bbl
