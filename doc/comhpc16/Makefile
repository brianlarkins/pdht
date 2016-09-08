SOURCES = $(shell find . -name '*.tex' -print)
BIBS    = $(shell find . -name '*.bib' -print)

TEX = paper.tex

CLS = IEEEtran.cls

STY = 

PACKAGE = 

PDF = ${TEX:%.tex=%.pdf}

PWD= ${shell pwd}

all: $(PDF)

%.pdf: %.tex $(CLS) $(STY) $(SOURCES) $(BIBS)
	pdflatex $*
	bibtex $*
	pdflatex $*
	while ( grep -q '^LaTeX Warning: Label(s) may have changed' $*.log)\
	do pdflatex $*; done

clean:
	$(RM) *.aux *.log  *.blg

distclean: clean
	$(RM) $(PDF) *.bbl
