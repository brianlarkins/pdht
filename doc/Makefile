SOURCES=abstract.tex conclusion.tex\
        implementation.tex introduction.tex\
        overview.tex paper.tex\
        progmodel.tex related.tex results.tex 

all: paper.pdf

paper.pdf:  $(SOURCES)
	pdflatex paper
	pdflatex paper
	pdflatex paper

clean:
	rm -f *~ *.dvi *.pdf *.aux *.bbl *.log *.blg *.synctex.gz *.fdb_latexmk *.fls
