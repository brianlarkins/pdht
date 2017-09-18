SOURCES=abstract.tex conclusion.tex\
        implementation.tex introduction.tex\
        overview.tex paper.tex\
        progmodel.tex background.tex results.tex 

all: paper.pdf

paper.pdf:  $(SOURCES)
	latexmk -pdf paper.tex

clean:
	rm -f *~ *.dvi *.pdf *.aux *.bbl *.log *.blg *.synctex.gz *.fdb_latexmk *.fls paper.out *.vim paper_aux
