
SHELL=/bin/sh
PAPER=suppinfo

.SUFFIXES:
.SUFFIXES: .bib .pdf .tex
.PHONY: clean
.PHONY: keepPDF

run: $(PAPER).pdf
	rm -rf *.aux *.lof *.log *.lot *.toc *.bbl *.blg
$(PAPER).pdf: $(PAPER).bbl $(PAPER).tex
	pdflatex $(PAPER).tex -draftmode
	pdflatex $(PAPER).tex 

$(PAPER).bbl: $(PAPER).aux
	bibtex $(PAPER).aux

$(PAPER).aux: $(PAPER).bib 
	pdflatex $(PAPER).tex -draftmode
	pdflatex $(PAPER).tex -draftmode

clean:
	rm -rf *.aux *.lof *.log *.lot *.toc *.bbl *.blg *pdf
keepPDF:
	rm -rf *.aux *.lof *.log *.lot *.toc *.bbl *.blg
