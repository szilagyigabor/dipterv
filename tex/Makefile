full: tartalom/package.tex tartalom/preamble.tex tartalom/xmp.tex mybib.bib dipterv.tex
	pdflatex -interaction=nonstopmode -halt-on-error dipterv.tex > output.txt
	-bibtex dipterv >> output.txt
	pdflatex -interaction=nonstopmode -halt-on-error dipterv.tex >> output.txt
	pdflatex -interaction=nonstopmode -halt-on-error dipterv.tex >> output.txt
	rm -f dipterv.aux dipterv.bbl dipterv.blg dipterv.out dipterv.log dipterv.toc dipterv.xmpdata pdfa.xmpi tartalom/*.aux texput.log
	-mv dipterv.pdf ../pdf
