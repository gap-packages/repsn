
GAP=gap4ws

gapdoc:
        @echo 'RequirePackage("gapdoc"); MakeGAPDocDoc( "./", "repsn", [], "repsnxml" );;' | $(GAP) -q -T

clean:
        @rm -f repsn.aux repsn.bbl repsn.blg repsn.idx repsn.ilg repsn.ind\
           repsn.log repsn.pnr repsn.tex repsn.toc repsn.brf
	@rm -f *~
        

spotless: clean
	@rm -f manual.dvi manual.html manual.pdf manual.ps manual.six
	@rm -f chap*.html chap*.txt
