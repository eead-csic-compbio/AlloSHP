
minimap2release = 2.24

install: install_Red install_Red2Ensembl install_Cgaln install_gnuplot
	#-sudo apt install -y git wget

install_Red:
	if [ ! -d "lib/Red" ]; then \
		cd lib && git clone https://github.com/EnsemblGenomes/Red.git && cd Red/src_2.0 && make bin && make && \
			cd .. && rm -f bin/*.o && rm -f bin/*/*.o; \
	fi

install_Red2Ensembl:
	if [ ! -e "utils/Red2Ensembl.py" ]; then \
		cd utils && wget https://raw.githubusercontent.com/Ensembl/plant-scripts/refs/heads/master/repeats/Red2Ensembl.py && \
			chmod +x Red2Ensembl.py; \
	fi

install_Cgaln:
	if [ ! -d "lib/Cgaln" ]; then \
		cd lib && git clone https://github.com/rnakato/Cgaln.git && cd Cgaln && make && rm -f *.fasta *.o; \
	fi

install_gnuplot:
	if [ ! -d "lib/gnuplot-5.2.8" ]; then \
		cd lib && wget https://unlimited.dl.sourceforge.net/project/gnuplot/gnuplot/5.2.8/gnuplot-5.2.8.tar.gz && \
			tar xfz gnuplot-5.2.8.tar.gz && cd gnuplot-5.2.8 && ./configure with_latex=no && make && \
			rm -rf demo docs && rm -f src/*.o && cd .. && rm -f gnuplot-5.2.8.tar.gz; \
        fi

# optional, tests
install_minimap2:
	if [ ! -d "lib/minimap2" ]; then \
		cd lib && wget https://github.com/lh3/minimap2/releases/download/v${minimap2release}/minimap2-${minimap2release}.tar.bz2 && \
			tar xfj minimap2-${minimap2release}.tar.bz2 && cd minimap2-${minimap2release} && make && cd .. && \
			rm -f minimap2-${minimap2release}.tar.bz2 && ln -fs minimap2-${minimap2release} minimap2; \
	fi
