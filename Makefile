
minimap2release = 2.24

install: install_Red install_Cgaln

install_Red:
	if [ ! -d "bin/Red" ]; then \
		cd bin && git clone https://github.com/EnsemblGenomes/Red.git && cd Red/src_2.0 && make bin && make; \
	fi

install_Cgaln:
	if [ ! -d "bin/Cgaln" ]; then \
		cd bin && git clone https://github.com/rnakato/Cgaln.git && cd Cgaln && make && rm -f *.fasta; \
	fi


# optional, tests
install_minimap2:
	if [ ! -d "bin/minimap2" ]; then \
		cd bin && wget https://github.com/lh3/minimap2/releases/download/v${minimap2release}/minimap2-${minimap2release}.tar.bz2 && \
			tar xfj minimap2-${minimap2release}.tar.bz2 && cd minimap2-${minimap2release} && make && cd .. && \
			rm -f minimap2-${minimap2release}.tar.bz2 && ln -fs minimap2-${minimap2release} minimap2; \
	fi
