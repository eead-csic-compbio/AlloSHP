
minimap2release = 2.24

install: install_minimap2 install_red install_Cgaln
	-sudo apt install -y wget 

install_minimap2:
	if [ ! -d "bin/minimap2" ]; then \
		cd lib && wget https://github.com/lh3/minimap2/releases/download/v${minimap2release}/minimap2-${minimap2release}.tar.bz2 && \
			tar xfj minimap2-${minimap2release}.tar.bz2 && cd minimap2-${minimap2release} && make && cd .. && \
			rm -f minimap2-${minimap2release}.tar.bz2 && ln -fs minimap2-${minimap2release} minimap2; \
	fi

install_Red:
	cd bin && git clone https://github.com/EnsemblGenomes/Red.git && cd Red/src_2.0 && make bin && make
	#in case you need to use an alternative g++ compiler
        #cd bin && git clone https://github.com/EnsemblGenomes/Red.git && cd Red/src_2.0 && make bin && make CXX=g++-10

install_Cgaln:
        cd bin && git https://github.com/rnakato/Cgaln.git && cd Cgaln && make
