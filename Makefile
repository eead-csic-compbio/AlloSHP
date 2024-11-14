
install: install_Red install_Red2Ensembl install_Cgaln install_GSAlign

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

install_GSAlign:
	if [ ! -d "lib/GSAlign" ]; then \
		cd lib && git clone https://github.com/hsinnan75/GSAlign.git && cd GSAlign && make && ls bin && rm -rf test; \
	fi

test_install:
	perl WGA -c

test:
	./WGA -A sample_data/Bdis.fna.gz -B sample_data/Bsta.fna.gz && \
		./WGA -A sample_data/Bdis.fna.gz -B sample_data/Bsta.fna.gz -g && \
                ./vcf2alignment -v sample_data/BdisBd2_BstaChr01.vcf.gz -c sample_data/config.tsv -l BdisBd2_BstaChr01.vcf.log -d 5 -m 3 &&\
                ./vcf2synteny -v sample_data/BdisBd2_BstaChr01.vcf.gz -c sample_data/config.synteny.tsv -l BdisBd2_BstaChr01.vcf.log \
                        -d 5 -m 3 -r Bdis -o BdisBd2_BstaChr01.DP5.M3.synteny.fasta; \

# Cgaln only, all logs saved
demo: 
	./WGA -A sample_data/Bdis.fna.gz -B sample_data/Bsta.fna.gz > BdisBd2_BstaChr01.log 2>&1 && \
		./vcf2alignment -v sample_data/BdisBd2_BstaChr01.vcf.gz -c sample_data/config.tsv -l BdisBd2_BstaChr01.vcf.log -d 5 -m 3 \
			> BdisBd2_BstaChr01.DP5.M3.log 2>&1 &&\
		./vcf2synteny -v sample_data/BdisBd2_BstaChr01.vcf.gz -c sample_data/config.synteny.tsv -l BdisBd2_BstaChr01.vcf.log \
			-d 5 -m 3 -r Bdis -o BdisBd2_BstaChr01.DP5.M3.synteny.fasta > BdisBd2_BstaChr01.DP5.M3.synteny.log 2>&1; \

clean:
	rm -rf Bdis* _Cgalnidx		



#cannot cope with raw barley chromosomes, see https://doi.org/10.1186/s13059-023-03071-z
#minimap2release = 2.24
#install_minimap2:
#	if [ ! -d "lib/minimap2" ]; then \
#		cd lib && wget https://github.com/lh3/minimap2/releases/download/v${minimap2release}/minimap2-${minimap2release}.tar.bz2 && \
#			tar xfj minimap2-${minimap2release}.tar.bz2 && cd minimap2-${minimap2release} && make && cd .. && \
#			rm -f minimap2-${minimap2release}.tar.bz2 && ln -fs minimap2-${minimap2release} minimap2; \
#	fi
