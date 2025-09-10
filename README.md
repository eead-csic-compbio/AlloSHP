# AlloSHP

A command-line tool for detecting and extracting single homeologous polymorphisms (SHPs) from the subgenomes of allopolyploid species. This tool integrates three main algorithms, WGA, VCF2ALIGNMENT and VCF2SYNTENY, and allows the detection of SHPs for the study of diploid-polyploid complexes with available diploid progenitor genomes. It requires FASTA of the reference genomes and [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) input files and produces multiple sequence alignments (MSA) of the SHP for each allopolyploid subgenome.

Rubén Sancho (1,2), Pilar Catalán (2), Bruno Contreras Moreira (1,3)

1. Estación Experimental de Aula Dei-CSIC, Zaragoza, Spain
2. Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
3. Fundación ARAID, Zaragoza, Spain

https://anaconda.org/bioconda/alloshp/badges/version.svg


[![Build Status](https://github.com/eead-csic-compbio/AlloSHP/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/eead-csic-compbio/AlloSHP/actions/workflows/ci.yml)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/alloshp/badges/version.svg)](https://anaconda.org/bioconda/alloshp)

## Installation

### Bioconda

A bioconda recipe is available at <https://anaconda.org/bioconda/alloshp> for linux-64 and osx-64 systems.
After activating a conda environment this software can be installed as follows:

    conda install bioconda::alloshp

Note the recipe includes sample data, which you can analyze in a few minutes with:

    make -f ${CONDA_PREFIX}/Makefile test_conda

### Local compilation

This protocol has been developed and tested on Linux x86_64 systems, but it might work on MacOS settings with tweaks.
It requires some standard Linux utilities (gzip, sort, perl, make, python3, g++, etc) and a few third-party dependencies
which can be installed locally as follows:

    git clone https://github.com/eead-csic-compbio/AlloSHP.git
    cd AlloSHP
    make install
    # optionally, takes a few minutes
    make test

### Dependencies 

The table shows the main dependencies of this package:

|software|flag|source|notes|
|:-------|:---|:-----|:----|
|Cgaln|  |https://github.com/rnakato/Cgaln|requires gcc compiler|
|GSAlign| -g |https://github.com/hsinnan75/GSAlign|requires g++ compiler|
|Red|skipped with -m|https://github.com/EnsemblGenomes/Red|requires g++ compiler|
|Red2Ensembl.py|   |https://github.com/Ensembl/plant-scripts|requires python3|
|gnuplot|  |http://www.gnuplot.info|required for dotplots in PDF format| 

Other Linux dependencies required include: `wget python3 g++ gnuplot-qt libdb-dev` 

See [ci.yml](https://github.com/eead-csic-compbio/vcf2alignment/blob/master/.github/workflows/ci.yml).

## Pipeline overview

<!-- flowchart -->

## 1) Input data 

These are the data files required to run this pipeline:

+ 2+ FASTA files of genome assemblies of diploid species of the taxa of interest, one file per species.
Our *Brachypodium* [sample data](https://github.com/eead-csic-compbio/vcf2alignment/tree/master/sample_data)
includes the GZIP compressed genomes of *Brachypodium distachyon* and *Brachypodium stacei*, 
which are named `Bdis.fna.gz` and `Bsta.fna.gz`. **Each file should have unique chromosome names**.

+ A VCF file with GBS/RADseq/RNAseq/WGS sequence reads of presumably polyploids mapped to the concatenated previously described FASTA files.
See file `BdisBd2_BstaChr01.vcf.gz` in [sample data](https://github.com/eead-csic-compbio/vcf2alignment/tree/master/sample_data)
for an example.

+ Configuration files (see [sample data](https://github.com/eead-csic-compbio/vcf2alignment/tree/master/sample_data)) for `vcf2alignment` (config.tsv) and `vcf2synteny` scripts (config.synteny.tsv (without outgroup) or config.outg.synteny.tsv (with outgroup)).


## 2) How to run

A typical analysis involves calling three scripts: [WGA](https://github.com/eead-csic-compbio/AlloSHP/blob/master/WGA), 
[vcf2alignment](https://github.com/eead-csic-compbio/AlloSHP/blob/master/WGA) and 
[vcf2synteny](https://github.com/eead-csic-compbio/AlloSHP/blob/master/WGA). 

Note that the steps below can be run in one go as follows:

    make test

### 2.1) Whole-genome alignments (WGA)

Script `WGA` must be computed to find syntenic segments among the diploid reference genomes available for read mapping.
By default, this uses [CGaln](https://github.com/rnakato/Cgaln),
which requires the input sequences to be [soft-masked](https://genomevolution.org/wiki/index.php/Masked) 
(this is also taken care of by the script).
The table shows the flags of this script:

|flag|note|
|:-------|:---|
|-h  | this message|
|-A  | FASTA file of genome A (example: -A speciesA.fna[.gz])|
|-B  | FASTA file of genome B example: -B speciesB.fna[.gz])|
|-o  | output folder (optional, default: speciesA.speciesB)|
|-l  | min contig length [Mbp] (optional, default: -l $minlenMb)|
|-m  | FASTA files already soft-masked (optional, default: masked with Red)|
|-n  | number of cores  (optional, some tasks only, default: $ncores)|
|-g  | use multithreaded GSAlign algorithm (optional, default: Cgaln)|
|-C  | parameters for Cgaln aligner (optional, default: -C '-X4000'), where X stands for X-drop-off threshold for gapped extension of HSP ([paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-224))|
|-N  | parameters for Cgaln indexer (optional, default: -N '-K11 -BS10000'), where K is k-mer size and BS block size ([paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-224))|
|-G  | parameters for GSAlign aligner (optional, default: -G '-no_vcf -one'), see [paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6569-1)|
|-M  | parameters for utils/mapcoords.pl (optional, default: -M '0.25 0.05'. 1st: max ratio of mapped positions in other blocks; 2nd: max ratio of coordinates with multiple positions in the same block)|
|-t  | path to dir for temp files (optional, default: -t /tmp)|
|-c  | print credits and checks install (recommended)|

In our example, we set *B. distachyon* as the master reference genome for being the best quality assembly at hand
(see section 2.3 below). We now find out syntenic segments on the other genomes, defined as secondary genomes (*B. stacei*). 
Each diploid reference is hence considered a **subgenome** to which reads map:

    ./WGA -A sample_data/Bdis.fna.gz -B sample_data/Bsta.fna.gz

This produces the following main output files (see below and square):

+ BED: Bdis.fna.gz.Bsta.fna.gz/Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X4000_0.25_0.05.bed
+ LOG: Bdis.fna.gz.Bsta.fna.gz/Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X4000_0.25_0.05.coords.log
+ PDF: Bdis.fna.gz.Bsta.fna.gz/Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X4000_0.25_0.05.dot.pdf

And this standard out in the terminal:

    ## ./WGA -A sample_data/Bdis.fna.gz -B sample_data/Bsta.fna.gz -o Bdis.fna.gz.Bsta.fna.gz -l 1 -m 1 -G 0 -N -K11 -BS10000 -C -X4000 -M 0.25 0.05 -n 4

    ## root: Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X4000_0.25_0.05

    # filter_FASTA_sequences: [passed 59130575 bp] Bd2
    # filter_FASTA_sequences: [passed 31564145 bp] Chr01

    # chrcode A example: Bd(\d+)
    # chrcode B example: Chr(\d+)

    ## soft-masking filtered sequences

    # filtered file A: sample_data/Bdis.fna.gz.sm.fasta
    # filtered file B: sample_data/Bsta.fna.gz.sm.fasta

    ## CGaln algorithm


    ## indexing masked, filtered sequences

    [...]

    ## computing and plotting Whole Genome Alignment

    [...]

    ## converting Whole Genome Alignment to BED

    ## output files:

    # BED: Bdis.fna.gz.Bsta.fna.gz/Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X4000_0.25_0.05.bed
    # LOG: Bdis.fna.gz.Bsta.fna.gz/Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X4000_0.25_0.05.coords.log
    # PDF: Bdis.fna.gz.Bsta.fna.gz/Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X4000_0.25_0.05.dot.pdf

    ## WGA summary: valid blocks: 54 unique positions: 8956568

    # time used (s): 119 memory used (Mb): 407.9

The most important result files are the 0-based **BED** list of syntenic positions, which will be used in the last step,
and the **PDF** dotplot, which requires `gnuplot` in your system, which must be inspected to assess the quality of
the WGA. Note that flags `-N` and `-C` can be used to tweak the WGA parameters after inspection of the dotplot.
Also, you can use `-o` to set your own **output folder**.

The BED file produced by `WGA` looks like this:

    # chrA  posA    endA    baseA   strandA chrB    posB    endB    baseB   block   SNP
    Bd2     44989686        44989687        A       -       Chr01   13393558        13393559        A       10      .
    Bd2     44989685        44989686        T       -       Chr01   13393559        13393560        T       10      .
    Bd2     44989684        44989685        C       -       Chr01   13393560        13393561        C       10      .
    Bd2     44989683        44989684        C       -       Chr01   13393561        13393562        C       10      .
    Bd2     44989682        44989683        T       -       Chr01   13393562        13393563        T       10      .
    Bd2     44989681        44989682        G       -       Chr01   13393563        13393564        G       10      .
    Bd2     44989680        44989681        C       -       Chr01   13393564        13393565        C       10      .
    Bd2     44989679        44989680        T       -       Chr01   13393565        13393566        T       10      .
    Bd2     44989678        44989679        G       -       Chr01   13393566        13393567        A       10      SNP
    Bd2     44989677        44989678        C       -       Chr01   13393567        13393568        C       10      .
    ...

    Where:
    + chrA: Name of the chromosome regarding the Master reference genome
    + posA: Start position regarding the Master reference genome
    + endA: End position regarding the Master reference genome
    + baseA: Nucleotide base regarding the Master reference genome
    + strandA: Defined as + (forward) or - (reverse) regarding the Master reference genome
    + chrB: Name of the chromosome regarding the Secondary reference genome
    + posB: Start position regarding the Secondary reference genome
    + endB: End position regarding the Secondary reference genome
    + baseB: Nucleotide base regarding the Secondary reference genome
    + block: Syntenic block according to CGaln
    + SNP: SNP between reference genomes
 

Alternatively, the GSAlign WGA algorithm can be invoked as follows, with flag `-g`:
 
    ./WGA -A sample_data/Bdis.fna.gz -B sample_data/Bsta.fna.gz -g

Note that you can change the default GSAlign settings with the optional flag `-G`. 
This might be required for your genomes of interest.

![whole-genome alignment plot](./pics/Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X12000_-fc_-cons.dot.png)

*Figure 1. WGA dotplot resulting from Cgaln alignment of two homologous chromosomes.*

![whole-genome alignment plot](./pics/dotplot.png)

*Figure 2. Multi-chromosome WGA dotplot.* A good WGA dotplot will have long diagonal runs of aligned genome regions instead of clouds.

### 2.2) Filtering valid positions in the VCF file

Script `vcf2alignment` parses an input VCF file, which might be GZIP/BZIP2 compressed, 
and produces a list of valid sites considering min read depth (-d) and max missing data (-m).
Note this requires a config file that matches sample names in the VCF file to their human-readable names
(see example [sample_data/config.tsv](https://github.com/eead-csic-compbio/vcf2alignment/blob/master/sample_data/config.tsv)).
A report log file with valid 1-based coordinates (-l) is saved to be used in the last step.

The table shows the flags of `vcf2alignment`:

|flag|note|
|:-------|:---|
|-h  | this message|
|-v  | input VCF file (example: -v data.vcf.gz)|
|-c  | input TSV config file (example: -c config.tsv)|
|-l  | output report file name, 1-based coordinates (example: -l vcf.report.log.gz)|
|-o  | output MSA file name (optional, example: -o out.fasta)|
|-d  | min read depth at each position for each sample (optional, example: -d 3, default -d 3, use -d 0 if VCF file lacks DP)|
|-m  | max missing samples (optional, example: -m 10, default -m 10)|
|-f  | output format (optional, example: -f nexus, default -f fasta)|
|-p  | take only polymorphic sites (optional, by default all sites, constant and SNPs, are taken)|
|-H  | take also heterozygous sites (optional, by default only homozygous are taken)|

The configuration file structure (TSV format) for `vcf2alignment`:

    # original_sample_header	final_sample_header	config_tag
	
    sample1.bam	sample1	real_name
    sample2.bam	sample2	real_name
    ...
    sampleN.bam	sampleN	real_name
	
    # original_sample_header (1st column): For each sample, the sample name as shown in the input VCF file.
    # final_sample_header (2nd column): For each sample, the user-chosen sample name to be displayed in downstream results.
    # config_tag (3rd column): For each sample, a mandatory tag (real_name) that must appear for `vcf2alignment` correct processing.

In our example, we use a toy VCF file that contains read-mapping positions on chromosomes Bd2 of *Brachypodium distachyon* and Chr01 of *Brachypodium stacei*:

    ./vcf2alignment -v sample_data/BdisBd2_BstaChr01.vcf.gz -c sample_data/config.tsv -l BdisBd2_BstaChr01.vcf.log.gz -d 5 -m 3

This produces the following main output file:

+ LOG: BdisBd2_BstaChr01.vcf.log.gz

Its contents look like this:

    # number of samples found=6
    # valid locus: Bd2_15643 3 A
    # valid locus: Bd2_15648 3 C
    # valid locus: Bd2_16138 2 T
    # ...
    # valid locus: Chr01_8301300 3 A
    # valid locus: Chr01_8301301 3 G
    # number of valid loci=1746620
    # number of polymorphic loci=22248

    # stats (#SNPs):
      ABR2.bam : 971194
      ...
      Bhyb26.bam : 1704849
      ...

    # stats per contig/chr (#SNPs):
      ABR2.bam        Bd2     962203
      ABR2.bam        Chr01   8991
      ...
      Bhyb26.bam      Bd2     945905
      Bhyb26.bam      Chr01   758944
      ...
    # stats (#N):
      ABR2.bam : 775426
      ...
      Bhyb26.bam : 41771
      ...

    # stats per contig/chr (#N):
      ABR2.bam        Bd2     3670
      ABR2.bam        Chr01   771756
      ...
      Bhyb26.bam      Bd2     19968
      Bhyb26.bam      Chr01   21803
      ...
    


### 2.3) Producing a multiple sequence alignment of polyploid subgenomes

Script `vcf2synteny` parses the VCF file obtained from reads mapped to multiple concatenated reference genomes, the LOG file of valid loci computed by `vcf2alignment`, and the synteny-based equivalent coordinates (BED file) computed by `WGA` to align the polymorphic loci referenced by syntenic positions, separating them on each reference genome, and defining them as SHPs. The resulting MSA will have as many subgenomes as reference genomes were used.

The table shows the flags of `vcf2synteny`:

|flag|note|
|:-------|:---|
|-h  |  this message|
|-v  |  input VCF file (example: -v data.vcf.gz)|
|-l  |  report file from vcf2alignment, 1-based (example: -l vcf.rport.log.gz)|
|-c  |  input TSV config file (example: -c config.synteny.tsv)|
|-o  |  output FASTA file name (example: -o out.fasta)|
|-r  |  master reference genome (example: -r Bdis)|
|-d  |  min depth of called SNPs (optional, example: -d 3, default -d 3)|
|-m  |  max missing samples (optional, example: -m 10, default -m 10)|
|-V  |  output VCF file name (optional, coordinates from -r genome, example: -f out.vcf)|
|-1  |  syntenic coords are 1-based (optional, 0-based/BED by default, WGA in -c config.synteny.tsv)|
|-p  |  take only polymorphic sites (optional, by default all sites, constant and SNPs, are taken)|
|-H  |  take also heterozygous sites (optional, by default only homozygous, requires -l report from vcf2alignment -H)|
|-N  |  new temp files, don't re-use (optional, by default temp files are re-used if available at -t)|
|-t  |  path to dir for temp file  (optional, default: -t $tmppath)|

The configuration file structure (TSV format) for `vcf2synteny` includes four blocks:

    # -----------------------------------------
	# ---------- BLOCK A (mandatory) ----------
    # -----------------------------------------
	
    # original_sample_header	final_sample_header	config_tag
	
    sample1.bam	sample1	real_name
    sample2.bam	sample2	real_name
    ...
	# -----------------------------------------
    # ---------- BLOCK B (mandatory) ----------
	# -----------------------------------------
 
	# ID_ref	BED_file	config_synteny_tag
 
    Ref_2	PATH/Ref_master.Ref_2.bed	WGA
	Ref_3	PATH/Ref_master.Ref_3.bed	WGA
    ...
	# -----------------------------------------
    # ---------- BLOCK C (mandatory) ----------
	# -----------------------------------------
 
	# ID_ref	Chr_code	config_synteny_tag
 
    Ref_master	Chr_code_master	chrcode
    Ref_2	Chr_code_2	chrcode
	Ref_3	Chr_code_3	chrcode
    ...
	# ------------------------------------------------------------
	# ---------- BLOCK D (Only if there is an outgroup) ----------
 	# ------------------------------------------------------------
  
    # ID_outgroup	BED_file/Chr_code_outgroup	config_synteny_tag
	
	outg	PATH/Ref_master.outg.bed	WGA
    outg	Chr_code_outg	chrcode

    # --------------------------------------- END ---------------------------------------
 
    # original_sample_header (1st column; Block A): For each sample, the sample name as shown in the input VCF file.
    # final_sample_header (2nd column; Block A): For each sample, the user-chosen sample name to be displayed in downstream results.
    # config_tag (3rd column; Block A): For each sample, a mandatory tag of block A (real_name) that must appear for `vcf2synteny` correct processing.
	
	# ID_ref (1st column; Block B): Abbreviation (3-4 letters; e.g., Aet for Aegilops tauschii or Bsta for Brachypodium stacei) for each secondary reference genome (and subgenomes in downstream output).
    # BED_file (2nd column; Block B): For each secondary reference genome, BED file (and path) obtained in step 2.1.
	# config_synteny_tag (3rd column; Block B): For each secondary reference genome, a mandatory tag of Block B (WGA) that must appear for `vcf2synteny` correct processing.
    
	# ID_ref (1st column; Block C): Abbreviation (3-4 letters) for master and each secondary reference genome (and subgenomes in downstream output).
    # Chr_code (2nd column; Block C): Regular expressions (e.g., Chr(\d+)) to match chromosome names from master and secondary reference genomes used in step 2.1, can use those proposed by WGA
	# config_synteny_tag (3rd column; Block C): For master and each secondary reference genome, a mandatory tag of Block B (chrcode) that must appear for `vcf2synteny` correct processing.

    # ID_outgroup (1st column; Block D): outg abbreviation for outgroup species.
	# BED_file/Chr_code_outgroup (2nd column; Block D): First line: For outgroup genome, BED file (and path) obtained in step 2.1.; Second line: Regular expression to match chromosome names from the outgroup genome used in step 2.1.
    # config_synteny_tag (3rd column; Block D): A mandatory tag of Block C, first line (WGA) and second line (chrcode), that must appear for `vcf2synteny` correct processing.
 
	

In our example, we use a toy VCF file that contains read-mapping positions on chromosomes Bd2 of *Brachypodium distachyon* and Chr01 of *Brachypodium stacei*, the LOG file computed in 2.2 section by `vcf2alignment` and the syntenic positions computed in 2.1 section by `WGA` (see config file for `vcf2synteny`):
 
    ./vcf2synteny -v sample_data/BdisBd2_BstaChr01.vcf.gz -c sample_data/config.synteny.tsv -l BdisBd2_BstaChr01.vcf.log.gz \
			-d 5 -m 3 -r Bdis -o BdisBd2_BstaChr01.DP5.M3.synteny.fasta

Note that a different config file is now used (see example 
[sample_data/config.synteny.tsv](https://github.com/eead-csic-compbio/vcf2alignment/blob/master/sample_data/config.synteny.tsv)),
which also contains: 

+ a path to the BED file obtained in step 2.1
+ regular expressions to match chromosome names from reference genomes used in step 2.1, can use those proposed by `WGA`

Note that `vcf2synteny` performs several sort operations. 
With large genomes, these might require significant disk space to save temporary results.
By default, these are stored in `/tmp`, but this can be changed with the flag `-t`. 
The examples in the [Makefile](./Makefile) use `-t` 
pointing to the same **output folder** used by `WGA`, so that all files are contained there and can be safely removed if needed.
 
Anyway, this script produces the following output:

    # computing Bdis.Bsta.coords.positions.tsv (3 steps)

    # master reference: Bdis
    # secondary references: Bsta
    # synteny files (SYNTENYZEROBASED=1): 
    # Bsta : _Bdis.Bsta.coords.positions.tsv...
    # total positions=762307

    # decompressing VCF file with GZIP
    # number of samples found=6
    # number of loci read from VCF: 10000
    ...
    # number of loci read from VCF: 760000
    # sorting SNPs by position ...
    # aligned position: Bd2_48166443 : Chr01_1718238,
    ...
    # aligned position: Bd2_58679964 : Chr01_384888,  
    # number of valid loci=762307
    # number of polymorphic loci=9512

    # Bdis_ABR2_Bdis variants: 0 / 762307
    # Bdis_ABR2_Bsta variants: 8471 / 762307
    # Bdis_Bd21Control_Bdis variants: 0 / 762307
    # Bdis_Bd21Control_Bsta variants: 4048 / 762307
    # Bhyb_Bhyb26_Bdis variants: 0 / 762307
    # Bhyb_Bhyb26_Bsta variants: 741493 / 762307
    # Bhyb_ABR113_Bdis variants: 0 / 762307
    # Bhyb_ABR113_Bsta variants: 677818 / 762307
    # Bsta_ABR114_Bdis variants: 0 / 762307
    # Bsta_ABR114_Bsta variants: 728488 / 762307
    # Bsta_TE4.3_Bdis variants: 0 / 762307
    # Bsta_TE4.3_Bsta variants: 754853 / 762307

The resulting multiple sequence alignment (MSA) has as many lines per sample as references, which are handled as subgenomes (and artifactual subgenomes).
The first 200 positions of the MSA derived from the sample data look as follows:

    >Bdis_ABR2_Bdis
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCNNNNANNCNNNNNNNNNNNNNNNCNTNNNNTNNNNNNNNNNNNNGNTNNNNNTNNGGCNNNNNNNNNNNNNNGANANCC
    >Bdis_ABR2_Bsta
    ATCTCGCGGCTGCCCCCCCGAGTTCGGAGGACCACCGGCCTCGCCGGGCCTCCACGTCGTGGCCACTGCTTCGGCAACGCACTGCTGACCTCACCGCTGCCGTCGCACTGGCAAACGGGTCAGCAAATCAAGCGCGCTGCGTGTCGTCNCGTCTCGGGCCATGCCGCTTTTCATCTGGCCGCCCTGGTTGCGCGACACCC
    >Bdis_Bd21Control_Bdis
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCNNNNANNCNNNNNNNNNNNNNNNCNTNNNNTNNNNNNNNNNNNNGNTNNNNNTNNGGCNNNNNNNNNNNNNNGANANCC
    >Bdis_Bd21Control_Bsta
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    >Bhyb_Bhyb26_Bdis
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCNNNNANNCNNNNNNNNNNNNNNNCNTNNNNTNNNNNNNNNNNNNGNTNNNNNTNNGGCNNNNNNNNNNNNNNGANANCC
    >Bhyb_Bhyb26_Bsta
    ATCTCGCGGCTGCCCCCCCGAGTTCGGAGGACCACCGGCCTCGCCGGGCCTCCACGTCGTGGCCACTGCTTCGGCAACGCACTGCTGACCTCACCGCTGCCGTCGCACTGGCAAACGGGTCAGCAAATCAAGCGCGCTGCGTGTCNTCGCGTCTCGGGCCATGCCGCTTTTCATCTGGCCGCCCTGGTTGCGCGACACCC
    >Bhyb_ABR113_Bdis
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    >Bhyb_ABR113_Bsta
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    >Bsta_ABR114_Bdis
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGNTNNNNNTNNGGCNNNNNNNNNNNNNNGANANCC
    >Bsta_ABR114_Bsta
    ATCTCGCGGCTGCCCCCCCGAGTTCGGAGGACCACCGGCCTCGCCGGGCCTCCACGTCGTGGCCACTGCTTCGGCAACGCACTGCTGACCTCACCGCTGCCGTCGCACTGGCAAACGGGTCAGCAAATCAAGCGCGCTGCGTGTCGTCGCGTCTCGGGCCATGCCGCTTTTCATCTGGCCGCCCTGGTTGCGCGACACCC
    >Bsta_TE4.3_Bdis
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    >Bsta_TE4.3_Bsta
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAATCAAGCGCGCTGCGTGTCGTCGCGTCTCGGGCCATGCCGCTTTTCATCTGGCCGCCCTGGTTGCGCGACACCC

NOTE: **Artifactual subgenomes** must be eliminated from the final multiple sequence alignment (MSA). If the ploidy, and therefore the number of subgenomes expected to be recovered, is unknown, we propose using a cross-validation criterion based on the percentage of SHPs recovered in diploid samples. In diploid samples, only one predominant genome/subgenome is expected to be recovered. For each subgenome, the highest SHP percentage obtained in the diploid samples from the non-specific mappings will be set as a threshold. *More details in the publication*.

This completes this protocol.


### 2.4) Other uses: get a multiple sequence alignment of samples in a VCF file

Working with collaborators, we have noticed that often you need to produce a multi-FASTA file out of the samples in a VCF file.
A typical reason for this is to compute a phylogenetic tree. Please note that `vcf2alignment` can be used exactly for this, as follows,
note the optional -o flag. If the input VCF file does not contain DP data, or you want to take all base calls regardless of their depth,
add also -d: 

    ./vcf2alignment -v sample_data/BdisBd2_BstaChr01.vcf.gz -c sample_data/config.tsv -l BdisBd2_BstaChr01.log.gz -d 0 -o BdisBd2_BstaChr01.fasta

## 3) Citation

Main paper:

Sancho R, Catalan P, Vogel JP, Contreras-Moreira B (2025) Deconvolution of Single Homeologous Polymorphism (SHP) drives phylogenetic analysis of allopolyploids.
bioRxiv 2025.07.17.665301; doi: https://doi.org/10.1101/2025.07.17.665301 

Key dependencies:

Lin HN, Hsu WL (2020) GSAlign: an efficient sequence alignment tool for intra-species genomes. BMC Genomics 21:182. https://doi.org/10.1186/s12864-020-6569-1

Nakato R, Gotoh O (2010) Cgaln: fast and space-efficient whole-genome alignment. BMC Bioinformatics 11:224. https://doi.org/10.1186/1471-2105-11-224

Girgis HZ (2015) Red: an intelligent, rapid, accurate tool for detecting repeats de-novo on the genomic scale. BMC Bioinformatics 16:227. https://doi.org/10.1186/s12859-015-0654-5

Contreras-Moreira B, Filippi CV, Naamati G, García Girón C, Allen JE, Flicek P (2021) Efficient masking of plant genomes by combining kmer counting and curated repeats Genomics. Plant Genome https://doi.org/10.1002/tpg2.20143

