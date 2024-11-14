# vcf2alignment

This protocol computes Whole Genome Alignments (WGA) to discover syntenic SNPs out of reads mapped to concatenated 
genome references. It requires FASTA and [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) input files and 
produces multiple sequence alignments of subgenomes that make up polyploids.

Rubén Sancho (1,2), Pilar Catalán (2), Bruno Contreras Moreira (1,3)

1. Estación Experimental de Aula Dei-CSIC, Zaragoza, Spain
2. Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
3. Fundación ARAID, Zaragoza, Spain

[![Build Status](https://github.com/eead-csic-compbio/vcf2alignment/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/eead-csic-compbio/vcf2alignment/actions/workflows/ci.yml)

## Installation

While this protocol was developed on Linux x86_64 systems, it should also work on MacOS settings.
It requires some standard Linux utilities (gzip, sort, perl, make, python3, g++, etc) and a few third-party dependencies
which can be installed locally as follows:

    git clone https://github.com/eead-csic-compbio/vcf2alignment.git
    cd vcf2alignment
    make install
    # optionally, takes a couple minutes
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


## 2) How to run

A typical analysis involves calling three scripts: [WGA](https://github.com/eead-csic-compbio/vcf2alignment/blob/master/WGA), 
[vcf2alignment](https://github.com/eead-csic-compbio/vcf2alignment/blob/master/WGA) and 
[vcf2synteny](https://github.com/eead-csic-compbio/vcf2alignment/blob/master/WGA). 

Note that the steps below can be run in one go as follows:

    make test

### 2.1) Whole-genome alignments (WGA)

WGAs must be computed to find syntenic segments among the reference genomes available for read mapping.
By default this uses [CGaln](https://github.com/rnakato/Cgaln),
which requires the input sequences to be [soft-masked](https://genomevolution.org/wiki/index.php/Masked) ahead.
In our example we set *B. distachyon* as the master reference genome for being the best quality assembly at hand
(see 2.3 below).
In this we now find out syntenic segments on the other genomes, defined as secondary genomes (*B. stacei*). 
Each individual reference is hence considered a **subgenome** to which reads map:

    ./WGA -A sample_data/Bdis.fna.gz -B sample_data/Bsta.fna.gz

This produces the following output:

    ## WGA -A sample_data/Bdis.fna.gz -B sample_data/Bsta.fna.gz -l 1 -m 1 -G 0 -I -K11 -BS10000 -C -X12000 -fc -cons -n 4

    ## root: Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X12000_-fc_-cons

    # filter_FASTA_sequences: [passed 59130575 bp] Bd2
    # filter_FASTA_sequences: [passed 31564145 bp] Chr01

    # chrcode A example: Bd(\d+)
    # chrcode B example: Chr(\d+)

    ## soft-masking filtered sequences

    ## CGaln algorithm


    ## indexing masked, filtered sequences

    [...]

    ## computing and plotting Whole Genome Alignment

    [...]

    ## converting Whole Genome Alignment to BED

    ## output files:

    # BED: Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X12000_-fc_-cons.bed
    # LOG: Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X12000_-fc_-cons.coords.log
    # PDF: Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X12000_-fc_-cons.dot.pdf

    ## WGA summary: valid blocks: 54 unique positions: 8956568

    # time used (s): 119 memory used (Mb): 407.9

The most important result files are the 0-based **BED** list of syntenic positions, which will be used in the last step,
and the **PDF** dotplot, which requires `gnuplot` in your syste, which must be inspected to assess the quality of
the WGA. Note that flags -I and -C can be used to tweak the WGA parameters after inspection of the dotplot.

Alternatively, the GSAlign WGA algorithm can be invoked as follows, adding optionally flag -G to
adapt the default parameters to your genomes of interest:
 
    ./WGA -A sample_data/Bdis.fna.gz -B sample_data/Bsta.fna.gz -g

![whole-genome alignment plot](./pics/dotplot.png)
*Figure 1. Example WGA dotplot.* A good WGA dotplot will have long diagonal runs of aligned genome regions instead of clouds.

### 2.2) Filtering valid positions in the VCF file

Script `vcf2alignment` parses an input VCF file, with might be GZIP/BZIP2 compressed, 
and produces a list of valid sites considering min read depth (-d) and max missing data (-m).
Note this requires a config file that matches sample names in the VCF file to their human-readable names
(see example [sample_data/config.tsv](https://github.com/eead-csic-compbio/vcf2alignment/blob/master/sample_data/config.tsv)).
A report log file with valid 1-based coordinates (-l) is saved to be used in the last step:

    ./vcf2alignment -v sample_data/BdisBd2_BstaChr01.vcf.gz -c sample_data/config.tsv -l BdisBd2_BstaChr01.vcf.log -d 5 -m 3


### 2.3) Producing a multiple sequence alignment of polyploid subgenomes

Script `vcf2synteny` puts it all together and produces an alignment in FASTA format:
 
    ./vcf2synteny -v sample_data/BdisBd2_BstaChr01.vcf.gz -c sample_data/config.synteny.tsv -l BdisBd2_BstaChr01.vcf.log \
			-d 5 -m 3 -r Bdis -o BdisBd2_BstaChr01.DP5.M3.synteny.fasta

Note that a different config file is now used (see example 
[sample_data/config.synteny.tsv](https://github.com/eead-csic-compbio/vcf2alignment/blob/master/sample_data/config.synteny.tsv)),
which also contains: 

+ a path to the BED file obtained in step 2.1
+ regular expressions to match chromosome names from reference genomes used in step 2.1, can use those proposed by `WGA`

The resulting multiple alignment has as many lines per sample as references, which are handled as subgenomes:
![Multiple alignment generated](./pics/MSA_subgenomes.sample.png)



## 3) Citation

this paper, to be written

Lin HN, Hsu WL (2020) GSAlign: an efficient sequence alignment tool for intra-species genomes. BMC Genomics 21:182. https://doi.org/10.1186/s12864-020-6569-1

Nakato R, Gotoh O (2010) Cgaln: fast and space-efficient whole-genome alignment. BMC Bioinformatics 11:224. https://doi.org/10.1186/1471-2105-11-224

Girgis HZ (2015) Red: an intelligent, rapid, accurate tool for detecting repeats de-novo on the genomic scale. BMC Bioinformatics 16:227. https://doi.org/10.1186/s12859-015-0654-5

Contreras-Moreira B, Filippi CV, Naamati G, García Girón C, Allen JE, Flicek P (2021) Efficient masking of plant genomes by combining kmer counting and curated repeats Genomics. Plant Genome https://doi.org/10.1002/tpg2.20143

