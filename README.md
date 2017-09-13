# vcf2alignment

Protocol to produce multiple sequence alignments out of [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)
files which can be used for phylogenetic tree construction. 

**Authors**
Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)

1. Estación Experimental de Aula Dei-CSIC, Zaragoza, Spain
2. Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
3. Fundación ARAID, Zaragoza, Spain

## Pipeline overview

<!-- flowchart -->

### Software dependencies

This protocol has been tested on Linux x86_64 systems, although it should also work on Mac-OSX settings.
It requires Perl5, which should be installed on all Linux environments, plus some standard programs (gzip, bzip2).


## 1) Input data 

These are the data required to run this pipeline:

+ 1+ (ideally more) reference genomes of species of the taxa of interest, **concatenated in a single FASTA file**.
In our *Brachypodium* benchmark we only used complete chromosome arms, and left out single contigs and centromeric parts. 
This might require renaming chromosomes to make sure that each species has unique names.

+ **FASTQ files of sequence reads** from the samples to be analyzed, which should belong to the same taxa, and perhaps some outgroups as well.

## 2) Simple mode: mapped reads

### 2.1) Read mapping 

<!-- Explicar los mapeos con BWA mem o Hisat2, segun sea

Explicar que habra muestras con depth of coverage mas limitada y otras mejores,
aparte de otras que se pueden definir como outgroups para los arboles poesteriores
 -->

### 2.2) Merging BAM files to produce a single non-redundant VCF file

<!-- Explicar los comandos para ir desde los multiples SAM a un solo VCF --> 

```{shell}
utils/rm_double_lines.pl RNAseq_Bd5_Chr10_chr10.raw.vcf > sample_data/RNAseq_Bd5_Chr10_chr10.vcf
bzip2 sample_data/RNAseq_Bd5_Chr10_chr10.vcf
```

### 2.3) Producing a multiple alignment file

Script [vcf2alignment.pl](./vcf2alignment.pl) ships with the following global variables which might be modified 
in the source code (with a text editor) to change the expected outcome:

| variable name | default value | definition |
|:-----|:---------------:|:-------|
| $MINDEPTHCOVERPERSAMPLE | 10 | natural, min number of reads mapped supporting a locus |
| $MAXMISSINGSAMPLES | 8 | natural, max number of missing samples accepted per locus |
| $ONLYPOLYMORPHIC | 1 | set to 0 to keep fixed loci, helps with sparse data |
| $OUTFILEFORMAT | fasta | currently can also take nexus and phylip format |

By default, the produced alignment is in FASTA format. Note that a logfile is also saved in this example,
which contains a list of valid loci and several statistics:
```{shell}
./vcf2alignment.pl sample_data/RNAseq_Bd5_Chr10_chr10.vcf.bz2 \
  sample_data/RNAseq_Bd5_Chr10_chr10.fna &> sample_data/RNAseq_Bd5_Chr10_chr10.log 
```

If the format is changed in the source to 'phylip', 
an [interleaved PHYLIP](http://evolution.genetics.washington.edu/phylip/doc/sequence.html) 
file is produced, with names shortened to 10 chars (see source code to choose from prefixes or suffixes):
```{shell}
./vcf2alignment.pl sample_data/RNAseq_Bd5_Chr10_chr10.vcf.bz2 \
  sample_data/RNAseq_Bd5_Chr10_chr10.phy &> sample_data/RNAseq_Bd5_Chr10_chr10.log
```

As mentioned, [NEXUS](https://en.wikipedia.org/wiki/Nexus_file) alignments can also be produced 
by editing the source to 'nexus':
```{shell}
./vcf2alignment.pl sample_data/RNAseq_Bd5_Chr10_chr10.vcf.bz2 \
  sample_data/RNAseq_Bd5_Chr10_chr10.nex &> sample_data/RNAseq_Bd5_Chr10_chr10.log
```

The produced multiple alignment should be rendered with appropriate software for visual quality check:
![Multiple alignment generated](./pics/MSA_simple.png)


## 3) Advanced mode: syntenic coordinates + mapped reads

### 3.1) Whole-genome alignments

Alignments must be computed to find syntenic segments among the reference genomes available for read mapping.
We selected [CGaln](http://www.iam.u-tokyo.ac.jp/chromosomeinformatics/rnakato/cgaln/index.html) for this task,
which requires the input sequences to be [soft-masked](https://genomevolution.org/wiki/index.php/Masked) ahead.
These files are not included here as they're bulky, but we do show how we use processed them in our benchmark:

```{shell, eval=FALSE}
# index individual references (n=3)
~/soft/Cgaln/maketable Bdistachyon_msk.fna
~/soft/Cgaln/maketable Bstacei_msk.fna
~/soft/Cgaln/maketable Bsylvaticum_msk.fna

# alignments with custom parameters

~/soft/Cgaln/Cgaln Bdistachyon_msk.fna Bstacei_msk.fna \
  -o Bdistachyon.Bstacei.block12K.aln.hq.fna -r -X12000 -fc -cons -otype2

~/soft/Cgaln/Cgaln Bdistachyon_msk.fna Bsylvaticum_msk.fna \
  -o Bdistachyon.Bsylvaticum.block12K.aln.hq.fna -r -X12000 -fc -cons -otype2 
```

We recomend that users visualize the alignments to make sure they make sense and to optimize CGaln parameters. 
This can be done producing *.dot* files instead of FASTA output:
```{shell, eval=FALSE}
~/soft/Cgaln/Cgaln Bdistachyon_msk.fna Bstacei_msk.fna \
  -o Bdistachyon.Bstacei.block12K.aln.hq.dot -r -X12000 -fc -cons 

~/soft/Cgaln/Cgaln Bdistachyon_msk.fna Bsylvaticum_msk.fna \
  -o Bdistachyon.Bsylvaticum.block12K.aln.hq.dot -r -X12000 -fc -cons 
```
These files can then be inspected with gnuplot:

```{shell, eval=FALSE}
gnuplot
gnuplot> plot "result.dot" with lines
```
This an example whole-genome alignment plot:
![whole-genome alignment plot](./pics/dotplot.png)

These alignments can then be compressed and equivalent/syntenic positions extracted as follows: 
```{shell, eval=FALSE}
utils/mapcoords.pl Bdistachyon.Bstacei.block12K.aln.hq.fna.gz Bdistachyon_msk.fna Bstacei_msk.fna \
  > Bdistachyon.Bstacei.coords.tsv 2> Bdistachyon.Bstacei.coords.log

utils/mapcoords.pl Bdistachyon.Bsylvaticum.block12K.aln.hq.fna.gz Bdistachyon_msk.fna Bsylvaticum_msk.fna \
  > Bdistachyon.Bsylvaticum.coords.tsv 2> Bdistachyon.Bsylvaticum.coords.log
```

### 3.2) Read mapping 

<!-- Explicar los mapeos con BWA mem o Hisat2, segun sea

Explicar que habra muestras con depth of coverage mas limitada y otras mejores,
aparte de otras que se pueden definir como outgroups para los arboles poesteriores
 -->

### 3.3) Merging BAM files to produce a single non-redundant VCF file

<!-- Explicar los comandos para ir desde los multiples SAM a un solo VCF 

```{shell}
utils/rm_double_lines.pl RNAseq_Bd5_Chr10_chr10.raw.vcf > sample_data/RNAseq_Bd5_Chr10_chr10.vcf
bzip2 sample_data/RNAseq_Bd5_Chr10_chr10.vcf
```
-->

### 3.4) Producing a multiple alignment file

