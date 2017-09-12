# vcf2alignment

Protocol to produce multiple sequence alignments out of [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)
files which can be used for phylogenetic tree construction. 

**Authors**
Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)

1. Estación Experimental de Aula Dei-CSIC, Zaragoza, Spain
2. Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
3. Fundación ARAID, Zaragoza, Spain

## 0) Software dependencies

This protocol has been tested on Linux x86_64 systems, although it should also work on Mac-OSX settings.
It requires Perl5, which should be installed on all Linux environments, plus some standard programs (gzip, bzip2).

## 1) Flowchart

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



## 3) Advanced mode: mapped reads + syntenic chromosome coordinates
