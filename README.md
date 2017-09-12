# vcf2alignment

Protocol to produce multiple sequence alignments out of [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)
files which can be used for phylogenetic tree construction. 

**Authors**
Ruben Sancho (1,2), Bruno Contreras Moreira (1,3)

1. Estación Experimental de Aula Dei-CSIC, Zaragoza, Spain
2. Escuela Politécnica Superior de Huesca, U.Zaragoza, Spain
3. Fundación ARAID, Zaragoza, Spain

## Flowchart

## Simple mode: mapped reads

### Read mapping 

<!-- Explicar los mapeos con BWA mem o Hisat2, segun sea -->

### Merging BAM files to produce a single non-redundant VCF file

<!-- Explicar los comandos para ir desde los multiples SAM a un solo VCF --> 

```{shell}
./_rm_double_lines.pl RNAseq_Bd5_Chr10_chr10.raw.vcf > sample_data/RNAseq_Bd5_Chr10_chr10.vcf
bzip2 sample_data/RNAseq_Bd5_Chr10_chr10.vcf
```



## Advanced mode: mapped reads + syntenic chromosome coordinates
