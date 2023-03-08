#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

# converts an input VCF file, with might be GZIP/BZIP2 compressed, into a  
# multiple sequence alignment in several supported formats 
# (check @validformats below)

# Note: Chromosome names of genomes in VCF file must be different.

# Bruno Contreras, Ruben Sancho EEAD-CSIC 2017-2023

my %IUPACdegen = (  
  'AG'=>'R', 'GA'=>'R', 'CT'=>'Y', 'TC'=>'Y',
  'CG'=>'S', 'GC'=>'S', 'AT'=>'W', 'TA'=>'W',
  'GT'=>'K', 'TG'=>'K', 'AC'=>'M', 'CA'=>'M' );

my %revcomp = ('A'=>'T', 'T'=>'A','G'=>'C','C'=>'G', 'N'=>'N');

my @validformats = qw( phylip nexus fasta );

# VCF filtering and output options, edit as required
my $MINDEPTHCOVERPERSAMPLE = 3; # natural, min number of reads mapped supporting a locus
my $MAXMISSINGSAMPLES      = 8;  # natural, max number of missing samples accepted per locus
my $ONLYPOLYMORPHIC        = 0;  # set to 0 to keep fixed loci, helps with sparse data
my $OUTFILEFORMAT          = 'fasta'; # can also take other formats in @validformats

# first guess of key VCF columns, adjusted in real time below
my $COLUMNFIRSTSAMPLE      = 9; # zero-based, VCF format http://www.1000genomes.org/node/101, format is previous one
my $GENOTYPECOLUMNFORMAT   = 0; # zero-based, initial guess, they are set for each line
my $DEPTHCOLUMNFORMAT      = 1; 

# external binaries, edit if not installed elsewhere and not in path
my $GZIPEXE  = 'gzip'; 
my $BZIP2EXE = 'bzip2';

my %opts;
my ($filename,$configfile,$outfilename) = ('','','');
my ($mindepth,$maxmissing,$only_polymorphic,$outformat) = 
  ($MINDEPTHCOVERPERSAMPLE,$MAXMISSINGSAMPLES,$ONLYPOLYMORPHIC,$OUTFILEFORMAT);

getopts('hpv:c:o:d:m:f:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-v input VCF file              (example: -v data.vcf.gz)\n";
  print "-c input TSV config file       (example: -c config.tsv)\n";
  print "-o output file name            (example: -o out.fasta)\n";
  print "-d min depth of called SNPs    (optional, example: -d 10, default -d $mindepth)\n";
  print "-m max missing samples         (optional, example: -m 10, default -m $maxmissing\n";
  print "-f output format               (optional, example: -f nexus, default -f $outformat)\n";
  print "-p take only polymorphic sites (optional, by default all sites are taken)\n";

  exit(0);
}

if(!defined($opts{'v'}))
{
  die "# ERROR: need input VCF file (-v)\n";
} else {
  $filename = $opts{'v'}
}

if(!defined($opts{'c'}))
{
  die "# ERROR: need input TSV config file (-c)\n";
} else {
  $configfile = $opts{'c'}
}

if(!defined($opts{'o'}))
{
  die "# ERROR: need output filename (-o)\n";
} else {
  $outfilename = $opts{'o'}
}

if(defined($opts{'d'}) && $opts{'d'} > 0)
{
  $mindepth = $opts{'d'}
}

if(defined($opts{'m'}) && $opts{'d'} >= 0)
{
  $maxmissing = $opts{'m'}
}

if(defined($opts{'f'}) && grep(/^$opts{'f'}$/,@validformats))
{
  $outformat = $opts{'f'}
}

if(defined($opts{'p'}))
{
  $only_polymorphic = 1
}

warn "# $0 -v $filename -c $configfile -o $outfilename -d $mindepth -m $maxmissing -f $outformat -p $only_polymorphic\n\n";

######################################################

my ($n_of_samples,$n_of_loci,$n_var_loci,$depthcover,$missing,$genotype,$allele) = (0,0,0);
my ($corr_coord,$sample,$lastsample,$idx,$lastsampleidx,$file);
my (@samplenames,@MSA,%MSAref,%stats,%refallele,%refstrand);
my ($snpname,$badSNP,$shortname,$magic,%contigstats);
my %vcf_real_names; # To shorten sample names in output alignment
my %genomic_samples;# Set samples which should not count as missing data.
                    # For instance, we used it to leave outgroups out of these calculations,
                    # as their WGS reads are significantly deeper than GBS/RNAseq samples 

## parse config file
open(CONFIG,"<$configfile") || die "# cannot read $configfile\n";
while(my $line = <CONFIG>)
{
  #./Bhyb127/Bhyb127.sort.q30.bam  Bhyb_Bhyb127    real_name
  #syl_Cor_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam        1       deep_sample  
  chomp($line);
  my @cdata = split(/\t/,$line);
  if($cdata[2] eq 'real_name')
  {
    $vcf_real_names{ $cdata[0] } = $cdata[1];
  } 
  elsif($cdata[2] eq 'deep_sample')
  {
    $genomic_samples{ $cdata[0] } = 1
  } 
  else
  {
    print "# unrecognized configuration: $line\n";
  }
}
close(CONFIG);

## check input file and choose right way to parse input lines
my ($genomic_samples,$gbs_samples); 
open(INFILE,$filename) || die "# cannot open input file $filename, exit\n";
sysread(INFILE,$magic,2);
close(INFILE);

if($filename =~ /\.gz$/ || $magic eq "\x1f\x8b") # compressed input
{
  printf(STDERR "# decompressing VCF file with GZIP\n");
  open(VCF,"$GZIPEXE -dc $filename |");
}
elsif($filename =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
{
  if(!open(VCF,"$BZIP2EXE -dc $filename |"))
  {
    die "# cannot read BZIP2 compressed $filename $!\n"
      ."# please check $BZIP2EXE is installed\n";
  }
}
else{ open(VCF,'<',$filename) }
  
while(<VCF>)   
{  
  chomp($_);
  my @rawdata = split(/\t/,$_);

  if($n_of_samples > 0)
  {
    #Bd1  346  .  A  C  999  PASS  DP=4008;VDB=5.239035e-01;...  GT:PL:DP:SP:GQ  1/1:255,255,0:185:0:99  ...

    # skip non-polymorphic sites if required
    next if($rawdata[4] eq '.' && $only_polymorphic);

    # skip indels
    next if($rawdata[3] =~ m/[A-Z]{2,}/ || $rawdata[4] =~ m/[A-Z]{2,}/); 
  
    # skip multiallelic SNPs
    #next if(length($rawdata[4]) > 1);

    # find out which data fields to parse
    my @sampledata = split(/:/,$rawdata[$COLUMNFIRSTSAMPLE-1]);
    foreach my $sd (0 .. $#sampledata)
    {
      if($sampledata[$sd] eq 'GT'){ $GENOTYPECOLUMNFORMAT = $sd; last }
      $sd++;
    }
    foreach my $sd (0 .. $#sampledata)
    {
      if($sampledata[$sd] eq 'DP'){ $DEPTHCOLUMNFORMAT = $sd; last }
      $sd++;
    } 

    my (@sequence);
    my %nts; # ($rawdata[3] => 1 ); # adds ref base call
    $sample=$depthcover=$missing=$badSNP=0;
    ($genomic_samples,$gbs_samples) = (0,0);
    foreach $idx ( $COLUMNFIRSTSAMPLE .. $lastsampleidx )
    {
      #0/0:0,255,255:93:0:99
      #1/1:255,255,0:185:0:99
      #0/1:236,0,237:66:7:9
      my @sampledata = split(/:/,$rawdata[$idx]); 
      $genotype = $sampledata[$GENOTYPECOLUMNFORMAT];
      $depthcover = $sampledata[$DEPTHCOLUMNFORMAT]; 
      $allele = 'N'; # default 
      
      if($genotype eq '0/0')
      {
        if($depthcover >= $mindepth){ 
          $allele = $rawdata[3]; 
          if(defined($genomic_samples{$samplenames[$sample]})){ $genomic_samples++ }
          else{ $gbs_samples++ }
        }
        else{ if(!$genomic_samples{$samplenames[$sample]}){ $missing++ } }
      }
      elsif($genotype eq '1/1')
      {
        if($depthcover >= $mindepth){ 
          $allele = (split(/,/,$rawdata[4]))[0]; 
          if(defined($genomic_samples{$samplenames[$sample]})){ $genomic_samples++ }
          else{ $gbs_samples++ }
        }
        else{ if(!$genomic_samples{$samplenames[$sample]}){ $missing++ } } 
      }
      elsif($genotype eq '2/2')
      {
        if($depthcover >= $mindepth){ 
          $allele = (split(/,/,$rawdata[4]))[1]; 
          if(defined($genomic_samples{$samplenames[$sample]})){ $genomic_samples++ }
          else{ $gbs_samples++ }
        }
        else{ if(!$genomic_samples{$samplenames[$sample]}){ $missing++ } }
      }
      elsif($genotype eq '3/3')
      {
        if($depthcover >= $mindepth){     
          $allele = (split(/,/,$rawdata[4]))[2]; 
          if(defined($genomic_samples{$samplenames[$sample]})){ $genomic_samples++ }
          else{ $gbs_samples++ }
        }
        else{ if(!$genomic_samples{$samplenames[$sample]}){ $missing++ } }
      }
      else # missing or htzg are treated as missin all the same
      {
        if(!$genomic_samples{$samplenames[$sample]}){ $missing++ }
      }
  
      if($missing > $maxmissing)
      {
        $badSNP = 1;
        last;
      }  
      
      $sequence[$sample] = $allele;
      if($allele ne 'N'){ $nts{$allele}++ }

      $sample++;
    } 

    # make sure genomic-only sites are skipped
    if($gbs_samples == 0) 
    {
      $badSNP = 1;
    } 

    # make sure monomorphic sites are skipped
    if(scalar(keys(%nts)) < 2 && $only_polymorphic){ $badSNP = 1 }
    
    if(!$badSNP)
    {
      if($sample != $n_of_samples)
      {
        die "# VCF line contains less samples than expected ($sample < $n_of_samples):\n$_\n";
      }
    
      $snpname = "$rawdata[0]_$rawdata[1]";
      
      printf(STDERR "# valid locus: $snpname $missing ".join(',',keys(%nts))."\n");
      
      foreach $sample (0 .. $lastsample)
      {
        $MSA[$sample] .= $sequence[$sample];

        # save stats of missing data
        if($sequence[$sample] eq 'N')
        {
          $stats{$sample}{'totalNs'}++;
          $contigstats{$sample}{$rawdata[0]}{'N'}++;
        }
        else
        {
          $stats{$sample}{'total'}++;
          $contigstats{$sample}{$rawdata[0]}{'SNP'}++;
        }  
      }

      if(scalar(keys(%nts)) > 1){ $n_var_loci++ }
      $n_of_loci++;
    }  
  }
  elsif($rawdata[0] eq '#CHROM')
  {
    #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1  sampleN
    push(@samplenames,@rawdata[9 .. $#rawdata]);
    $n_of_samples = scalar(@samplenames);
    $lastsample = $n_of_samples-1;
    $lastsampleidx = $#rawdata;
    printf(STDERR "# number of samples found=$n_of_samples\n");
  }
}  
close(VCF);

printf(STDERR "# number of valid loci=$n_of_loci\n");
if(!$only_polymorphic)
{
  printf(STDERR "# number of polymorphic loci=$n_var_loci\n");  
} warn "\n";

open(OUTFILE,">$outfilename") || die "# cannot create output file $outfilename, exit\n";


if($outformat eq 'nexus')
{
  print OUTFILE "#NEXUS\nBEGIN DATA;\nDIMENSIONS  NTAX=$n_of_samples NCHAR=$n_of_loci;\n";
  print OUTFILE "FORMAT DATATYPE=DNA  MISSING=N GAP=-;\nMATRIX\n";
}
elsif($outformat eq 'phylip')
{
  print OUTFILE "$n_of_samples    $n_of_loci\n";    
}

foreach $sample (0 .. $lastsample)
{
  if($vcf_real_names{$samplenames[$sample]}){
    $shortname = $vcf_real_names{$samplenames[$sample]} 
  } 
  else{ $shortname = $samplenames[$sample]; }

  if($outformat eq 'nexus')
  {
    print OUTFILE "$shortname $MSA[$sample]\n";
  }
  elsif($outformat eq 'phylip')
  {
    if(length($shortname)>10){ 
      $shortname = substr($shortname,0,9).'_'; # prefix
      #$shortname = '_'.substr($shortname,-9); # suffix
      printf(STDERR "# phylip sample name shortened: $samplenames[$sample] -> $shortname\n");
    }
    print OUTFILE "$shortname    $MSA[$sample]\n";
  }
  elsif($outformat eq 'fasta')
  {
    print OUTFILE ">$shortname\n";
    print OUTFILE "$MSA[$sample]\n";
  }
} 

if($outformat eq 'nexus')
{
  print OUTFILE ";\nEND;\n";
}

close(OUTFILE);

printf(STDERR "\n\n# stats (#SNPs):\n");

foreach $sample (0 .. $lastsample)
{
  printf(STDERR "$samplenames[$sample] : $stats{$sample}{'total'}\n");
} 

printf(STDERR "\n# stats per contig/chr (#SNPs):\n");

foreach $sample (0 .. $lastsample)
{
  foreach my $contig (sort keys(%{$contigstats{$sample}}))
  {
    printf(STDERR "%s\t%s\t%d\n",
      $samplenames[$sample],$contig,$contigstats{$sample}{$contig}{'SNP'} || 0);
  }
}

printf(STDERR "\n# stats (#N):\n");

foreach $sample (0 .. $lastsample)
{
  printf(STDERR "$samplenames[$sample] : $stats{$sample}{'totalNs'}\n");
}

printf(STDERR "\n# stats per contig/chr (#N):\n");

foreach $sample (0 .. $lastsample)
{
  foreach my $contig (sort keys(%{$contigstats{$sample}}))
  {
    printf(STDERR "%s\t%s\t%d\n",
      $samplenames[$sample],$contig,$contigstats{$sample}{$contig}{'N'} || 0);
  }
}


