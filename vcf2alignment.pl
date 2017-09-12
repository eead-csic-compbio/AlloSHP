#!/usr/bin/perl -w
use strict;

# converts an input VCF file, with might be GZIP compressed, into a  
# multiple sequence alignment in several supported formats 
# (check @validformatsbelow)

# Input: VCF file with reads mapped to concatenated Bd,Bstacei & Bsylvaticum references, might be gzipped 

# Bruno Contreras, Ruben Sancho EEAD-CSIC 2017

my %vcf_real_names = (
'arb_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.arbuscula',
'B422_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.phoenicoides_B422',
'hyb_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.hybridum_BdTR6g',
'mex_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.mexicanum',
'pho_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.phoenicoides',
'pin_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.pinnatum_2x',
'boi_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.boissieri',
'ret_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.retusum',
'rup_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.rupestre',
'sta_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.stacei_TE4.3',
'syl_Cor_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.sylvaticum_Cor',
'syl_Esp_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.sylvaticum_Esp',
'syl_Gre_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.sylvaticum_Gre',
'Bdistachyon_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.distachyon_Bd21',
'Oryza_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'Oryza sativa',
'Sorghum_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'Sorghum bicolor'
);

my %genomic_samples = (
        'syl_Cor_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1,
        'syl_Esp_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1,
        'syl_Gre_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1,
	'Bdistachyon_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1,
	'Oryza_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1,
	'Sorghum_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1
);

# VCF filtering and output options, edit as required
my $MINDEPTHCOVERPERSAMPLE = 10; #10
my $MAXMISSINGSAMPLES      = 8;
my $ONLYPOLYMORPHIC        = 1; # set to 0 to keep fixed loci, helps with sparse data
my $OUTFILEFORMAT          = 'fasta'; 

# VCF format, should be ok but change if needed
my $COLUMNFIRSTSAMPLE      = 9; # zero-based, VCF format http://www.1000genomes.org/node/101, format is previous one
my $GENOTYPECOLUMNFORMAT   = 0; # zero-based, initial guess, they are set for eahc line
my $DEPTHCOLUMNFORMAT      = 1; 

# external binaries, edit if not installed elsewhere and not in path
my $GZIPEXE  = 'gzip'; 
my $BZIP2EXE = 'bzip2';

if(!$ARGV[1]){ die "# usage: $0 <infile.vcf[.gz|.bz2]> <outfile>\n" }

my ($filename,$outfilename) = @ARGV;

my ($n_of_samples,$n_of_loci,$depthcover,$missing,$genotype,$allele) = (0,0);
my ($corr_coord,$sample,$lastsample,$idx,$lastsampleidx,$file);
my (@samplenames,@MSA,%MSAref,%stats,%refallele,%refstrand);
my ($snpname,$badSNP,$shortname,$magic,%contigstats);

my %IUPACdegen = (  'AG'=>'R', 'GA'=>'R', 'CT'=>'Y', 'TC'=>'Y',
                    'CG'=>'S', 'GC'=>'S', 'AT'=>'W', 'TA'=>'W',
                    'GT'=>'K', 'TG'=>'K', 'AC'=>'M', 'CA'=>'M' );	

my %revcomp = ('A'=>'T', 'T'=>'A','G'=>'C','C'=>'G', 'N'=>'N');

my @validformats = qw( phylip nexus fasta );
										
#######################################################

if(!grep(/^$OUTFILEFORMAT$/,@validformats))
{
	die "# unsupported output format, please select one of: ".join(', ',@validformats)."\n";
}

# check input file and choose right way to parse input lines
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
	
printf(STDERR "# input VCF file: $filename\n");
printf(STDERR "# MINDEPTHCOVERPERSAMPLE=$MINDEPTHCOVERPERSAMPLE\n");
printf(STDERR "# MAXMISSINGSAMPLES=$MAXMISSINGSAMPLES\n");
printf(STDERR "# ONLYPOLYMORPHIC=$ONLYPOLYMORPHIC\n");
printf(STDERR "# OUTFILEFORMAT=$OUTFILEFORMAT\n\n");
	
while(<VCF>)   
{  
	chomp($_);
	my @rawdata = split(/\t/,$_);

	if($n_of_samples > 0)
	{
		#Bd1	346	.	A	C	999	PASS	DP=4008;VDB=5.239035e-01;...	GT:PL:DP:SP:GQ	1/1:255,255,0:185:0:99	...

		# skip non-polymorphic sites if required
		next if($rawdata[4] eq '.' && $ONLYPOLYMORPHIC);

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
		my %nts = ( $rawdata[3] => 1 ); # adds ref base call
		$sample=$depthcover=$missing=$badSNP=0;
		($genomic_samples,$gbs_samples) = (0,0);
		foreach $idx ( $COLUMNFIRSTSAMPLE .. $lastsampleidx )
		{
			#0/0:0,255,255:93:0:99
			#1/1:255,255,0:185:0:99
			#0/1:236,0,237:66:7:9
			my @sampledata = split(/:/,$rawdata[$idx]); #print "$rawdata[$idx]\n";
			$genotype = $sampledata[$GENOTYPECOLUMNFORMAT];
			$depthcover = $sampledata[$DEPTHCOLUMNFORMAT]; 
			$allele = 'N'; # by default is unknown
      if($genotype eq '0/0')
      {
        if($depthcover >= $MINDEPTHCOVERPERSAMPLE){ 
				  $allele = $rawdata[3]; 
				  if(defined($genomic_samples{$samplenames[$sample]})){ $genomic_samples++ }
					else{ $gbs_samples++ }
				}
        else{ if(!$genomic_samples{$samplenames[$sample]}){ $missing++ } }
      }
      elsif($genotype eq '1/1')
      {
        if($depthcover >= $MINDEPTHCOVERPERSAMPLE){ 
					$allele = (split(/,/,$rawdata[4]))[0]; 
					if(defined($genomic_samples{$samplenames[$sample]})){ $genomic_samples++ }
          else{ $gbs_samples++ }
				}
        else{ if(!$genomic_samples{$samplenames[$sample]}){ $missing++ } } 
      }
      elsif($genotype eq '2/2')
      {
        if($depthcover >= $MINDEPTHCOVERPERSAMPLE){ 
					$allele = (split(/,/,$rawdata[4]))[1]; 
					if(defined($genomic_samples{$samplenames[$sample]})){ $genomic_samples++ }
          else{ $gbs_samples++ }
				}
        else{ if(!$genomic_samples{$samplenames[$sample]}){ $missing++ } }
      }
      elsif($genotype eq '3/3')
      {
        if($depthcover >= $MINDEPTHCOVERPERSAMPLE){ 		
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
	
			if($missing > $MAXMISSINGSAMPLES)
			{
				$badSNP = 1;
				last;
			}	
			
			$sequence[$sample] = $allele;
			if($allele ne 'N'){ $nts{$allele}++ }

			$sample++;
		} 

		# poor sites are skipped
		if($gbs_samples == 0) 
		{
			$badSNP = 1;
		} 

		# make sure monomorphic sites are skipped
		if(scalar(keys(%nts)) < 2 && $ONLYPOLYMORPHIC){ $badSNP = 1 }
		
		if(!$badSNP)
		{
			if($sample != $n_of_samples)
			{
				die "# VCF line contains less samples than expected ($sample < $n_of_samples):\n$_\n";
			}
		
			$snpname = "$rawdata[0]_$rawdata[1]";
			
			printf(STDERR "# valid locus: $snpname $missing\n");
			
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

			$n_of_loci++;
		}	
	}
	elsif($rawdata[0] eq '#CHROM')
	{
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sampleN
		push(@samplenames,@rawdata[9 .. $#rawdata]);
		$n_of_samples = scalar(@samplenames);
		$lastsample = $n_of_samples-1;
		$lastsampleidx = $#rawdata;
		printf(STDERR "# number of samples found=$n_of_samples\n");
	}
}  
close(VCF);

printf(STDERR "# number of valid loci=$n_of_loci\n");

open(OUTFILE,">$outfilename") || die "# cannot create output file $outfilename, exit\n";


if($OUTFILEFORMAT eq 'nexus')
{
	print OUTFILE "#NEXUS\nBEGIN DATA;\nDIMENSIONS  NTAX=$n_of_samples NCHAR=$n_of_loci;\n";
	print OUTFILE "FORMAT DATATYPE=DNA  MISSING=N GAP=-;\nMATRIX\n";
}
elsif($OUTFILEFORMAT eq 'phylip')
{
	print OUTFILE "$n_of_samples    $n_of_loci\n";		
}

foreach $sample (0 .. $lastsample)
{
	if($OUTFILEFORMAT eq 'nexus')
	{
		print OUTFILE "$vcf_real_names{$samplenames[$sample]} $MSA[$sample]\n";
	}
	elsif($OUTFILEFORMAT eq 'phylip')
	{
		$shortname = $vcf_real_names{$samplenames[$sample]};
		if(length($shortname)>10)
		{ 
			$shortname = '_'.substr($samplenames[$sample],-9); 
			printf(STDERR "# phylip sample name shortened: $samplenames[$sample] -> $shortname\n");
		}
		print OUTFILE "$shortname    $MSA[$sample]\n";
	}
	elsif($OUTFILEFORMAT eq 'fasta')
	{
		print OUTFILE ">$vcf_real_names{$samplenames[$sample]}\n";
		print OUTFILE "$MSA[$sample]\n";
	}
} 

if($OUTFILEFORMAT eq 'nexus')
{
	print OUTFILE ";\nEND;\n";
}

close(OUTFILE);

printf(STDERR "# stats (#SNPs):\n");

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


