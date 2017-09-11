#!/usr/bin/perl -w
use strict;

# converts an input VCF file, with might be GZIP compressed, into a  
# multiple sequence alignment in several supported formats 
# (check @validformatsbelow)

# Input1: VCF file with reads mapped to concatenated Bd,Bstacei & Bsylvaticum references, might be gzipped 
# Input2: SNPs from genomic alignments are also added to the final MSA 

# Bruno Contreras, Ruben Sancho EEAD-CSIC 2017

my $SYNTENYZEROBASED = 1; # set to 1 if synteny coords are 0-based

# alignments of reference outgroups genomes
# these are used to add outgroups to alignment
#my %ref_files = (
#'ref_RiceNipponBare'=>'/home/contrera/brachy/sintenia/Bdistachyon.rice.ref.coords.tsv',
#'ref_Sorgumbicolor'=>'/home/contrera/brachy/sintenia/Bdistachyon.sorgum.ref.coords.tsv'
#);

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
my $OUTFILEFORMAT          = 'fasta'; 

# VCF format, should be ok but change if needed
my $COLUMNFIRSTSAMPLE      = 9; # zero-based, VCF format http://www.1000genomes.org/node/101, format is previous one

my $GENOTYPECOLUMNFORMAT   = 0; # zero-based, initial guess, they are set for eahc line
my $DEPTHCOLUMNFORMAT      = 1; 

# external binaries
my $GZIPEXE = 'gzip'; # edit if not installed elsewhere and not in path

if(!$ARGV[1]){ die "# usage: $0 <infile.vcf[.gz]> <outfile>\n" }

my ($filename,$outfilename) = @ARGV;

my ($n_of_samples,$n_of_loci,$depthcover,$missing,$genotype,$allele) = (0,0);
my ($corr_coord,$sample,$lastsample,$idx,$lastsampleidx,$file);
my (@samplenames,@MSA,%MSAref,%stats,%refallele,%refstrand);
my ($snpname,$badSNP,$shortname,$magic,%contigstats);

my %IUPACdegen = ( 'AG'=>'R', 'GA'=>'R', 'CT'=>'Y', 'TC'=>'Y',
					  'CG'=>'S', 'GC'=>'S', 'AT'=>'W', 'TA'=>'W',
					  'GT'=>'K', 'TG'=>'K', 'AC'=>'M', 'CA'=>'M' );	# by default only biallellic heterozygotes 

my %revcomp = ('A'=>'T', 'T'=>'A','G'=>'C','C'=>'G', 'N'=>'N');

#my @refVCFnames = qw( ref_distachyon ref_Bstacei ref_Bsylvaticum );
my %refVCFchrcodes = (
'ref_distachyon' =>qr/Bd\d+/,
'ref_Bstacei'    =>qr/Chr\d+/,
'ref_Bsylvaticum'=>qr/chr\d+/
);
										
my @validformats = qw( phylip nexus fasta );
										
#######################################################

if(!grep(/^$OUTFILEFORMAT$/,@validformats))
{
	die "# not valid output format, please set it one of: ".join(', ',@validformats)."\n";
}

# read alignments of reference genomes
#foreach $file (sort keys(%ref_files))
#{
#   warn "# reading reference $ref_files{$file}...\n";
#   open(REF,$ref_files{$file}) || die "# cannot read $ref_files{$file}\n";
#   while(<REF>)
#   {
#      #Bd3 3725899 G reverse chr4  32091475  G 64
#      chomp;
#      my @rawdata = split(/\t/,$_);
#
#      $corr_coord = $rawdata[1];
#      if($SYNTENYZEROBASED)
#      {
#         $corr_coord++;
#      }
#
#      $snpname = "$rawdata[0]_$corr_coord"; # -> position in Bdistachyon
#
#      if($rawdata[3] eq 'forward')
#      {
#         $refallele{$snpname}{$file} = $rawdata[6];
#         $refstrand{$snpname}{$file} = 1;
#      }
#      else
#      {
#         $refallele{$snpname}{$file} = $revcomp{$rawdata[6]};
#         $refstrand{$snpname}{$file} = -1;
#      }
#		
#      #print "$snpname $file $refallele{$snpname}{$file}\n";
#   }
#   close(REF);
#}

# check input file and choose right way to parse input lines
my ($genomic_samples,$gbs_samples); 
open(INFILE,$filename) || die "# cannot open input file $filename, exit\n";
sysread(INFILE,$magic,2);
close(INFILE);

if($magic eq "\x1f\x8b") # compressed input
{
	printf(STDERR "# decompressing VCF file with GZIP\n");
	open(VCF,"$GZIPEXE -dc $filename |");
} 
else{ open(VCF,$filename) }
	
printf(STDERR "# input VCF file: $filename\n");
printf(STDERR "# MINDEPTHCOVERPERSAMPLE=$MINDEPTHCOVERPERSAMPLE\n");
printf(STDERR "# MAXMISSINGSAMPLES=$MAXMISSINGSAMPLES\n");
#printf(STDERR "# GENOTYPECOLUMNFORMAT=$GENOTYPECOLUMNFORMAT\n");
#printf(STDERR "# DEPTHCOLUMNFORMAT=$DEPTHCOLUMNFORMAT\n");
printf(STDERR "# OUTFILEFORMAT=$OUTFILEFORMAT\n");
	
while(<VCF>)   
{  
	chomp($_);
	my @rawdata = split(/\t/,$_);

	if($n_of_samples > 0)
	{
		#Bd1	346	.	A	C	999	PASS	DP=4008;VDB=5.239035e-01;...	GT:PL:DP:SP:GQ	1/1:255,255,0:185:0:99	...

		#next if($rawdata[0] ne /^Bd1/); # debugging
		#next if($rawdata[0] ne /^Chr01/); # debugging
	
		# skip non-chr contigs
		#next if($rawdata[0] =~ m/scaffold/ || $rawdata[0] =~ m/entromer/ || $rawdata[0] =~m /Segkk/);

		# skip non-polymorphic sites
		#next if($rawdata[4] eq '.');

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
			$depthcover = $sampledata[$DEPTHCOLUMNFORMAT]; #die "$genotype $depthcover\n";
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

		# make sure genomic-only sites are skipped
		if($gbs_samples == 0) # && $genomic_samples > 0)
		{
			$badSNP = 1;
			#warn "# VCF line contains only genomic samples:\n$_\n";
		} #else { warn"# $gbs_samples + $missing + $genomic_samples ($badSNP)\n" }

		#make sure monomorphic sites are skipped
		#if(scalar(keys(%nts)) < 2){ $badSNP = 1 }
		
		if(!$badSNP)
		{
			if($sample != $n_of_samples)
			{
				die "# VCF line contains less samples than expected ($sample < $n_of_samples):\n$_\n";
			}
		
			#print "$_\n";
	
			$snpname = "$rawdata[0]_$rawdata[1]";
			
			printf(STDERR "# valid locus: $snpname $missing\n");
			
			foreach $sample (0 .. $lastsample)
			{
				$MSA[$sample] .= $sequence[$sample];
				#print "$sequence[$sample]\n";

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

			# add SNPs from reference genomic alignments (outgroups)
			#foreach $file (sort keys(%ref_files))
			#{
			#	if($refallele{$snpname}{$file})
			#	{
			#		$allele = $refallele{$snpname}{$file};
			#		#$refstrand{$snpname}{$file} = -1;				
			#		#$revcomp{$rawdata[6]};
			#	}
			#	else
			#	{ 
			#		$allele = 'N';
			#	}

			#	$MSAref{$file} .= $allele;
			#}

			# add nucleotides from reference sequence in VCF file
			#foreach $file (@refVCFnames)
			#{
			#	if($rawdata[0] =~ m/$refVCFchrcodes{$file}/)
			#	{
			#		$allele = $rawdata[3];
			#	}
			#	else{ $allele = 'N' }
			#	
			#	$MSAref{$file} .= $allele;
			#}

			#last; #debug
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

# outgroup references
#foreach $file (sort keys(%ref_files))
#{
#	if($OUTFILEFORMAT eq 'nexus')
#   {
#      print OUTFILE "$file $MSAref{$file}\n";
#   }
#   elsif($OUTFILEFORMAT eq 'phylip')
#   {
#      $shortname = $file;
#      if(length($shortname)>10)
#      {
#         $shortname = '_'.substr($samplenames[$sample],-9);
#         printf(STDERR "# phylip sample name shortened: $samplenames[$sample] -> $shortname\n");
#      }
#      print OUTFILE "$shortname    $MSAref{$file}\n";
#   }
#   elsif($OUTFILEFORMAT eq 'fasta')
#   {
#      print OUTFILE ">$file\n";
#      print OUTFILE "$MSAref{$file}\n";
#   }
#}

# VCF references
#foreach $file (@refVCFnames)
#{  
#   if($OUTFILEFORMAT eq 'nexus')
#   {  
#      print OUTFILE "$file $MSAref{$file}\n";
#   }
#   elsif($OUTFILEFORMAT eq 'phylip')
#   {
#      $shortname = $file;
#      if(length($shortname)>10)
#      {  
#         $shortname = '_'.substr($samplenames[$sample],-9);
#         printf(STDERR "# phylip sample name shortened: $samplenames[$sample] -> $shortname\n");
#      }
#      print OUTFILE "$shortname    $MSAref{$file}\n";
#   }
#   elsif($OUTFILEFORMAT eq 'fasta')
#   {  
#      print OUTFILE ">$file\n";
#      print OUTFILE "$MSAref{$file}\n";
#   }
#}

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


