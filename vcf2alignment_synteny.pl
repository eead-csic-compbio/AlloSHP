#!/usr/bin/perl -w
use strict;

# takes an input VCF file with reads mapped to three concatenated brachy references (Bd,Bstacei & Bsylvaticum),
# with might be GZIP compressed, and uses synteny-based equivalent coordinates to map Bstacei & Bsylvaticum
# back to Bd chromosome positions, in an effort to separate subgenomes, which are multiply aligned in FASTA format

# Input1: VCF file with reads mapped to concatenated Bd,Bstacei & Bsylvaticum references, might be gzipped 
# Input2: SNPs from genomic alignments are also added to the final MSA 

# Bruno Contreras, Ruben Sancho EEAD-CSIC 31May2017

my $SYNTENYZEROBASED = 1; # set to 1 if synteny coords are 0-based

# synteny-based equivalent positions
# these files are used to translate Bstacei & Bsylv positions to Bd coordinates
my %synfiles = (
'Bstacei'=>'Bdistachyon.Bstacei.coords.SNP.tsv',
'Bsylvaticum'=>'Bdistachyon.Bsylvaticum2.coords.SNP.tsv'
);

my %chrcodes = (
'Bdis' =>qr/Bd\d+/,
'Bsta' =>qr/Chr\d+/,
'Bsyl' =>qr/chr\d+/
);

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
'Oryza_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'Oryza_sativa',
'Sorghum_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'Sorghum_bicolor'
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
my ($coord,$subgenome,$code,$corr_snpname,$chr);
my (%MSA,%corr2snpname,%synmap_fwd,%synmap_rev);

my %IUPACdegen = ( 'AG'=>'R', 'GA'=>'R', 'CT'=>'Y', 'TC'=>'Y',
					  'CG'=>'S', 'GC'=>'S', 'AT'=>'W', 'TA'=>'W',
					  'GT'=>'K', 'TG'=>'K', 'AC'=>'M', 'CA'=>'M' );	# by default only biallellic heterozygotes 

my %revcomp = ('A'=>'T', 'T'=>'A','G'=>'C','C'=>'G', 'N'=>'N');

#my @refVCFnames = qw( ref_distachyon ref_Bstacei ref_Bsylvaticum );
#my %refVCFchrcodes = (
#'ref_distachyon' =>qr/Bd\d+/,
#'ref_Bstacei'    =>qr/Chr\d+/,
#'ref_Bsylvaticum'=>qr/chr\d+/
#);
										
my @validformats = qw( phylip nexus fasta );
										
#######################################################

if(!grep(/^$OUTFILEFORMAT$/,@validformats))
{
	die "# not valid output format, please set it one of: ".join(', ',@validformats)."\n";
}

# read syntenic mapped position
warn "\n# SYNTENYZEROBASED=$SYNTENYZEROBASED\n";
foreach my $file (keys(%synfiles))
{
   warn "# parsing $synfiles{$file}...\n";
   my $n_of_positions = 0;
   open(SYNFILE,$synfiles{$file}) || die "# cannot read $synfiles{$file}\n";
   while(<SYNFILE>)
   {
      #Bd1  27606666  A forward Chr01 2828357 A 6 
      my @rawdata = split(/\t/,$_);

      $corr_coord = $rawdata[1];
      $coord      = $rawdata[5];

      if($SYNTENYZEROBASED)
      {
         $corr_coord++;
         $coord++;
      }

      $snpname = "$rawdata[4]_$coord"; # -> position in Bstacei or Bsylvaticum

      if($rawdata[3] eq 'forward'){ $synmap_fwd{$snpname} = "$rawdata[0]_$corr_coord" }
      else{ $synmap_rev{$snpname} = "$rawdata[0]_$corr_coord" }

     $n_of_positions++;
   }
   close(SYNFILE);
   warn "# total positions=$n_of_positions\n";
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

			# check to which subgenome this SNP belongs
                        $subgenome = 'NA';
                        foreach $code (keys(%chrcodes))
                        {
                                if($rawdata[0] =~ /$chrcodes{$code}/)
                                {
                                        $subgenome = $code;
                                        last;
                                }
                        }

		
			#print "$_\n";
	
			$snpname = "$rawdata[0]_$rawdata[1]";
			$corr_snpname = $snpname;

                        if($subgenome ne 'Bdis')
                        {
				# translate coords to common Bd frame using syntenic positions from CGaln
                                if($synmap_fwd{$snpname}){ $corr_snpname = $synmap_fwd{$snpname} }
                                elsif($synmap_rev{$snpname})
                                {
                                        $corr_snpname = $synmap_rev{$snpname};
                                        foreach $sample (0 .. $lastsample)
                                        {
                                                $sequence[$sample] = $revcomp{$sequence[$sample]}
                                        }
                                }
		                else{ next }
			}

			#printf(STDERR "# valid locus: $snpname $corr_snpname $missing\n");

			$chr = (split(/_/,$corr_snpname))[0];
                        foreach $sample (0 .. $lastsample)
                        {
                                $MSA{$corr_snpname}{$subgenome}[$sample] .= $sequence[$sample];

                                # save stats of missing data
                                if($sequence[$sample] eq 'N')
                                {
                                        $contigstats{$subgenome}{$chr}{$sample}{'N'}++;
                                }
                                else
                                {
                                        $contigstats{$subgenome}{$chr}{$sample}{'SNP'}++;
                                }
                        }

                        if($corr_snpname ne $snpname)
                        {
                                $corr2snpname{$corr_snpname} .= "$snpname,";
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

# print some stats
print STDERR "species\tsubgenome";
foreach $subgenome (keys(%contigstats))
{
        foreach $chr (sort keys(%{$contigstats{$subgenome}})){ printf(STDERR "\t%s",$chr) }
        print STDERR "\n";
        last;
}

foreach $sample (0 .. $lastsample)
{
        foreach $subgenome (keys(%contigstats))
        {
                printf(STDERR "%s\t%s",$vcf_real_names{$samplenames[$sample]},$subgenome);
                foreach $chr (sort keys(%{$contigstats{$subgenome}}))
                {
			printf(STDERR"\t%d",$contigstats{$subgenome}{$chr}{$sample}{'SNP'} || 0);
                        #printf(STDERR"\t%d(%s)",$contigstats{$subgenome}{$chr}{$sample}{'SNP'} || 0,$chr);
                }
                print STDERR "\n";
        }
}
print STDERR "\n\n";

# print sorted valid loci
my @sorted_snpnames =
        map  { $_->[0] }
   sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] }
   map  { [$_ , /^Bd(\d+)/ , /_(\d+)/ ] } keys(%MSA);

foreach $snpname (@sorted_snpnames)
{
        printf(STDERR "# aligned position: %s : %s\n",
                $snpname,$corr2snpname{$snpname}||'');
}

printf(STDERR "# number of valid loci=$n_of_loci\n");

# print sorted VCF SNPs
open(OUTFILE,">$outfilename") || die "# cannot create output file $outfilename, exit\n";

foreach $sample (0 .. $lastsample)
{
        foreach $subgenome (sort keys(%chrcodes))
        {
                my ($total,$noNs) = (0,0);

                print OUTFILE ">$vcf_real_names{$samplenames[$sample]}_$subgenome\n";

                foreach $snpname (@sorted_snpnames)
                {
                        $allele = $MSA{$snpname}{$subgenome}[$sample] || 'N';
                        if($allele ne 'N'){ $noNs++ }

                        print OUTFILE $allele;
                        $total++;
                }
                print OUTFILE "\n";

                printf(STDERR "# %s_%s variants: %d / %d (%1.3f)\n",
                        $vcf_real_names{$samplenames[$sample]},$subgenome,
                        $noNs,$total,$noNs/$total);
        }
}
warn "\n";

close(OUTFILE);



