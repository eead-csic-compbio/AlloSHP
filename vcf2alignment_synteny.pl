#!/usr/bin/perl -w
use strict;

# Takes an input VCF file with reads mapped to several concatenated references,
# and uses synteny-based equivalent coordinates to map polymorphism back to the 
# main/user-defined chromosome positions, in an effort to separate subgenomes, 
# which are then multiply aligned.
# Only FASTA format is produced to facilitate collapsing sequences afterwards.

# Input: VCF file with reads mapped to concatenated references, might be compressed

# Bruno Contreras, Ruben Sancho EEAD-CSIC 2017

# Synteny-based equivalent positions produced with utils/mapcoords.pl
# These files are used to translate Bstacei & Bsylv positions to Bd coordinates
my %synfiles = (
  'Bstacei'     =>'sample_data/Bdistachyon.Bstacei.coords.SNP.tsv',
  'Bsylvaticum' =>'sample_data/Bdistachyon.Bsylvaticum.coords.SNP.tsv'
);

# Prefixes/regular expressions that connect chromosomes to reference genomes
my %chrcodes = (
  'Bdis' =>qr/Bd\d+/,
  'Bsta' =>qr/Chr\d+/,
  'Bsyl' =>qr/chr\d+/
);

# Please edit and uncomment to shorten sample names in output alignment
my %vcf_real_names = (
  #'arb_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.arbuscula',
  #'hyb_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 'B.hybridum_BdTR6g',
);

# Please edit and uncomment to set samples which should not count as missing data.
# For instance, we used it to leave outgroups out of these calculations, as their
# reads were WGS reads with significantly more depth than GBS/RNAseq samples
my %genomic_samples = (
  'syl_Cor_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1,
  'syl_Esp_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1,
  'syl_Gre_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1,
  'Bdistachyon_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1,
  'Oryza_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1,
  'Sorghum_map_Bd_Bs_Bsyl_no_contigs.sorted.q30.bam' => 1
);

# VCF filtering and output options, edit as required
my $SYNTENYZEROBASED       = 1;  # set to 1 if synteny coords are 0-based, as those produced by CGaln
my $MINDEPTHCOVERPERSAMPLE = 10; # natural, min number of reads mapped supporting a locus
my $MAXMISSINGSAMPLES      = 8;  # natural, max number of missing samples accepted per locus
my $ONLYPOLYMORPHIC        = 0;  # set to 0 to keep fixed loci, helps with sparse data
my $OUTFILEFORMAT          = 'fasta'; # can also take other formats in @validformats



# first guess of key VCF columns, adjusted in real time below
my $COLUMNFIRSTSAMPLE      = 9; # zero-based, VCF format http://www.1000genomes.org/node/101, format is previous one
my $GENOTYPECOLUMNFORMAT   = 0; # zero-based, initial guess, they are set for eahc line
my $DEPTHCOLUMNFORMAT      = 1; 

# external binaries, edit if not installed elsewhere and not in path
my $GZIPEXE  = 'gzip'; 
my $BZIP2EXE = 'bzip2';

if(!$ARGV[1]){ die "# usage: $0 <infile.vcf[.gz|.bz2]> <output alignment file>\n" }

my ($filename,$outfilename) = @ARGV;

my ($n_of_samples,$n_of_loci,$depthcover,$missing,$genotype,$allele) = (0,0);
my ($corr_coord,$sample,$lastsample,$idx,$lastsampleidx,$file);
my (@samplenames,@MSA,%MSAref,%stats,%refallele,%refstrand);
my ($snpname,$badSNP,$shortname,$magic,%contigstats);
my ($coord,$subgenome,$code,$corr_snpname,$chr);
my (%MSA,%corr2snpname,%synmap_fwd,%synmap_rev);

my %IUPACdegen = (  'AG'=>'R', 'GA'=>'R', 'CT'=>'Y', 'TC'=>'Y',
                    'CG'=>'S', 'GC'=>'S', 'AT'=>'W', 'TA'=>'W',
                    'GT'=>'K', 'TG'=>'K', 'AC'=>'M', 'CA'=>'M' );	 

my %revcomp = ('A'=>'T', 'T'=>'A','G'=>'C','C'=>'G', 'N'=>'N');

my @validformats = qw( fasta );
										
#######################################################

if(!grep(/^$OUTFILEFORMAT$/,@validformats))
{
	die "# not valid output format, please set it one of: ".join(', ',@validformats)."\n";
}

# read syntenic mapped position
warn "\n# SYNTENYZEROBASED=$SYNTENYZEROBASED\n";
foreach my $file (keys(%synfiles))
{
   warn "# parsing syntenic positions $synfiles{$file}...\n";
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
warn "\n";


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
printf(STDERR "# OUTFILEFORMAT=$OUTFILEFORMAT\n");
	
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



