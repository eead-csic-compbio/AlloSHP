#!/usr/bin/env perl
use strict;
use warnings;

# parses GSAlign whole-genome alignments of A & B species and 
# produces:
# i) a FASTA file matching the format of Cgaln
# ii) a dot file matching the format of Cgaln
#
# B Contreras-Moreira, R Sancho EEAD-CSIC 2024

my $MINHSPLENGTH = 19; # as in Cgaln

my ($alnfile,$fastafileA,$fastafileB,$outfastafile,$outdotfile);

if(!$ARGV[2]) { 
  die "# usage: $0 <GSalign file.aln> <infileA.fasta> <infileB.fasta> <outfile.fasta> <outfile.dot>\n" 

}else { 
  ($alnfile,$fastafileA,$fastafileB,$outfastafile,$outdotfile) = @ARGV 
}

warn "# aln file: $alnfile\n";
warn "# infastaA: $fastafileA\n";
warn "# infastaB: $fastafileB\n";
warn "# FASTA file: $outfastafile\n";
warn "# dot file: $outdotfile\n";
warn "# \$MINHSPLENGTH = $MINHSPLENGTH\n\n";

my ($block,$hsp,$n_of_chrs,$length) = (0,0,0,0);
my ($strandA,$chrA,$chrB,$cumulscore,$pos,$seq);
my ($startA,$endA,$startB,$endB,$seqA,$seqB);
my ($baseA,$baseB,$coordA,$coordB,$hspA,$hspB);
my ($hsp2print,$dot2print);
my (%chr2numA,%chr2numB,@hsps,@dots);

## 0) create outfiles
open(FASTA,">",$outfastafile) || die "# ERROR: cannot create $outfastafile\n";
open(DOT,">",$outdotfile) || die "# ERROR: cannot create $outdotfile\n";

## 1) read FASTA files to get real chr names
$n_of_chrs = 0;
open(AFASTA,"$fastafileA") || die "# ERROR: cannot read $fastafileA\n";
while(<AFASTA>) {
  if(/^>(\S+)/) {
    $n_of_chrs++; 
    $chr2numA{ $1 } = $n_of_chrs;
  }
}
close(AFASTA);
warn "# number of chrs in A FASTA file: $n_of_chrs\n";

$n_of_chrs = 0;
open(BFASTA,"$fastafileB") || die "# ERROR: cannot read $fastafileB\n";
while(<BFASTA>) {
  if(/^>(\S+)/) {
    $n_of_chrs++; 
    $chr2numB{ $1 } = $n_of_chrs;        
  }
}
close(BFASTA);
warn "# number of chrs in B FASTA file: $n_of_chrs\n";

## 2) parse genome alignment blocks

# blocks start with:   
#Identity = 1324 / 1440 (91.90%) Orientation = Forward

# alignments are consecutive pairs of lines like:
#ref.Bd2          24148001    ttcaagaaggcaaaggGTAAGAAGGGAAATTTCAAGAAGGGTGGCAAA
#qry.Chr01        16470050    ttcaagaaggcaaagg--------------ttcaaaaag---------

# blocks end with:
#**********************

# dot file should look like this;
#####{ (forward) genomeA-fasta1 genomeB-fasta1  Bl #1
#HSP number: 105, length: 34, score: 68, score_cumulative: 12390
#58679688        389883
#58679721        389916

open(ALN,"<",$alnfile) || die "# ERROR: cannot read $alnfile\n";
while(<ALN>) {

  if(/^#Identity = \d+ \/ \d+ \(\S+\) Orientation = (\S+)/) {
    $strandA = $1;

    #initialize block 
    $block++;
    ($seqA,$seqB,$hspA,$hspB,$hsp,$cumulscore) = ('','','','',0,0);
    ($coordA,$coordB) = (-1,-1);    
    @hsps = ();
    @dots = ();

  } elsif(/^ref\.(\S+)\s+(\d+)\s+(\S+)/) {
    ($chrA,$pos,$seq) = ($1,$2,$3);
    $seqA .= $seq;
    if($coordA < 0){ $coordA = $pos }

  } elsif(/^qry\.(\S+)\s+(\d+)\s+(\S+)/) {
    ($chrB,$pos,$seq) = ($1,$2,$3);
    $seqB .= $seq;
    if($coordB < 0){ $coordB = $pos }
    
  } elsif(/^\*\*\*/) { # process HSPs in this block

    if(length($seqA) != length($seqB)) {
      die "# ERROR: failed due to aligned sequences with different lenght\n";
    }

    if($strandA eq 'Forward') {
      printf(FASTA "#####{ (forward) genomeA-fasta%d genomeB-fasta%d  Bl #%d\n",
        $chr2numA{ $chrA },
        $chr2numB{ $chrB },	
        $block); 

      printf(DOT "#####{ (forward) genomeA-fasta%d genomeB-fasta%d  Bl #%d\n",
        $chr2numA{ $chrA },
        $chr2numB{ $chrB },
        $block);

    } elsif($strandA eq 'Reverse') {
      printf(FASTA "#####{ (reverse) genomeA-fasta%d genomeB-fasta%d  Bl #%d\n",
        $chr2numA{ $chrA },
        $chr2numB{ $chrB },
        $block);

      printf(DOT "#####{ (reverse) genomeA-fasta%d genomeB-fasta%d  Bl #%d\n",
        $chr2numA{ $chrA },
        $chr2numB{ $chrB },
        $block);

    } else {
      next; # skip any other cases
    } 

    # split multi-line pairwsie alignment in gapless HSPs,
    # this mimics Cgaln output
    $pos = 0;
    $length = 0;
    ($startA,$startB) = ($coordA, $coordB);
    ($endA,$endB) = ($coordA, $coordB); 
    while($pos <= length($seqA)) {

      $baseA = substr($seqA,$pos,1);
      $baseB = substr($seqB,$pos,1); 

      if($baseA eq '-' || $baseB eq '-' || $pos == length($seqA)) {

	if($length > 0) {      
          $hsp++;

          # correct end coords
	  $endB--;
	  if($strandA eq 'Forward') {
	    $endA--;
	  } else {
	    $endA++;
	  }

	  # process previous HSP when indel found, both FASTA and DOT, provided is long enough
	  if($length > $MINHSPLENGTH) {

            $hsp2print = sprintf(">A_fst%d:%d-%d:HSP number %d:score %d:score_cumulative %d\n",
              $chr2numA{ $chrA },$startA,$endA,$hsp,$length,$cumulscore); 
            $hsp2print .= "$hspA\n";
            $hsp2print .= sprintf(">B_fst%d:%d-%d:HSP number %d:score %d:score_cumulative %d\n",
              $chr2numB{ $chrB },$startB,$endB,$hsp,$length,$cumulscore); 
            $hsp2print .= "$hspB\n"; 
 
            $dot2print = sprintf("#HSP number: %d, length: %d, score: %d, score_cumulative: %d\n",
              $hsp,$length,$length,$cumulscore);
            $dot2print .= sprintf("%d\t%d\n%d\t%d\n",
              $startA,$startB,$endA,$endB); 		 

            push(@hsps, $hsp2print);
            push(@dots, $dot2print); 	    
          }

	  # init next HSP
	  $length = 0;
          ($hspA,$hspB) = ('','');
          if($strandA eq 'Forward') {
            ($startA,$startB) = ($endA+1,$endB+1);
            ($endA,$endB) = ($startA,$startB);
          } else {
            ($startA,$startB) = ($endA-1,$endB+1);
            ($endA,$endB) = ($startA,$startB);		  
          }
	}
         
        if($baseA eq '-') {
          $startB++;
          $endB++;

	} elsif($baseB eq '-') {
          if($strandA eq 'Forward') {
            $startA++;
            $endA++;

	  } else {
            $startA--;
            $endA--;
	  }
	}
      } else { # grow current HSP

        $hspA .= $baseA;
        $hspB .= $baseB;	
        $length++;
	$cumulscore++;
	$endB++;
	if($strandA eq 'Forward') {
          $endA++;
        } else {
          $endA--;
	}
      }	      

      $pos++;	    
    }

    # print HSPs in reverse order, as in Cgaln
    foreach $hsp2print (reverse(@hsps)) {
      print FASTA "$hsp2print\n";
    }
    print FASTA "#####}\n";

    foreach $dot2print (reverse(@dots)) {
      print DOT "$dot2print\n";
    }
    print DOT "#####}\n";
  }
}
close(ALN);

close(DOT);
close(FASTA);

printf(STDERR "# total blocks: %d\n\n", $block);
