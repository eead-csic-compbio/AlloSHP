#!/usr/bin/env perl
use strict;
use warnings;

# parses a FASTA file obtained with Cgaln and produces
# the corresponding dot file
#
# B Contreras-Moreira, R Sancho EEAD-CSIC 2024

my ($fastafile,$lenfile,$outdotfile);

if(!$ARGV[2]) { 
  die "# usage: $0 <infile.fasta> <lengths.tsv> <outfile.dot>\n" 

} else { 
  ($fastafile,$lenfile,$outdotfile) = @ARGV 
}

warn "# infasta: $fastafile\n";
warn "# lenfile: $lenfile\n";
warn "# dotfile: $outdotfile\n";

my ($chr,$hsp,$score,$length,$cumulscore,$isrev);
my ($startA,$endA,$startB,$endB,$n_of_chrs);
my (%len,%totalen,%chrs);

## 1) parse chr lengths
open(LEN,"<",$lenfile) || die "# ERROR: cannot read $lenfile\n";
while(<LEN>) {
  #A	Bd1	75071545
  chomp;
  my @data = split(/\t/,$_);

  if(!defined($len{$data[0]})){
    $n_of_chrs = 1;
    $totalen{$data[0]} = 0;
    $chrs{$data[0]} = 0;
  } else {
    $n_of_chrs++;
  }

  $len{$data[0]}{$n_of_chrs} = $data[2];
  $totalen{$data[0]} += $data[2];
  $chrs{$data[0]}++;   
}
close(LEN); 


## 2) parse FASTA and write to DOT file
#####{ (forward) genomeA-fasta1	genomeB-fasta1	Bl #1
#>A_fst1:58679688-58679721:HSP number 105:score 68:score_cumulative 12390
#CAAAAAGATGGCTCTGCGCCAAGACGGTGGAGCT
#>B_fst1:389883-389916:HSP number 105:score 68:score_cumulative 12390
#CAAGAAAACGACTCTGCGCTAAGACGGTGGAGCT
#####}

open(DOT,">",$outdotfile) || die "# ERROR: cannot create $outdotfile\n";

open(FASTA,"$fastafile") || die "# ERROR: cannot read $fastafile\n";
while(<FASTA>) {
  if(/^####/){ 
    print DOT; 

    if(/genomeA-fasta(\d+)/) {
      $n_of_chrs = $1;      
    }

    if(/\(reverse/) { $isrev = '(revcom)' }
    else{ $isrev = '' }  

  } elsif(/^>A_fst\d+[^:]*:(\d+)-(\d+):HSP number (\d+):score (\d+):score_cumulative (\d+)/) {

    ($startA,$endA,$hsp,$score,$cumulscore) = ($1,$2,$3,$4,$5);
    $length = $score/2; # ~ https://github.com/rnakato/Cgaln/blob/f006f56a88f334056263c58a0641f0fcb27645a6/NA.h#L12

    foreach $chr (1 .. $n_of_chrs-1) {	    
      $startA += $len{'A'}{$chr};
      $endA += $len{'A'}{$chr};      
    }

    printf(DOT "#HSP%s number: %d, length: %d, score: %d, score_cumulative: %d\n",
      $isrev, $hsp, $length, $score, $cumulscore);

  } elsif(/^>B_fst\d+:(\d+)-(\d+)/) { 
    ($startB,$endB) = ($1,$2);
    printf(DOT "%d\t%d\n%d\t%d\n\n",
      $startA,$startB,$endA,$endB);

  } else { }
}
close(FASTA);


## 2.1) print chr sizes for grid

# A
my $cumul = 0;
foreach $chr (1 .. $chrs{'A'}) {
  
  printf(DOT "%d\t%d\n",
    $len{'A'}{$chr}+$cumul, 0);
  printf(DOT "%d\t%d\n\n\n",
    $len{'A'}{$chr}+$cumul, $totalen{'B'});

  $cumul += $len{'A'}{$chr};
}

# B
$cumul = 0;
foreach $chr (1 .. $chrs{'B'}) {

  printf(DOT "%d\t%d\n",
    0, $len{'B'}{$chr}+$cumul);
  printf(DOT "%d\t%d\n\n\n",
    $totalen{'A'}, $len{'B'}{$chr}+$cumul);

  $cumul += $len{'B'}{$chr};
}


close(DOT);



