#!/usr/bin/env perl
use strict;
use warnings;

# parses a FASTA file obtained with Cgaln and produces
# the corresponding dot file
#
# B Contreras-Moreira, R Sancho EEAD-CSIC 2024

my ($fastafile,$outdotfile);

if(!$ARGV[1]) { 
  die "# usage: $0 <infile.fasta> <outfile.dot>\n" 

}else { 
  ($fastafile,$outdotfile) = @ARGV 
}

warn "# infasta: $fastafile\n";
warn "# dot file: $outdotfile\n";

my ($hsp,$length,$cumulscore);
my ($startA,$endA,$startB,$endB);

## 0) create outfile
open(DOT,">",$outdotfile) || die "# ERROR: cannot create $outdotfile\n";

## 1) parse FASTA file
#####{ (forward) genomeA-fasta1	genomeB-fasta1	Bl #1
#>A_fst1:58679688-58679721:HSP number 105:score 68:score_cumulative 12390
#CAAAAAGATGGCTCTGCGCCAAGACGGTGGAGCT
#>B_fst1:389883-389916:HSP number 105:score 68:score_cumulative 12390
#CAAGAAAACGACTCTGCGCTAAGACGGTGGAGCT
#####}

open(FASTA,"$fastafile") || die "# ERROR: cannot read $fastafile\n";
while(<FASTA>) {
  if(/^####/){ 
    print DOT; 

  } elsif(/^>A_fst\d+[^:]*:(\d+)-(\d+):HSP number (\d+):score (\d+):score_cumulative (\d+)/) {

    ($startA,$endA,$hsp,$length,$cumulscore) = ($1,$2,$3,$4,$5);
    printf(DOT "#HSP number: %d, length: %d, score: %d, score_cumulative: %d\n",
      $hsp, $length, $length, $cumulscore);

  } elsif(/^>B_fst\d+:(\d+)-(\d+)/) { 
    ($startB,$endB) = ($1,$2);
    printf(DOT "%d\t%d\n%d\t%d\n\n",
      $startA,$startB,$endA,$endB);

  } else { }
}
close(FASTA);

close(DOT);
