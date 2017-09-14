#!/usr/bin/perl -w
use strict;

# Takes an input FASTA aligned with vcf2alignment_synteny and collapses subgenome lines (Bdis, Bsta and Bsyl)
# in one line 

# Bruno Contreras, Ruben Sancho EEAD-CSIC June2017

if(!$ARGV[1]){ die "# usage: $0 <FASTA file with 1 subgenome/sequence> <FASTA file with single consensus sequence>\n" }

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

my ($seq,$pos,$base,@sequences);
my $sequence_length = 0;
my $n_of_sequences = 0;

# check input file and choose right way to parse input lines
open(FASTA,$infile) || die "# cannot open input file $infile, exit\n";
while(my $line = <FASTA>)
{
	next if($line =~ m/^>/);
	# por tanto solo quedan las secuencias o lineas en blanco

	chomp($line); # elimino salto de linea final
	$seq = uc($line);
	$sequence_length = length($seq);
	$sequences[$n_of_sequences] = [ split(//,$seq) ];

	$n_of_sequences++;
}
close(FASTA);

# abre fichero de salida
open(OUTFASTA,">",$outfile);
print OUTFASTA ">collapsed\n";

# recorre el bloque de sequencias de izq a derecha
foreach $pos (0 .. $sequence_length-1) # en Perl los indices se cuentan desde 0, como en muchos otros lenguajes
{
	# recorre cada columna alineada del bloque y cuenta las bases observadas
	my %bases_alineadas;
	foreach $seq (0 .. $#sequences) # $#sequences es el ultimo indice valido de la lista @sequences
	{
		$base = $sequences[$seq][$pos];
		$bases_alineadas{$base}++;
	}

	# keys(%bases_alineadas) devuelve una lista con las bases, que son las llaves de la tabla %bases_alineadas
	my @bases_columna = keys(%bases_alineadas);
	if(scalar(@bases_columna) == 1) # bases conservadas
	{		
		print OUTFASTA "$bases_columna[0]";
	}
	elsif(scalar(@bases_columna) == 2 && defined($bases_alineadas{'N'})) # una base conservada + Ns
        {
		foreach $base (@bases_columna)
		{            
			next if($base eq 'N');   
                	print OUTFASTA "$base";
		}
        }
	else
	{
		print OUTFASTA "N";
	}
}

print OUTFASTA "\n"; # salto de linea

close(OUTFASTA);	






