#! /usr/bin/perl

my (@seq);
my $index1 = -1;
my $index2;
while (my $line = <STDIN>)
{
	chomp $line;
	if ($line =~ /^>/)
	{
		$index1++;
		$index2 = -1;
	}
	else
	{
		my @temp = split (//, $line);
		for (my $i = 0; $i <= $#temp; $i++)
		{
			$index2++;
			$seq[$index1][$index2] = $temp[$i];
		}
	}
}

my $count = 0;
my $diff = 0;
for (my $i = 0; $i <= $index2; $i++)
{
	my $nucl_1 = &maiuscolo ($seq[0][$i]);
	my $nucl_2 = &maiuscolo ($seq[1][$i]);
	if ( ($nucl_1 ne 'N') && ($nucl_1 ne '-') && ($nucl_2 ne 'N') && ($nucl_2 ne '-') )
	{
		$count++;
		if ($nucl_1 ne $nucl_2) {$diff++;}
	}
}
my $value = sprintf("%.6f", 1 - ($diff / $count) );
print STDOUT $value."\n";

sub maiuscolo
{
	my ($seq);
	($seq) = @_;
	my $seq_m = $seq;
	if ( ($seq eq 'A') || ($seq eq 'a') ) {$seq_m = 'A';}
	elsif ( ($seq eq 'C') || ($seq eq 'c') ) {$seq_m = 'C';}
	elsif ( ($seq eq 'G') || ($seq eq 'g') ) {$seq_m = 'G';}
	elsif ( ($seq eq 'T') || ($seq eq 't') ) {$seq_m = 'T';}
	elsif ( ($seq eq 'N') || ($seq eq 'n') ) {$seq_m = 'N';}
	elsif ( ($seq eq 'R') || ($seq eq 'r') ) {$seq_m = 'R';}
	elsif ( ($seq eq 'Y') || ($seq eq 'y') ) {$seq_m = 'Y';}
	elsif ( ($seq eq 'S') || ($seq eq 's') ) {$seq_m = 'S';}
	elsif ( ($seq eq 'W') || ($seq eq 'w') ) {$seq_m = 'W';}
	elsif ( ($seq eq 'K') || ($seq eq 'k') ) {$seq_m = 'K';}
	elsif ( ($seq eq 'M') || ($seq eq 'm') ) {$seq_m = 'M';}
	elsif ( ($seq eq 'B') || ($seq eq 'b') ) {$seq_m = 'B';}
	elsif ( ($seq eq 'D') || ($seq eq 'd') ) {$seq_m = 'D';}
	elsif ( ($seq eq 'H') || ($seq eq 'h') ) {$seq_m = 'H';}
	elsif ( ($seq eq 'V') || ($seq eq 'v') ) {$seq_m = 'V';}
	return ($seq_m);
}
