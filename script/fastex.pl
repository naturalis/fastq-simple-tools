#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# process command line arguments
my ( $fasta, $coord );
GetOptions(
	'fasta=s' => \$fasta,
	'coord=s' => \$coord,
);

# parse the coordinate
my ( $id, $start, $stop );
if ( $coord =~ /^(.+?):(\d+)-(\d+)$/ ) {
	( $id, $start, $stop ) = ( $1, $2, $3 );
}
else {
	die "Couldn't parse coordinate: $coord";
}

# open the file handle
open my $fh, '<', $fasta or die $!;

# iterate over lines
my $seq;
my $flag;
LINE: while(<$fh>) {
	chomp;
	if ( /^>(.+)$/ ) {
		my $current = $1;
		if ( $current eq $id ) {
			$flag = 1;
		}
		else {
			$flag = 0;
			last LINE if $seq;
		}
	}
	else {
		$seq .= $_ if $flag;
	}
}

# extract substring
my ($absstart) = sort { $a <=> $b } $start, $stop;
my $revcom     = $stop < $start;
my $length     = abs( $stop - $start );
my $substring  = substr $seq, $absstart - 1, $length;

# reverse complement, if needed
if ( $revcom ) {
	$substring = reverse $substring;
	$substring =~ tr/ACGT/TGCA/;
}

# print output
print '>', $id, "\n", $substring, "\n";