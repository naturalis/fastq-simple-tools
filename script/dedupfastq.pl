#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# let's be clear about this
warn "Note: this script assumes single-line FASTQ (as produced by the IonTorrent)\n";

# process command line arguments
my ( $infile, $help );
GetOptions(
	'infile=s'         => \$infile,
	'help+'            => \$help,
);

# give usage message
if ( $help ) {
	die "Usage: $0 --infile=<FASTQ>\n";
}

# this is a simple counter that should allow us to keep track of where
# we are in the FASTQ record
my $line = 0;

# this will be progressively populated to finally create a single 
# FASTQ record, which is printed out if it matches the primer, or
# flushed if it doesn't.
my @record;

# start iterating over file or STDIN
my $fh;
if ( $infile ) {
	open $fh, '<', $infile or die $!;
}
else {
	$fh = \*STDIN;
}

# keep track of how many duplicates 
my %seen;

# flag to see if raw seq is dup
my $isdup;

while(<$fh>) {
	chomp;

	# this is the ID line
	if ( $line % 4 == 0 ) {
		if ( /^\@/ ) {
			@record  = ( $_ );
			$isdup = 0;
		}
		else {
			die "ID line should begin with \@, not $_";
		}
	}

	# this is the sequence line
	elsif ( $line % 4 == 1 ) {
		$isdup = $seen{$_}++;
		push @record, $_;
	}

	# this is the + line
	elsif ( $line % 4 == 2 ) {
		if ( /^\+(?:$record[0])?$/ ) {
			push @record, $_;
		}
		else {
			die "separator should be + (and optional ID), not $_";
		}
	}

	# this is the phred line
	elsif ( $line % 4 == 3 ) {
		push @record, $_;
		print join("\n", @record), "\n" if not $isdup;
	}
	$line++;
}