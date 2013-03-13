#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use String::Approx 'adist';

# let's be clear about this
warn "Note: this script assumes single-line FASTQ (as produced by the IonTorrent)\n";

# process command line arguments
my $mismatch = 4;
my ( $infile, $help, $primerfile, $with_mid_labels, $test, $labels );
GetOptions(
	'primerfile=s'     => \$primerfile,
	'infile=s'         => \$infile,
	'help+'            => \$help,
	'with_mid_labels+' => \$with_mid_labels,
	'test+'            => \$test,
	'mismatch=i'       => \$mismatch,
);

# ion torrent-specific prefix
my $prefix = 'CCATCTCATCCCTGCGTGTCTCCGACTCAG';

# maybe strip these out?
my %midlabels = (
	IonXpress_001 => 'CTAAGGTAAC',
	IonXpress_002 => 'TAAGGAGAAC',
	IonXpress_003 => 'AAGAGGATTC',
	IonXpress_004 => 'TACCAAGATC',
	IonXpress_005 => 'CAGAAGGAAC',
	IonXpress_006 => 'CTGCAAGTTC',
	IonXpress_007 => 'TTCGTGATTC',
	IonXpress_008 => 'TTCCGATAAC',
	IonXpress_009 => 'TGAGCGGAAC',
	IonXpress_010 => 'CTGACCGAAC',
	IonXpress_011 => 'TCCTCGAATC',
	IonXpress_012 => 'TAGGTGGTTC',
	IonXpress_013 => 'TCTAACGGAC',
	IonXpress_014 => 'TTGGAGTGTC',
	IonXpress_015 => 'TCTAGAGGTC',
	IonXpress_016 => 'TCTGGATGAC',
	IonXpress_017 => 'TCTATTCGTC',
	IonXpress_018 => 'AGGCAATTGC',
	IonXpress_019 => 'TTAGTCGGAC',
	IonXpress_020 => 'CAGATCCATC',
	IonXpress_021 => 'TCGCAATTAC',
	IonXpress_022 => 'TTCGAGACGC',
	IonXpress_023 => 'TGCCACGAAC',
	IonXpress_024 => 'AACCTCATTC',
	IonXpress_025 => 'CCTGAGATAC',
	IonXpress_026 => 'TTACAACCTC',
	IonXpress_027 => 'AACCATCCGC',
	IonXpress_028 => 'ATCCGGAATC',
	IonXpress_029 => 'TCGACCACTC',
	IonXpress_030 => 'CGAGGTTATC',
	IonXpress_031 => 'TCCAAGCTGC',
	IonXpress_032 => 'TCTTACACAC',
);

# give usage message
if ( $help or not $primerfile or not $infile ) {
	die "Usage: $0 --primerfile=<TSV> --infile=<FASTQ>\n";
}

# read the different primers
my ( %primers, %handles );
{
	open my $fh, '<', $primerfile or die $!;
	while(<$fh>) {
		chomp;
		my ( $name, $pattern ) = split /\t/, $_;
		$pattern =~ s/^$prefix//; # IonTorrent specific
		if ( not $with_mid_labels ) {
			for my $label ( keys %midlabels ) {
				if ( $pattern =~ /$midlabels{$label}/ ) {
					warn "stripping $label from $name\n";
					$pattern =~ s/$midlabels{$label}//;
				}
			}
		}
		warn "$name\t$pattern\n";
		$primers{$name} = $pattern;
		my $outfile = $infile;
		$outfile =~ s/.fastq$/-$name.fastq/;
		if ( not $test ) {
			open $handles{$name}, '>>', $outfile or die $!;
		}
	}
}

# this is a simple counter that should allow us to keep track of where
# we are in the FASTQ record
my $line = 0;

# this will be progressively populated to finally create a single 
# FASTQ record, which is printed out if it matches the primer, or
# flushed if it doesn't.
my @record;

# this will continue one or more output handles to which the matching
# record is printed
my @handles;

# start iterating over file
open my $fh, '<', $infile or die $!;

# keep track of how many hits 
my %seen;

while(<$fh>) {
	chomp;

	# this is the ID line
	if ( $line % 4 == 0 ) {
		if ( /^\@/ ) {
			@record  = ( $_ );
			@handles = ();
		}
		else {
			die "ID line should begin with \@, not $_";
		}
	}

	# this is the sequence line
	elsif ( $line % 4 == 1 ) {
		for my $name ( keys %primers ) {
			my $distance = adist( $primers{$name}, $_ );
			if ( $distance <= $mismatch ) {
				push @handles, $handles{$name};
				$seen{$name}++;
			}
		}
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
		for my $handle ( @handles ) {
			print $handle join("\n", @record), "\n" if not $test;
		}
	}
	$line++;
}
use Data::Dumper;
print Dumper(\%seen);
