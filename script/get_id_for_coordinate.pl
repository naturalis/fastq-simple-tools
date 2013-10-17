#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my ( $gff3, $coordinate );
GetOptions(
	'verbose+'     => \$verbosity,
	'gff3=s'       => \$gff3,
	'coordinate=s' => \$coordinate,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# parse coordinate
my ( $chromo, $start, $end );
if ( $coordinate =~ /^([^:]+):(\d+)-(\d+)$/ ) {
	( $chromo, $start, $end ) = ( $1, $2, $3 );
	$log->info("will look for chromosome $chromo, start=$start, end=$end");
}
else {
	die "Can't parse coordinate $coordinate\n";
}

# open file handle
$log->info("going to scan GFF3 file $gff3 for IDs that span coordinate $coordinate");
open my $fh, '<', $gff3 or die $!;
while(<$fh>) {
	chomp;
	my @record = split;
	
	# column numbers in gff3
	my $chr_idx   = 0;
	my $start_idx = 3;
	my $end_idx   = 4;
	my $desc_idx  = 8;
	my $type_idx  = 2;	
	
	# we're on the right chromosome
	if ( $record[$chr_idx] eq $chromo ) {
		
		# we're in the right interval
		if ( $record[$start_idx] < $start and $record[$end_idx] > $end ) {
			$log->debug("found annotation that spans interval $start..$end");
			
			# we've found a valid annotation
			if ( $record[$desc_idx] and $record[$desc_idx] =~ /ID=([^;]+);/ ) {
				my $id = $1;
				print $coordinate, "\t", $id, "\t", $record[$type_idx], "\n";
			}
		}
	}
}