#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util 'sum';

use Bio::SeqIO;
use Bio::DB::Sam;
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $window = 1000;
my $verbosity = WARN;
my $bam;   # input file, BAM/SAM format
my $fasta; # input file, FASTA format
my $cover; # flag: write coverage in window
my $phred; # flag: write phred score in window
my $map;   # flag: write mapping quality in window
my $snp;   # flag: write SNPs in window
GetOptions(
	'window=i' => \$window,
	'bam=s'    => \$bam,
	'verbose+' => \$verbosity,
	'fasta=s'  => \$fasta, 
	'cover'    => \$cover,
	'phred'    => \$phred,
	'map'      => \$map,
	'snp'      => \$snp,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level'   => $verbosity,
	'-class'   => 'main',
);
my $sam = Bio::DB::Sam->new(
	'-bam'     => $bam,
	'-fasta'   => $fasta,
);
my $seqio = Bio::SeqIO->new(
	'-format'  => 'fasta',
	'-file'    => $fasta,
);

# iterate over sequences in the reference genome
while ( my $seq = $seqio->next_seq ) {
	my $id     = $seq->id;
	my $length = $seq->length; 

	# iterate over windows
	for ( my $start = 0; $start <= $length; $start += $window ) {

		# don't want to run beyond the length of the chromosome
		# for the end coordinate
		my $end = $start + $window < $length ? $start + $window : $length;

		# get all mapped reads in the window
		my @alignments = $sam->get_features_by_location(
			'-seq_id' => $id,
			'-start'  => $start,
			'-end'    => $end,
		);

		# create record of results
		my @record = ( $id, $start, $end );

		# add the coverage within this window
		if ( $cover ) {
			push @record, scalar @alignments;
		}

		# add the average phred score within this window
		if ( $phred ) {
			my @phreds = map { $_->qscore } @alignments;
			push @record, sum(@phreds) / scalar(@phreds);
		}

		# add the average mapping quality within this window
		if ( $map ) {
			my @quals = map { $_->qual } @alignments;
			push @record, sum(@quals) / scalar(@quals);
		}

		# add the SNP count within this window
		if ( $snp ) {
		
			# create callback to pass into pileup
			my %matches;
			my $fetch_back = sub {
				my ( $seqid, $pos, $p ) = @_;

				# get the reference base
				my $r = $sam->segment($seqid,$pos,$pos)->dna;

				# iterate over pileups
				for my $pileup ( @{ $p } ) {
					my $a    = $pileup->alignment;
					my $qpos = $pileup->qpos;
					my $dna  = $a->query->dna;
					my $base = $pileup->indel == 0 ? substr($dna,$qpos,1)
						: $pileup->indel >  0 ? '*'
						: '-';
					$matches{matched}++ if $r eq $base;
					$matches{total}++;
				}
    			};
			$sam->pileup( "${id}:${start}-${end}", $fetch_back );
			push @record, $matches{matched} / $matches{total};
		}
		
		# done, write output
		print join( "\t", @record ), "\n";

	}
}
