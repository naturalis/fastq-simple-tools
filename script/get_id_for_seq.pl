#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Tools::Run::StandAloneBlastPlus;

# globals
my $queryfile; # FASTA file with loci to search for
my $dbfile;    # local BLAST database file
my $gff3;      # GFF3 genome annotation file
my $log;       # logger
my $range = 5; # range by which the intersection may be offset
my $feature = 'gene'; # feature type of interest
my $Verbosity = 1; # logging level

# basic logging functionality
sub LOG ($$) {
	my ($msg,$method) = @_;
	my ( $package, $file1up, $line1up, $sub ) = caller( 2 );
	my ( $pack0up, $file, $line, $sub0up )    = caller( 1 );
	my $log = sprintf( "%s %s [line %s] - %s\n", uc $method, $sub || '', $line, $msg );
	print STDERR $log;
}
sub DEBUG ($) { LOG shift, 'DEBUG' if $Verbosity >= 3 }
sub INFO ($)  { LOG shift, 'INFO'  if $Verbosity >= 2 }
sub WARN ($)  { LOG shift, 'WARN'  if $Verbosity >= 1 }

# process command line arguments
sub process_args {

	my $help; # flag to display help message	
	GetOptions(
		'queryfile=s' => \$queryfile,
		'dbfile=s'    => \$dbfile,
		'verbose+'    => \$Verbosity,
		'help'        => \$help,
		'gff3=s'      => \$gff3,
		'feature=s'   => \$feature,
		'range=i'     => \$range,
	);

	# too complicated to use, need help message:
	if ( $help ) {
		print <<"HELP";
	Usage: $0 -query <fasta file> -db <local blast file> -gff3 <gff3 file> \
		[-type <gene|exon|CDS|mRNA|intro|five_prime_UTR|three_prime_UTR>]

	The query file contains one or more fasta sequences for interesting loci.
	These are then blast searched against the reference genome so that we
	can reconstruct the coordinates of the loci. We then search for these 
	coordinates in the GFF3 file and report back the nearest identifier(s).
HELP
		exit 0;
	}
}

# iterate over sequences in input file, dispatch each to blast
sub read_query_sequences {
	my $queryfile = shift;
	
	# open file reader
	DEBUG "instantiating file reader for $queryfile in FASTA format";
	my $io = Bio::SeqIO->new(
		'-file'   => $queryfile,
		'-format' => 'fasta',
	);	

	# assemble results here
	my @result;

	# iterate over sequences in query file. the query file contains loci
	# whose coordinates into the reference genome we are going to fetch
	while( my $seq = $io->next_seq ) {
		my $id = $seq->id;
		INFO "going to BLAST sequence '$id'";
		my @report = run_blast($id,$seq);
		push @result, \@report;
	}
	return @result;
}

# run blast, dispatch each hit to GFF3 scanner
sub run_blast {
	my ($id,$seq) = @_;
	
	# open BLAST wrapper
	DEBUG "instantiating standalone BLAST+ with database $dbfile";
	my $blast = Bio::Tools::Run::StandAloneBlastPlus->new(
		'-db_data' => $dbfile,
	);
	
	# run BLAST
	my $result = $blast->blastn( '-query' => $seq );
	my $resultname = $result->query_name; # should be same as $id
	INFO "found result for '$resultname'";

	# this will be serialized to the final report
	my @report;

	# iterate over hits
	while( my $hit = $result->next_hit ) {
		my $hitname = $hit->name;
		$hitname =~ s/^lcl\|//; # strip prefix
		INFO "found hit '$hitname'"; # should be chromosome

		# iterate over high scoring pairs
		while( my $hsp = $hit->next_hsp ) {
			my $start = $hsp->start("subject");
			my $end   = $hsp->end("subject");
			INFO "HSP range is $start..$end (".($end-$start)." bp)";

			# having found the coordinates of the focal
			# hit into the reference we will now
			# query the GFF3 file for the same coordinates
			my @result = scan_gff3 (
				'start' => $start,
				'end'   => $end,
				'chr'   => $hitname, 
				'id'    => $id,
			);
	
			# add result to report
			push @report, {
				'region'  => "${hitname}:${start}-${end}",
				'score'   => $hit->raw_score,
				'matches' => \@result,
				'query'   => $id,
			} if @result;
		}
	}
	return @report;
}

# scan GFF3 file, return all intersecting features of the right type
sub scan_gff3 {
	my %args = @_;
	INFO "going to search for IDs on $args{chr} between $args{start} and $args{end}";
	
	# column numbers in gff3
	my $chr   = 0;
	my $type  = 2;	
	my $start = 3;
	my $end   = 4;
	my $desc  = 8;
	
	# collect IDs here
	my @result;
	
	# open handle to GFF3
	open my $fh, '<', $gff3 or die $!;
	while(<$fh>) {
		my @record = split;
		
		# we are on the right chromosome
		if ( $record[$chr] eq $args{'chr'} ) {
			
			# we are within range
			if ( ($record[$start]-$range) < $args{'start'} && ($record[$end]+$range) > $args{'end'} ) {
				DEBUG "found region $args{'start'}..$args{'end'}";
				
				# the feature is of the requested type
				if ( $record[$type] eq $feature ) {		

					# parse the ID from the description
					if ( $record[$desc] =~ m/ID=([^;]+)/ ) {
						my $feature_id = $1;
						push @result, $feature_id;
						INFO "found $feature_id $record[$chr]:$record[$start]-$record[$end] ";						
					}
				}
				else {
					DEBUG "focal feature is wrong type: $record[$type] ne $feature";
				}
			}
		}
	}
	return @result;
}

# initialize program, dispatch, report results
sub main {
	process_args();
	my @result = read_query_sequences($queryfile);
	print Dumper(\@result);
}
main();
