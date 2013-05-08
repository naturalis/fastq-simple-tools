#!/usr/bin/perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Path 'make_path';

=pod

=head1 TITLE

genome_chopper.pl - splits a scaffolded genome into contigs for GenBank submission

=head1 SYNOPSIS

 genome_chopper.pl -split -fasta <infile>

=head1 ARGUMENTS AND OPTIONS

=over

=item B<-fasta "infile">

Input FASTA file containing scaffolds.

=item B<-workdir "workdir">

Directory to which to write the contigs and AGP coordinates. Optional, default is 
a directory called C<submission> below the current working directory.

=item B<-info "info file">

An input file with extra metadata (key/value pairs, one on each line) that are written 
to the FASTA headers of the contigs. Example:

 organism=Ophiophagus hanna

This will result in FASTA headers that contain C<[organism=Ophiophagus hanna]>.
Optional. Refer to the NCBI documentation for more info about possible key/value pairs:
L<http://www.ncbi.nlm.nih.gov/Sequin/modifiers.html>

=item B<-prefix "four letter code">

All generated placeholder identifiers for contigs will start with this code. Optional,
default is C<OPHA>.

=item B<-split>

Writes the contigs as individual FASTA files in the workdir. Optional. The default is
to write the FASTA data to one big file. Both appear to be legal for GenBank submissions,
but the big file might be difficult for C<tbl2asn> (see below) to handle.

=item B<-outfile "outfile.fsa">

Writes the contigs as a single, concatenated FASTA file in the workdir. Optional,
default is C<submission.fsa>.

=item B<-agp "outfile.agp">

Writes the coordinates of contigs and gaps to an AGP2.0 file in the workdir. Optional,
default is C<submission.agp>.

=item B<-minlength "199">

Specifies the length that contigs must I<exceed> in order to be included in the output.
Optional, default is 199.

=item B<-verbose>

Increases the loging verbosity. Use once for informational messages and twice for info
as well as debug messages (which is very verbose!). Optional.

=item B<-help>

Prints this help message and quits.

=back

=head1 DESCRIPTION

GenBank submissions of whole-genome shotgun sequencing ("WGS") projects need to consist
of FASTA records of all individual contigs (with no NNNs) and tabular data that places
these contigs and the gaps between them on their containing scaffolds. This script takes
a big FASTA file (presumably a genome consisting of scaffolds) and chops it up into 
submission-ready contigs (FASTA files with the right header formatting) and a coordinates
file in "AGP2.0" format. After running this script you should be able to then run NCBI's
C<tbl2asn> tool to create the submission files, e.g. by invoking it as:

 tbl2asn -p submission -t template.sbt -a s -V v -Z discrep

Where C<-p submission> refers to the workdir to which your FASTA data was written, and
C<-t template.sbt> refers to an ASN.1 file with submission metadata that can be generated
either using Sequin or using the web form at 
L<http://www.ncbi.nlm.nih.gov/WebSub/template.cgi>

=head1 SEE ALSO

=over

=item The C<tbl2asn> documentation: L<http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/>

=item The WGS submission docs: L<http://www.ncbi.nlm.nih.gov/genbank/wgs.submit>

=item Eukaryotic Genome Submission Guide: 
L<http://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission>

=back

=cut

# process command line arguments
my $fasta;                           # input fasta file with concatenated scaffolds
my $agp       = 'submission.agp';    # AGP2.0 file
my $verbosity = 1;                   # logging level, default is WARN
my $minlength = 199;                 # minimum contig length
my $workdir   = 'submission';        # directory to write to
my $info      = 'info.ini';          # extra metadata file with key/value pairs
my $prefix    = 'OPHA';              # four-letter prefix
my $outfile   = 'submission.fsa';    # output file of contigs
my $split;                           # split into separate fasta files
GetOptions(
	'fasta=s'     => \$fasta,
	'agp=s'       => \$agp,
	'verbose+'    => \$verbosity,
	'minlength=i' => \$minlength,
	'workdir=s'   => \$workdir,
	'info=s'      => \$info,
	'outfile=s'   => \$outfile,
	'split'       => \$split,
	'help|?' => sub { pod2usage( '-verbose' => 2 ) },
);

# basic logging functionality
sub LOG ($$) {
	my ($msg,$method) = @_;
	my ( $package, $file1up, $line1up, $sub ) = caller( 2 );
	my ( $pack0up, $file, $line, $sub0up )    = caller( 1 );
	my $log = sprintf( "%s %s [%s, %s] - %s\n", uc $method, $sub || '', $file, $line, $msg );
	print STDERR $log;
}
sub DEBUG ($) { LOG shift, 'DEBUG' if $verbosity >= 3 }
sub INFO ($)  { LOG shift, 'INFO'  if $verbosity >= 2 }
sub WARN ($)  { LOG shift, 'WARN'  if $verbosity >= 1 }

# read metadata
my %info   = read_info( $info );
my $object = 1; # global counter for generated contig IDs

# create workdir if not exist
if ( not -d $workdir ) {
	INFO "going to create working directory $workdir";	
	make_path( $workdir );
}

# create output handles
open my $fastaFH, '>', "${workdir}/${outfile}" or die $!;
my $agpFH = create_agp( %info );

# initializes AGP file with timestamped header, returns file handle
sub create_agp {
	my %info = @_;
	my $date = create_timestamp();
	open my $agpFH, '>', "${workdir}/${agp}" or die $!;
	print $agpFH <<"HEADER";
##agp-version	2.0
# ORGANISM: $info{organism}
# ASSEMBLY DATE: ${date}
# DESCRIPTION: AGP specifying the assembly of scaffolds from WGS contigs
HEADER

return $agpFH;
}

# creates the current time stamp in the format
# suggested by the AGP2.0 format spec, i.e.
# dd-Month-year
sub create_timestamp {
	my @months = qw(
		January
		February
		March
		April
		May
		June
		July
		August
		September
		October
		November
		December
	);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $day = sprintf '%02d', $mday;
	$year += 1900;
	my $month = $months[$mon];
	return "${day}-${month}-${year}";
}

# read the data
read_fasta( $fasta );

# read file with key/value pairs, returns hash 
sub read_info {
	my $info = shift;
	
	# an info file was provided
	if ( -e $info ) {
		INFO "going to read metadata from $info";
		my %result;
		
		# open a file handle
		open my $fh, '<', $info or die $!;
		while(<$fh>) {
			chomp;
			
			# assume simple INI-style key/value pairs
			my ( $key, $value ) = split /=/, $_;
			$result{$key} = $value;
		}
		return %result;
	}
	else {
		INFO "no additional metadata provided";
	}
}

# read large FASTA file, dispatch each encountered record to the scaffold handler
sub read_fasta {
	my $file = shift;
	INFO "going to read FASTA data from $fasta";
	
	# these will iteratively be populated with record names and data
	my ( $name, $seq );
	
	# open a file handle
	open my $fh, '<', $file or die $!;
	
	# start reading from the handle
	while(<$fh>) {
		chomp;
		
		# line contains a FASTA defline
		if ( />(\S+)/ ) {
			my $newname = $1;
			
			# handle the last seen record, if there is one
			handle_scaffold( $name => $seq ) if $name and $seq;
			
			# start a new record
			$name = $newname;
			$seq = '';
		}
		else {
		
			# grow the current record
			$seq .= $_;
		}
	}
	
	# handle the last record, this because no next defline
	# was seen so the handler wasn't triggered for the last record
	handle_scaffold( $name => $seq ) if $name and $seq;
}

# write FASTA data. with the -split flag this will write to a separate
# file each time, otherwise to one big file. Includes the key/value pairs
# from the info file in the FASTA header
sub write_fasta {
	my ( $name, $seq ) = @_;
	
	# create or re-use handle depending on whether the
	# split argument was used
	my $outFH;
	if ( $split ) {
		open $outFH, '>', "${workdir}/${name}.fsa" or die $!;
	}
	else {
		$outFH = $fastaFH;
	}
	
	# start defline
	print $outFH ">$name";
	
	# add extra info
	if ( %info ) {
		print $outFH ' ', join ' ', map { "[$_=$info{$_}]" } keys %info;
	}
	
	# line break
	print $outFH "\n";
	
	# split the seq into 80-char lines
	my @lines = unpack("(A80)*", $seq);
	
	# write the lines
	print $outFH join("\n", @lines), "\n";

	DEBUG "wrote FASTA record $name";
}

# splits a scaffold into NNN-separated contigs. dispatches the contig seq data
# to the FASTA writer and the coordinates to the AGP writer.
sub handle_scaffold {
	my ( $name, $seq ) = @_;
	my $original_length = length($seq);
	INFO "going to handle scaffold $name with $original_length base pairs";
	
	# check minimum length
	if ( length($seq) > $minlength ) {

		# we collect all the contigs and their coordinates into this array,
		# all stretches between them will be written as known missing
		my @contigs;
		
		# split the seq on bases, this so we can reconstruct the 
		# original coordinates
		my @bases = split //, $seq;
		my $start = 0;
		
		# in principle we iterate over all bases in the outer loop
		CONTIG: while ( $start < @bases ) {
		
			# we have arrived at a contig
			if ( $bases[$start] ne 'N' ) {
			
				# seek for the end of the contig
				my $stop = $start + 1;
				while( $stop < @bases and $bases[$stop] ne 'N' ) {
					$stop++;
				}
				
				# $stop has counted to the first N inclusive, so 
				# go one step back
				my $contigseq  = join '', @bases[ $start .. $stop - 1];
				
				if ( length($contigseq) > $minlength ) {
					my $contigname = $prefix . sprintf('%08d', $object++) . '.1';
					write_fasta( $contigname => $contigseq );
					
					# store as pseudo-object
					push @contigs, { 
						'name'   => $contigname, 
						'start'  => $start + 1, 
						'stop'   => $stop,
						'length' => length($contigseq), 
					};
				}
				else {
					DEBUG "contig is too short, skipping";
				}
				
				$start = $stop;
				next CONTIG;
			}
			$start++;
		}
		
		# now write the AGP coordinates
		write_agp($name, @contigs);
	}
	else {
		DEBUG "scaffold $name is too short: $original_length <= $minlength";
	}	
}

sub write_agp {
	my ( $name, @contigs ) = @_;
	my $component = 1;
	for my $i ( 0 .. $#contigs ) {
	
		# print the contig itself
		print $agpFH 
			$name, "\t",                    # the scaffold name
			$contigs[$i]->{'start'}, "\t",  # the contig's start within the scaffold
			$contigs[$i]->{'stop'}, "\t",   # the contig's end within the scaffold
			$component++, "\t",             # counts contigs and gaps
			'W', "\t",                      # component type: WGS contig
			$contigs[$i]->{'name'}, "\t",   # locally unique contig ID
			1, "\t",                        # the contig's local start
			$contigs[$i]->{'length'}, "\t", # the contig's local end
			'?', "\n";                      # orientation unknown
		
		# print subsequent gap, if any
		if ( $i < $#contigs and $contigs[$i]->{'stop'} < $contigs[$i+1]->{'start'} ) {
			my $start = $contigs[$i]->{'stop'} + 1;
			my $stop  = $contigs[$i+1]->{'start'}  - 1;
			print $agpFH 
				$name, "\t",                    # the scaffold name
				$start, "\t",                   # the gap's start within the scaffold
				$stop, "\t",                    # the gap's end within the scaffold
				$component++, "\t",             # counts contigs and gaps
				'N', "\t",                      # component type: gap with specified size
				( $stop - $start ) + 1, "\t",   # gap length
				'scaffold', "\t",               # a gap between two contigs in a scaffold
				'yes', "\t",                    # evidence of linkage between lines
				'paired-ends', "\n";            # linkage due to paired ends		
		}
	}	
}

