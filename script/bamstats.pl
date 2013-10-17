#!/usr/bin/perl
use strict;
use warnings;
use List::Util 'sum';
use Getopt::Long;
use Bio::DB::Sam;

# global variable: verbosity
my $Verbosity = 1;

# basic logging functionality
sub LOG ($$) {
    my ($msg,$method) = @_;
    my ( $package, $file1up, $line1up, $sub ) = caller( 2 );
    my ( $pack0up, $file, $line, $sub0up )    = caller( 1 );
    my $log = sprintf "%s %s [%s %s] - %s\n", uc $method, $sub || '', 'line', $line, $msg;
    print STDERR $log;
}
sub DEBUG ($) { LOG shift, 'DEBUG' if $Verbosity >= 3 }
sub INFO ($)  { LOG shift, 'INFO'  if $Verbosity >= 2 }
sub WARN ($)  { LOG shift, 'WARN'  if $Verbosity >= 1 }
sub ERROR ($) { LOG shift, 'ERROR'; exit 1 }

# process command line arguments
sub process_args {
    my ( $bam, $ref, $coordinate );
    GetOptions(
        'bam=s'        => \$bam,
        'verbose+'     => \$Verbosity,
        'ref=s'        => \$ref,
        'coordinate=s' => \$coordinate,
    );
    
    # communicate settings
    INFO "BAM file: $bam";
    INFO "REF fasta: $ref";
    INFO "Interval: $coordinate";

    # instantiate helper objects
    my $db = Bio::DB::Sam->new(
        '-bam'   => $bam,
        '-fasta' => $ref,
    );
    
    # parse coordinate
    my ( $chromo, $start, $end );
    if ( $coordinate =~ /^([^:]+):(\d+)-(\d+)$/ ) {
        ( $chromo, $start, $end ) = ( $1, $2, $3 );
    }
    else {
        ERROR "Can't parse coordinate $coordinate";
    }
    
    # return result object
    return 
        'ref'    => $ref, 
        'db'     => $db, 
        'chromo' => $chromo,
        'start'  => $start,
        'end'    => $end;
}

# iterates over reads inside interval, computes 
# coverage and mapping quality
sub compute_coverage_and_quality {
    my %args = @_;
    
    # make interval for get_features_by_location
    my %loc = (
        '-seq_id' => $args{'chromo'},
        '-start'  => $args{'start'},
        '-end'    => $args{'end'},      
    );  

    # we compute coverage by counting the number of aligned reads with $cover
    # and computing their average length by summing the lengths
    my $reads  = 0;
    my $length = 0;

    # we record all distinct starting positions of aligned reads with %starts. 
    # subsequently we will calculate the base composition one base upstream from 
    # that on the reference. if there is fragmentation this should show an AT 
    # bias on the reference.
    my %starts;

    # we record all qualities in @quals, then average that
    my @quals;

    # iterate over alignments within interval
    for my $aln ( $args{'db'}->get_features_by_location(%loc) ) {

        # start and end coordinates of the alignment on the reference
        my ( $start, $end ) = ( $aln->start, $aln->end );
        $starts{ $start } = 1;
    
        # running tally of qualities
        push @quals, $aln->qual;
    
        # running tally of alignment lengths
        $length += ( $end - $start );
    
        # number of reads
        $reads++;
        INFO "Read: $reads" unless $reads % 100000;
    }

    # compute coverage
    my $cover = ( $reads * ( $length / $reads ) ) / ( $loc{'-end'} - $loc{'-start'} );

    # compute average quality
    my $qual = sum(@quals) / scalar(@quals);
    
    # communicate results
    INFO "Total reads: $reads";
    INFO "Average cover: $cover";
    INFO "Average mapping quality: $qual";
    INFO "Going to lookup ".scalar(keys(%starts))." upstream bases";
    return $qual, $cover, keys %starts;
}

# given a fasta file name, chromosome name and a list of locations, computes
# compositional bias one base upstream from the locations
sub compute_upstream_bias {
    my %args = @_;
    INFO "computing compositional bias in seq $args{chromo} from file $args{ref}";
    
    # open file handle
    open my $fh, '<', $args{'ref'} or ERROR $!;
    
    # concatenate record
    my $seq;
    while(<$fh>) {
        chomp;
        if ( /^>$args{chromo}\s*$/ ) {
            $seq = '';
            next;
        }
        $seq .= $_ if defined $seq;
    }
    
    # compute content. locations are 1-based so decrement to go 1 upstream
    my %content;
    for my $i ( @{ $args{'locations'} } ) {
        my $base = substr $seq, $i - 1, 1;
        $content{uc $base}++;
    }
    
    # calculate fraction of A's and T's on the total
    my $at = sum(grep { defined $_ } @content{qw[A T]})/sum(values(%content));
    INFO "AT bias on reference, one base upstream: $at";
    return $at;
}

sub main {
    my %args = process_args();
    my ( $qual, $cover, @locations ) = compute_coverage_and_quality( %args );
    my $bias = compute_upstream_bias( %args, 'locations' => \@locations );
    print <<"RESULT";
Average cover: $cover
Average mapping quality: $qual
AT bias on reference: $bias
RESULT
}
main();

