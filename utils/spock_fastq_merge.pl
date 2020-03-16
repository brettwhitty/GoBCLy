#!/usr/bin/env perl

use strict;
use warnings;

$| = 1;

use Cwd qw{ abs_path };
use File::Temp qw{ tmpnam };
use Getopt::Long;
use PerlIO::gzip;
use POSIX qw{ mkfifo };

my $fastq_barcodes_file = '';
my $fastq_reads_file = '';
my $named_pipe = 0;
my $tmp_file = 0;
my $gzip = 0;
my $output_file = '';
my $force = 0;

GetOptions(
    'barcodes|b=s'  =>  \$fastq_barcodes_file,
    'reads|r=s'     =>  \$fastq_reads_file,
    'pipe|p!'       =>  \$named_pipe,
    'tmp|t!'        =>  \$tmp_file,
    'gzip|z!'       =>  \$gzip,
    'output|o=s'    =>  \$output_file,
    'force|f!'      =>  \$force,
);

if ($output_file ne '' && -f $output_file) {
    die "Specified output file '$output_file' already exists, use --force/-f flag to overwrite!" unless $force;
}

unless ($fastq_barcodes_file ne '' && -s $fastq_barcodes_file && $fastq_reads_file ne '' && -s $fastq_reads_file) {
    print STDERR $0." --barcodes/-b </PATH/TO/BARCODES_R1.FASTQ.GZ> --reads/-r </PATH/TO/READS_R1.FASTQ.GZ> [--pipe/-p]\n";
    exit 1;
}

## expand these in case relative or symlinks
$fastq_barcodes_file = abs_path($fastq_barcodes_file);
$fastq_reads_file   = abs_path($fastq_reads_file);


my $output_file_name = $output_file || '/dev/stdout';
my $write_mode = '>';
if ($tmp_file) {
    $output_file_name = tmpnam();
    print $output_file_name;
}
if ($named_pipe) {
    $output_file_name = tmpnam();
    mkfifo($output_file_name, 0700)
        || die "mkfifo $output_file_name failed: $!";
    $write_mode = '>>';
    print $named_pipe;
}
if ($gzip) {
    $write_mode .= ':gzip';
}

do_merge($fastq_barcodes_file, $fastq_reads_file, $output_file_name, $write_mode);

sub do_merge {
    my ($fastq_barcodes_file, $fastq_reads_file, $output_file_name, $write_mode) = @_;

    open my $outfh, $write_mode, $output_file_name
        or die $!;

    open my $infh_barcodes, '<:gzip', $fastq_barcodes_file
        or die $!;

    open my $infh_reads, '<:gzip', $fastq_reads_file
        or die $!;

    while (1) {
        my $barcode = get_record($infh_barcodes);
        my $read = get_record($infh_reads);
        if (!defined($barcode) && !defined($read)) {
            last;
        }
        print $outfh join("\n", (
                $read->[0] . ' ' . $barcode->[1],
                $barcode->[1] . $read->[1],
                $read->[2],
                $barcode->[3] . $read->[3],
        ))."\n";
    }
    close $outfh;
}

sub get_record {
    my ($fh) = @_;

    my @lines = ();
    for (my $i = 0; $i < 4; $i++) {
        my $line = <$fh>;
        chomp $line;
        #if (eof()) { return undef; }
        push(@lines, $line);
    }

    return \@lines;
}
