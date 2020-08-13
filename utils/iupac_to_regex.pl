#!/usr/bin/env perl

$| = 1;

use strict;
use warnings;

## Quick script to convert IUPAC ambiguity codes into PCRE syntax.
##
## For use with dbSNP-masked hg19, Motto, etc.
##
## Brett Whitty <brettwhitty@gmail.com>, all rights retained

use Getopt::Long;

my $fasta = 1;
my $to_uc = 0;
my $five_rev = 0;
my $three_rev = 1;
my $do_subseq = 0;
my $subseq_len = 60;
my $subseq_step = 20;
my $subseq_minlen = 40;
my $hyperscan_format = 0;
my $zero_pad = 0;
my $re_edit_dist = 0; # 2; ## TODO: fix flag support in the hikeeba-personafi dev code
my $re_hamm_dist = 0;
my $re_min_match = 0;
my $re_or = 0;
my $re_flags = '';
my $re_ext_flags = '';
my $int_idx = 1;

GetOptions(
    'fasta|f!'      =>  \$fasta,            ## fasta output
    'uppercase|u!'  =>  \$to_uc,            ## to uppercase
    'five|5=i'      =>  \$five_rev,         ## Assuming input is 5' => 3' (+ strand), controls reversing seq
                                            ##     (0 = forward, 1 = reverse, [^01] = both)
    'three|3=i'     =>  \$three_rev,        ## minus-strand complement and reversed
                                            ##     (0 = complement, 1 = reverse, [^01] = both)
    'subseq|s!'     =>  \$do_subseq,        ## create regexes for subsequences
    'len|l=i'       =>  \$subseq_len,       ## subsequence length
    'step|t=i'      =>  \$subseq_step,      ## step value for subsequences
    'hs!'           =>  \$hyperscan_format, ## generate a hyperscan-compatible regex list
                                            ## instead of FASTA-style output
    'zero|z!'       =>  \$zero_pad,         ## zero pad re idx instead of adding 10 to chr id number
    'int|I!'        =>  \$int_idx,          ## integer index of output; don't parse input seq IDs for tokens 
    'edit=i'        =>  \$re_edit_dist,     ## edit distance
    'hamming=i'     =>  \$re_hamm_dist,     ## hamming distance
    'match=i'       =>  \$re_min_match,     ## min match
    'or!'           =>  \$re_or,            ## concatenate patterns from subseqs with |
    'flags=s'       =>  \$re_flags,         ## additional flags for regexes (eg: 'i', 'o')
    'ext-flags=s'   =>  \$re_ext_flags,     ## extended flags
##    'bool!'         =>  \$bool,	    ## boolean or
);

## if user has provided flags, check they're valid
if (length($re_flags) > 0) { 
	$re_flags = check_flags($re_flags); 
}

print STDERR "## WARNING!!! This script expects header & sequence to each be on a single line\n";

my $input_header = '';
my @seq = ();

my $re_idx_counter = {};
my $re_idx = '';

while (<>) {

    chomp;

    if (/^>/) {

        $input_header = $_;

	if ($int_idx) {
	        $re_idx = ++$re_idx_counter->{'_'};
	} else {
	        $input_header =~ /^>([^:]+)/;
	        my $chr_id = $1;

	        ## initialize regular expression index for doing hyperscan compatible output
	        if (! defined($re_idx_counter->{$chr_id})) {
	            $re_idx_counter->{$chr_id} = chr_to_idx($chr_id);
        	}

	        ## id used when writing regex with --hs flag
	        $re_idx = ++$re_idx_counter->{$chr_id};
	}

        ## append regex family id
        $input_header .= " idx_re=$re_idx";

    } else {

        my $five_fwd  = $_;
        my $five_rvs = reverse($five_fwd);
        my $three_fwd = complement_dna($five_fwd);
        my $three_rvs = reverse($three_fwd);

        #if (0) {

        my @orients = ();

        if ($five_rev == 0 || $five_rev >= 2)  {
            push(@orients, "fwd\t$five_fwd");
        }
        if ($five_rev == 1 || $five_rev >= 2) {
            push(@orients, "rvs\t$five_rvs");
        }

        if ($three_rev == 0 || $three_rev >= 2)  {
            push(@orients, "comp_fwd\t$three_fwd");
        }
        if ($three_rev == 1 || $three_rev >= 2) {
            push(@orients, "comp_rvs\t$three_rvs");
        }

        #my @parens = ();
        #foreach my $seq(@orients) {
        #    push(@parens, do_iupac(parenthesize($seq)));
        #}

        #print '('.join('|', @parens).')'."\n";
        #}

        foreach my $i (@orients) {
            my ($orient, $seq) = split(/\t/, $i, 2);

            my $orient_header = $input_header;
            $orient_header =~ s/^(>\S+)/$1.$orient/;

            my @outseqs = ();
            if ($do_subseq) {
                my $workseq = $seq;
                while (length($workseq) >= $subseq_len) {
                    push(@outseqs, substr($workseq, 0, $subseq_len));
                    $workseq = substr($workseq, $subseq_len);
                }
                if (length($workseq) > 0 && length($workseq) <= $subseq_minlen) {
                    push(@outseqs, $workseq);
                }
            } else {
                push(@outseqs, $seq);
            }
            my @subseq_patterns;
            my $subseq_idx = 0;
            foreach my $outseq(@outseqs) {
                if ($do_subseq && scalar(@outseqs) > 1) {
                    my $subseq_header = $orient_header;
                    ## more human-friendly
                    my $subseq_id = $subseq_idx + 1;
                    $subseq_header =~ s/^(>\S+)/$1.$subseq_id/;
                    print $subseq_header." offset=".($subseq_idx * $subseq_len)." subseq_len=".length($outseq) unless $hyperscan_format;
                    $subseq_idx++;
                } else {
                    print $orient_header unless $hyperscan_format;
                }
                print "\n" unless $hyperscan_format;
                my $pattern = do_iupac($outseq);

                ## TODO: add sensible flags here
                #
                # i = case insensitive
                # s = . matches newlines
                # m = multiline anchoring
                # o* = report match ID at most once
                # e* = allow patterns that can match empty buffers
                # u* = UTF-8 mode
                # p* = unicode property support
                # f* = prefiltering mode
                # l* = leftmost start of match reporting
                # C = logical combination of patterns
                # Q = quiet at matching
                #
                ## * this is 'H' in the hyperscan manual
                ## * this is 'V' in the hyperscan manual
                ## ...
                ##
                ##
                #
                # ... plus extended parameters:
                #
                # /pattern/flags{key=value,...}
                #
                #   flags = flags governing which of the other fields in the structure are used (?!?)
                #   min_offset = the minimum end offset in the data stream at which the expression should match
                #   max_offset = the maximum end offset in the data stream at which the expression should match
                #   min_length = minimum length from start to end required to successfully match the expression
                #   edit_distance = Levenshtein distance
                #   hamming_distance = Hamming distance
                #
                #my $re_flags = 'iH';
                ### 'i' shouldn't be necessary for FASTQ
                #$re_flags = 'o';

                ## TODO: cleanup
                ## setting these locally from flags for no particular reason
                my $levenshtein_distance = $re_edit_dist;
                my $hamming_distance = $re_hamm_dist;
                my $min_match = $re_min_match;

                ## TODO: parser will choke on these at the moment if used
                $re_ext_flags = join('', (
                    '{',
                    join(','), (
                       ($levenshtein_distance) ?
                         "edit_distance=$levenshtein_distance"
                         : '',
                       ($hamming_distance) ?
                         "hamming_distance=$hamming_distance"
                         : '',
                       ($min_match) ?
                         "min_length=$min_match"
                         : '',
                     ),
                    '}',
                ));

                ## if empty, don't bother
                $re_ext_flags =~ s/^\{\}$//;

                if ($hyperscan_format) {
                    print join('', (
                            $re_idx,
                            ':',
                            '/',
                            $pattern,
                            '/',
                            $re_flags,
                            $re_ext_flags,
                    ))."\n" unless $re_or;

                    push(@subseq_patterns, $pattern) if $re_or;

                } else {
                    print $pattern."\n"
                }
            }
            if ($hyperscan_format && $re_or) {
                ##
                ## reordering these here from left => right to "spiral" from center
                ##
                ## TODO: these should be optimized for search efficiency,
                ## eg: allow user to specify nested ORs from left => right, right => left;
                ## shrinking/expandig subseq window from outside in / inside out, etc.
                my @temp = ();
                while (@subseq_patterns) {
                    if (@subseq_patterns) {
                        push(@temp, shift @subseq_patterns);
                    }
                    if (@subseq_patterns) {
                        push(@temp, pop @subseq_patterns);
                    }
                }
                @subseq_patterns = reverse @temp;
                print join('', (
                    $re_idx,
                    ':',
                    '/',
                    '('.join('|', @subseq_patterns).')',
                    '/',
                    $re_flags,
                    $re_ext_flags,
                ))."\n";
            }
        }
    }
}

## returns complement of a DNA sequence, including IUPAC ambiguity codes;
## preserves upper/lower case
sub complement_dna {
    my ($r_seq) = @_;

    $r_seq =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb/TtGgCcAaKkYyWwSsRrMmBbDdHhVv/;

    return $r_seq;
}

## converts strings with IUPAC codes to regex strings with character sets, etc.
sub do_iupac {

    my ($seq) = @_;

    $_ = $seq;

    ## assume fasta line

    ## replace IUPAC ambiguity codes with regex patterns

    ## protect the word boundaries
    s/\\b/\$/g;
    s/\\B/\_/g;

    #M
    s/M/[AC]/g; s/m/[ac]/g;
    #R
    s/R/[AG]/g; s/r/[ag]/g;
    #W
    s/W/[AT]/g; s/w/[at]/g;
    #S
    s/S/[CG]/g; s/s/[cg]/g;
    #Y
    s/Y/[CT]/g; s/y/[ct]/g;
    #K
    s/K/[GT]/g; s/k/[gt]/g;
    #V
    s/V/[ACG]/g; s/v/[acg]/g;
    #H
    s/H/[ACT]/g; s/h/[act]/g;
    #D
    s/D/[AGT]/g; s/d/[agt]/g;
    #B
    s/B/[CGT]/g; s/b/[cgt]/g;
    ##. => not going to support this one
    #s/\./[^ATGCatgc]/g;
    #N/X
    s/[NnXx-]/./g;

    if ($to_uc) {
        $_ = uc($_);
    }


    ## replace word boundary placeholder
    ## !!! AFTER UPPERCASING !!!
    s/\$/\\b/g;
    s/\_/\\B/g;


    return $_;
}


## carried forward from an alternative version of this script;
## not currently used
sub parenthesize {
    my ($s) = @_;

    my @chars = undef;
    my @stack = undef;

    ## from the left
    @chars = (split(//, $s));
    @stack = ();
    while (scalar(@chars) >= $subseq_minlen) {
        push(@stack, shift @chars);
        #    splice(@stack, scalar(@stack), 0, splice(@chars, 0, $step));
    }
    ## @ = \b
    my $left = '@'.join('', @chars);
    while (my $char = shift @stack) {
    #while (my @seg = splice(@stack, 0, $step)) {
        push(@chars, $char);
        $left = join('', @chars).'@'.'|'.$left;
    }
    $left = "(".$left.')';

    @chars = split(//, $s);
    @stack = ();
    while (scalar(@chars) >= $subseq_minlen) {
        push(@stack, pop @chars);
    }
    my $right = join('', @chars).'@';
    while (my $char = shift @stack) {
        unshift(@chars, $char);
        #$right .= '|'.'^'.join('', @chars);
        $right .= '|'.'@'.join('', @chars);
    }
    $right = '('.$right.')';

    return '('.$left.'|'.$right.')';
}

## TODO: generalize to handle any non-numerics
sub chr_to_idx {
    my ($chr) = @_;

    $chr =~ s/^chr//i;
    $chr =~ s/^X$/23/i;
    $chr =~ s/^Y$/24/i;
    $chr =~ s/^M$/25/i;

    ## not sure if HS will choke on this
    if ($zero_pad) {
        return sprintf("%04d", int($chr * 100));
    } else {
        return int(($chr + 10) * 100);
    }
}

sub check_flags {
	my ($flags_opt) = @_;
        
	my $valid = {
		'i' => 'case insensitive',
	 	's' => 'matches newlines',
		'm' => 'multiline anchoring',
		'o' => 'report match ID at most once',
		'e' => 'allow patterns that can match empty buffers',
                'u' => 'UTF-8 mode',
                'p' => 'unicode property support',
                'f' => 'prefiltering mode',
                'l' => 'leftmost start of match reporting',
                'C' => 'logical combination of patterns',
                'Q' => 'quiet at matching',
	};
	## mutually exclusive flags
	my $incompatible = {
		'o'	=>	{ 'l' => 1 },
		'l'	=>	{ 'o' => 1 },
	};
 
	my @flags = split(//, $flags_opt);

	my $is_set = {};
	foreach my $flag(@flags) {
		if (! defined($valid->{$flag})) {
			confess("Provided flag character '$flag' is not valid!");
		}
		if (defined($incompatible->{$flag})) {
			foreach my $icflag(keys %{$incompatible->{$flag}}) {
				if (defined($is_set->{$icflag})) {
					confess("Provided flag character '$flag' is incompatible with '$icflag'!");
				}
			}
		}
						
		$is_set->{$flag} = 1;
	}

	return join('', keys %{$is_set});	
}
