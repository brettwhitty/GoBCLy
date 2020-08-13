#!/usr/bin/env perl

use strict;
use warnings;

$| = 1;

## Brett Whitty <brettwhitty@gmail.com>, all rights retained

my $patterns = {};
my $pflags = {};
my $pextflags = {};

while (<>) {
    chomp;

    if (/^(\d+):\/([^\/]+)\/(\S+)$/) {
        my ($id, $pattern, $flags) = ($1, $2, $3);

	unless (defined($patterns->{$id})) {
		$patterns->{$id} = [];
	}
        push (@{$patterns->{$id}}, $pattern);

        if ($flags =~ /\{([^\}]+)\}/) {
            my $extflags = $1;
            foreach my $i(split(/,/, $extflags)) {
                $pextflags->{$id}->{$i} = 1;
            }
            $flags =~ s/\{([^\}]+)\}//;
        }
        foreach my $i(split(//, $flags)) {
            $pflags->{$id}->{$i} = 1;
        }
    } else {
        print "$_\n";
    }
}
foreach my $id(sort keys %{$patterns}) {
    my $pattern = $id
        .':'
        .'/('
        .join('|', @{$patterns->{$id}})
        .')/'
        .join('', keys %{$pflags->{$id}})
        .'{'
        .join(',', keys %{$pextflags->{$id}})
        .'}';
    ## lazy handling of empty extended flags
    $pattern =~ s/\{\}$//;
    print $pattern."\n";
}
