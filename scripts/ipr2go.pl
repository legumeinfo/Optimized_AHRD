#!/usr/bin/env perl

use strict;
use English;
$OFS="\t";
$ORS="\n";
my $score=0.5;

my %p2go;
while (<>) {
    chomp;
    my @data = split /\t/;
    if ($data[$#data] !~ /NULL/) {
        my $p = $data[0];
        my @go = split /,/, $data[$#data];
        foreach my $go (@go) {
            my ($aspect, $id, $name) = ($go =~ /(Cellular Component|Molecular Function|Biological Process):(.*) \((GO:[0-9]+)\)/);
            if (defined $id) {
                $p2go{$p}->{$id} = $name;
            }
        }
    }
}

foreach my $p (keys %p2go) {
    foreach my $id (keys %{$p2go{$p}}) {
        print $p, $p2go{$p}->{$id}, $id;
    }
}
