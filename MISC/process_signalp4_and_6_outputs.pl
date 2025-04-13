#!/usr/bin/perl

# JFR::Fasta is available https://github.com/josephryan/JFR-PerlModules

# $SIG4 is output from SignalP-4.1 (cmd: signalp -n hsym -l sig.log Hsym_primary_v1.0.aa > sig.out

# $SIG6 is output from SignalP-6.0 (cmd: 
# signalp6 --model_dir=signalp6_slow_sequential/signalp-6-package/models --fastafile=Hsym_primary_v1.0.aa --output_dir=outdir2 --organism=eukarya --format=txt --bsize=100 --write_procs=20 --mode=slow-sequential > s6.out 2> s6.err
# perl -ne '@ff = split; print "$ff[0]\n" if ($ff[1] eq 'SP')' outdir2/prediction_results.txt > sp_ids.txt
# )

# mkdir fasta
# perl process_signalp4_and_6_outputs.pl > process_signalp4_and_6_outputs.out

use strict;
use warnings;
use JFR::Fasta;
use Data::Dumper;

our $FILE = '../00-DATA/Hsym_primary_v1.0.aa';
our $SIGP4 = '../01-SIGNALP/sig.out';
our $SIGP6 = '../03-SIGNALP6/outdir2/prediction_results.txt';

MAIN: {
    my $rh_t = get_tms($SIGP4);  
    my $rh_s = get_sigp($SIGP6);
    my $fp = JFR::Fasta->new($FILE);
    my $onesies = 0;
    print "#ID,NUM,NUM_COMPARISONS,SCORE,NORMALIZED_SCORE,AAS_ADJACENT_TO_CUTSITE\n";
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        next if ($rh_t->{$id});
        next unless ($rh_s->{$id});
        next unless ($rec->{'seq'} =~ m/G[KR][KRED]/);
        open OUT, ">fasta/$id.fa" or die "cannot open fasta/$id.fa:$!";
        if ($rec->{'seq'} =~ m/G[KR][KRED].*G[KR][KRED]/) {
            my @bits = split /G[KR][KRED]/, $rec->{'seq'};
            my $last = pop @bits;
            my %seen = ();
            my $score = 0;
            my $num = 0;
            my $num_cmps = 0;
            my @sixlets = ();
            for (my $i = 0; $i < @bits; $i++) {
                next unless ($bits[$i]) =~ m/(......)$/;
                my $i_bit = $1;
                push @sixlets, $i_bit;
                $num++;
                for (my $j = 0; $j < @bits; $j++) {
                    next unless ($bits[$j] =~ m/(......)$/);
                    my $j_bit = $1;
                    next if ($i == $j);
                    my $key  = $i . '-' . $j;
                    my $ikey = $j . '-' . $i;
                    next if ($seen{$key});
                    $score += get_score($i_bit,$j_bit);
                    $seen{$key}++;
                    $seen{$ikey}++;
                    $num_cmps++;
                }
            }
            my $bit_string = join ':', @sixlets;
            my $fa_string = join "\n>", @sixlets;
            print OUT ">$fa_string\n";
            my $denom = ($num * ($num - 1)) / 2;
            my $norm_score = 0;
            $norm_score = $score / $denom if ($denom);
            print "$id,$num,$num_cmps,$score,$norm_score,$bit_string\n";
#            print "$id,$num,$bit_string,$score\n";
        } else {
            $onesies++;
        }
    }
    print "ONESIES: $onesies\n";
}

sub get_sigp {
    my $file = shift;
    my %data = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        next if ($line =~ m/^#/);
        my @ff = split /\s+/, $line;
        $data{$ff[0]} = 1 if ($ff[1] eq 'SP');
    }
    return \%data;
}

sub get_tms {
    my $file = shift;
    my %data = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        next if ($line =~ m/^#/);
        my @ff = split /\s+/, $line;
        $data{$ff[0]} = 1 if ($ff[11] eq 'SignalP-TM');
    }
    return \%data;
}

sub get_score {
    my ($seq1, $seq2) = @_;
    my $matches = 0;
    for (my $i = 0; $i < length($seq1); $i++) {
        $matches++ if substr($seq1, $i, 1) eq substr($seq2, $i, 1);
    }
    return $matches;
}
