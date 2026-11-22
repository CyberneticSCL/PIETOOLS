#!/usr/bin/perl
# report_table.pl TAG STARTS [TAG STARTS ...] - the result table of the
# settings rerun of 11/21/2026 (PLAN.md) from exs_<TAG>_final.csv in the
# worktree harness folder: one row per example and executive, named by
# example file and library index (maintainer, 11/21: no bare problem
# numbers), with the start, the final settings and the path that accepted
# them (direct: the gate of 11/20; post: the fixed-storage check).
# MMP, 11/21/2026.
use strict; use warnings;
my $H = 'C:/wt1121r/PIETOOLS/testing_codes/private/baselines/2026-11-18_studies/private';
my $CK = "C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS/cleaning/private/settings_rerun_1121/catalog_keys.csv";
sub csv { my ($l) = @_; $l =~ s/\r?\n$//; my @v = ($l =~ /("[^"]*"|[^,]*)(?:,|$)/g); pop @v if @v && $v[-1] eq '' && $l !~ /,$/; s/^"|"$//g for @v; return @v; }
my (%file, %exec);                         # catalog key -> example file, executive (catalog_keys.csv)
open(my $fl, "<", $CK) or die "$CK: $!"; <$fl>;
while (<$fl>) { my @v = csv($_); $file{$v[0]} = $v[3]; ($exec{$v[0]} = $v[4]) =~ s/^PIETOOLS_//; }
close $fl;
sub num { my ($x) = @_; return (defined $x && $x =~ /^-?[0-9.eE+]+$/ && $x ne 'NaN') ? $x : undef; }
sub fmt { my ($x) = @_; my $n = num($x); return '–' unless defined $n; 1 while $n =~ s/^(\d+)(\d{3})/$1,$2/; return $n; }
print "| Example | Index | Executive, form | Start | Variables: start -> final | Path | Smallest direct | Final st (from the tier) |\n";
print "|---|---|---|---|---|---|---|---|\n";
while (@ARGV) {
    my ($tag, $starts) = (shift, shift);
    my %start;
    open(my $fs, '<', $starts) or die "$starts: $!"; <$fs>;
    while (<$fs>) { my @v = csv($_); $start{$v[0]} //= [@v]; }
    close $fs;
    # the start that certified, from the per-build log (outcome 'start'):
    # the starts file lists every alternative, the first is not always it
    my %cert;
    if (open(my $fb, '<', "$H/exs_${tag}.csv")) { <$fb>;
        while (<$fb>) { my @v = csv($_); next unless defined $v[4] && $v[4] eq 'start';
            my ($f, $t, $m) = ($v[3] =~ /^start\d+ (\S+) (\S+) ?(.*)$/); $cert{$v[1]} = [$f, $t, $m] if defined $f; }
        close $fb; }
    my $ff = "$H/exs_${tag}_final.csv";
    open(my $fh, '<', $ff) or die "$ff: $!"; <$fh>;
    while (<$fh>) {
        my ($key,$form,$tier,$stm,$n0,$m0,$n,$m,$builds,$reason,$path,$nd,$md) = csv($_);
        my ($lib,$id,$ex) = ($key =~ m{^([A-Z]+):([^/]+)(?:/(.+))?$});
        my ($exf, $idx) = ($file{$key} // "?", $id);
        $exf =~ s/\.m$//;
        my $exn = $exec{$key} // "?";
        # the form reaches the stability, gain and H2-norm executives only
        # (lb_exec calls the estimator and controller executives without one)
        $exn .= ", $form" unless $exn =~ /control|estimator/;
        my $s = $start{$key};
        my $sd = $s ? sprintf('%s %s', $s->[3], ($s->[4] eq 'examplest' ? 'example st' : ($s->[4] eq '' ? 'tier' : 'listed mods'))) : '?';
        if (my $c = $cert{$key}) {
            $sd = sprintf('%s %s %s', $c->[0], $c->[1], ($c->[2] eq 'examplest' ? 'example st' : ($c->[2] eq '' ? 'tier' : 'listed mods')));
        }
        my $var = sprintf('%s -> %s', fmt($n0), fmt($n));
        $var .= " ($reason)" if $reason ne 'converged';
        my $pth = $path eq '' ? 'none' : $path;
        my $dir = (num($nd) ? fmt($nd) : '–');
        $stm =~ s/\|/\\|/g;
        print "| $exf | $idx | $exn | $sd | $var | $pth | $dir | `$stm` |\n";
    }
    close $fh;
}
