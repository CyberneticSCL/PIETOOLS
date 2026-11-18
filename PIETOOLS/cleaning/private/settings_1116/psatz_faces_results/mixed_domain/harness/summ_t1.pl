#!/usr/bin/perl
# summ_t1.pl LOG: a markdown table of the T1 reach run: one row per system,
# one column per term set (names from the '@@ T1' line), entries the reach
# gam_min/lam or the failing status at 1.1 lam. SCRATCH, MMP, 11/16/2026.
use strict; use warnings;
my ($log) = @ARGV;
open(my $fh, '<', $log) or die $!;
my (@names, %row, @sys, %lam);
my $cur = '';  my $i = 0;
while (my $l = <$fh>) {
    if ($l =~ /^\@\@ T1 (\S+) terms=(\S+)/) { $cur = $1; @names = split /,/, $2; $i = 0; push @sys, $cur; next; }
    if ($l =~ /^\@\@ POINC (\S+) lam=(\S+)/) { $lam{$1} = $2; next; }
    if ($l =~ /^\@\@ REACH (\S+) .*reach=(\S+) last_fail=(\S*) solves=(\d+)/) {
        my ($s, $r, $f) = ($1, $2, $3);
        my $nm = $names[$i] // "set$i";  $i++;
        $row{$s}{$nm} = ($r eq 'NaN') ? "$f at 1.1" : $r;
    }
}
close $fh;
print "| PIE | lam | ", join(" | ", @names), " |\n";
print "|---|---|", join("", map { "---|" } @names), "\n";
for my $s (@sys) {
    print "| $s | ", ($lam{$s} // ''), " | ", join(" | ", map { $row{$s}{$_} // '' } @names), " |\n";
}
