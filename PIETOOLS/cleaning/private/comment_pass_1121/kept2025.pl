#!/usr/bin/perl
# kept2025.pl FILE... - for each file of the pass, the lines that HEAD shares with the 2025
# version (lines of the prior release) and that the working copy no longer has: the pass keeps
# every line of the 2025 release in a condensed file. Prints MISSING lines, or "KEPT  FILE".
use strict; use warnings;
my $REPO = '/c/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS';
my $S    = '/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/cpass1121';
my %p25;
open(my $c, '<', "$S/census.tsv") or die; while (<$c>) { chomp; my @a = split /\t/; $p25{$a[0]} = $a[1] } close $c;
for my $f (@ARGV) {
    my @old  = map { s/\r?\n$//r } `git -C $REPO show "origin/PIETOOLS-2025:PIETOOLS/$p25{$f}"`;
    my @head = map { s/\r?\n$//r } `git -C $REPO show "HEAD:PIETOOLS/$f"`;
    my @work = do { open(my $h, '<:raw', "$REPO/PIETOOLS/$f") or die; map { s/\r?\n$//r } <$h> };
    my %R = map { $_ => 1 } @old;
    my (%need, %have);
    $need{$_}++ for grep { $R{$_} && /\S/ } @head;
    $have{$_}++ for @work;
    my @miss = grep { ($have{$_} // 0) < $need{$_} } sort keys %need;
    if (@miss) { print "MISSING in $f:\n", map { "    $_\n" } @miss } else { print "KEPT  $f\n" }
}
