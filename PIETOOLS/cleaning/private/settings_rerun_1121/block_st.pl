#!/usr/bin/perl
# block_st.pl ID ... - for each DDE/NDS/DDF block id of the settings rerun
# (11/21/2026), the st and form lines the CURRENT library block runs. The
# harness copy blocks/<id>.m (lb_blocks, made before the 11/19-11/20 edits
# of the libraries) supplies the system lines; the block is located in the
# library by its first three system lines in order (within 15 lines), and
# its st and form lines are read up to the next block separator.
# MMP, 11/21/2026.
use strict; use warnings;
my $root = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
my $hb = "$root/testing_codes/private/baselines/2026-11-18_studies/private/blocks";
my %lib = (DDE => "$root/PIETOOLS_examples/examples_DDE_library_PIETOOLS.m",
           NDS => "$root/PIETOOLS_examples/examples_NDSDDF_library_PIETOOLS.m",
           DDF => "$root/PIETOOLS_examples/examples_NDSDDF_library_PIETOOLS.m");
sub norm { my ($s) = @_; $s =~ s/\r//; $s =~ s/^[\s%]+//; $s =~ s/%.*$//; $s =~ s/\s+//g; return $s; }
for my $id (@ARGV) {
    my ($pre) = ($id =~ /^([A-Z]+)/);
    open(my $b, '<', "$hb/$id.m") or do { print "$id: no harness block\n"; next; };
    my @sys;
    while (<$b>) { my $n = norm($_); push @sys, $n if $n =~ /^(DDE|DDF|NDS)\.(A0|Ai|tau|B1|C1|A)\b/ && !grep { $_ eq $n } @sys; last if @sys == 3; }
    close $b;
    open(my $l, '<', $lib{$pre}) or die; my @L = <$l>; close $l;
    my $at = -1;
    for my $i (0..$#L) {
        next unless norm($L[$i]) eq $sys[0];
        my ($k, $j) = (1, $i+1);
        while ($k < @sys && $j <= $#L && $j < $i+15) { $k++ if norm($L[$j]) eq $sys[$k]; $j++; }
        if ($k == @sys) { $at = $i; last; }
    }
    if ($at < 0) { print "$id: block not found by [@sys]\n"; next; }
    my @out;
    for my $j ($at..$#L) {
        last if $j > $at && $L[$j] =~ /^\s*%\s*%*\s*%{20,}/;
        my $t = $L[$j]; $t =~ s/\r?\n$//; $t =~ s/\s+% MMP, \d\d\/\d\d\/\d{4}( \(was\))?\s*$//;
        push @out, $t if $t =~ /^\s*%?\s*st\s*=|^\s*%?\s*st\.|_form\s*=|^\s*%?\s*(Hinf_gain|Hinf_gain_dual|stability|stability_dual|H2_norm|H2_norm_dual|Hinf_control|Hinf_estimator)\s*=\s*1/;
    }
    printf "%s (line %d): %s\n", $id, $at+1, join(' | ', map { s/^\s+//r } @out);
}
