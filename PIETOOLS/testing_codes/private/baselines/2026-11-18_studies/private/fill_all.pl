#!/usr/bin/perl
# fill_all.pl BASE OUT FILE:TIER ...  study E of PLAN_1118_studies.md: copies the merged pass BASE
# (results_all_1118.csv) to OUT, replacing or adding the rows of TIER from each FILE (the gap runs
# E1str, E2l, E3vh; their seeded rows of other tiers are ignored), and prints per tier the counts
# of rows by verdict and the summed seconds. MMP, 11/18/2026.
use strict; use warnings;
my ($B,$O,@G) = @ARGV;
my (@order,%row);
open(my $fh,'<',$B) or die "$B"; my $hdr = <$fh>; $hdr =~ s/\r?\n$//;
while (<$fh>) { s/\r?\n$//; next unless length; my @c = split /,/, $_, 4; my $k = "$c[0],$c[1],$c[2]";
  push @order,$k unless exists $row{$k}; $row{$k} = $_; }
close $fh;
my $nrep = 0; my $nadd = 0;
for my $g (@G) { my ($f,$t) = $g =~ /^(.*):([^:]+)$/;   # the last colon: paths carry a drive letter
  open(my $gh,'<',$f) or die "$f"; <$gh>;
  while (<$gh>) { s/\r?\n$//; next unless length; my @c = split /,/, $_, 4; next unless $c[2] eq $t;
    my $k = "$c[0],$c[1],$c[2]"; if (exists $row{$k}) { $nrep++ } else { push @order,$k; $nadd++ } $row{$k} = $_; }
  close $gh; }
open(my $oh,'>',$O) or die; print $oh "$hdr\n"; print $oh "$row{$_}\n" for @order; close $oh;
print "rows replaced $nrep, added $nadd, total ".scalar(@order)."\n";
my %s;
for my $k (@order) { my @c = split /,/, $row{$k}, 10; my $t = $c[2]; $s{$t}{rows}++; $s{$t}{$c[3]}++;
  $s{$t}{sec} += $c[7] if $c[7] =~ /^[\d.]+$/; }
my @st = qw(certified converged inaccurate infeasible unknown lo_fail hi_pass error skipped crash);
print "| tier | rows | ".join(' | ',@st)." | seconds |\n|".('---|' x (@st+3))."\n";
for my $t (qw(extreme stripped light heavy veryheavy)) { next unless $s{$t};
  print "| $t | $s{$t}{rows} | ".join(' | ', map { $s{$t}{$_} // 0 } @st).sprintf(" | %.0f |\n",$s{$t}{sec}); }
