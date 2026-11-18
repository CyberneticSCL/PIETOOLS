#!/usr/bin/perl
# lb_compare.pl NEW.csv [NEW2.csv ...] -- compares the 11/17 tier-pass rows with the
# 10/28 sopvar-path reference (sop_results.csv): per (library, id, settings) the
# status, result, decision variables, m (from the note) and seconds. Prints a
# summary per tier and the rows that differ (status, size, or value beyond 1e-6
# relative), plus the rows new to the 11/17 pass. MMP, 11/17/2026.
use strict; use warnings;
my $ref = $ENV{LB_REF} // 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS/testing_codes/private/baselines/2026-10-28_example_library/sop_results.csv';
my %R; my %N;
sub load { my ($f,$h) = @_; open(my $fh,'<',$f) or die "$f: $!"; my $hdr = <$fh>;
  while (<$fh>) { chomp; s/\r$//; my @c = split /,/, $_, 10; next if @c < 10;
    my ($lib,$id,$set,$st,$res,$doc,$ns,$sec,$ndv,$note) = @c;
    my ($m) = $note =~ /\bm=(\d+)/; $m //= 'NaN';
    $h->{"$lib|$id|$set"} = {st=>$st,res=>$res,ns=>$ns,sec=>$sec,ndv=>$ndv,m=>$m,note=>$note,lib=>$lib,id=>$id,set=>$set}; }
  close $fh; }
load($ref,\%R);
load($_,\%N) for @ARGV;
my %sum; my @diff; my @new;
for my $k (sort keys %N) { my $n = $N{$k}; my $t = $n->{set};
  $sum{$t}{rows}++; $sum{$t}{$n->{st}}++; $sum{$t}{sec} += ($n->{sec} =~ /^[\d.]+$/ ? $n->{sec} : 0);
  if (!exists $R{$k}) { push @new, $n; $sum{$t}{new}++; next; }
  my $r = $R{$k}; my @why;
  push @why, "status $r->{st}->$n->{st}" if $r->{st} ne $n->{st};
  push @why, "ndv $r->{ndv}->$n->{ndv}" if $r->{ndv} ne $n->{ndv};
  push @why, "m $r->{m}->$n->{m}" if $r->{m} ne $n->{m};
  if ($r->{res} =~ /^[-\d.eE+]+$/ && $n->{res} =~ /^[-\d.eE+]+$/) {
    my $d = abs($r->{res}-$n->{res}); my $s = abs($r->{res}) > 0 ? abs($r->{res}) : 1;
    push @why, sprintf("value %s->%s (%.1e rel)",$r->{res},$n->{res},$d/$s) if $d/$s > 1e-6; }
  elsif ($r->{res} ne $n->{res}) { push @why, "value $r->{res}->$n->{res}"; }
  if (@why) { push @diff, [$n,$r,join('; ',@why)]; $sum{$t}{diff}++; } else { $sum{$t}{same}++; }
}
print "## Summary per tier (11/17 pass vs 10/28 sopvar reference)\n| tier | rows | same | diff | new | certified | infeasible | unknown | error | skipped | seconds |\n|---|---|---|---|---|---|---|---|---|---|---|\n";
for my $t (sort keys %sum) { my $s = $sum{$t};
  printf "| %s | %d | %d | %d | %d | %d | %d | %d | %d | %d | %.0f |\n", $t, $s->{rows}//0, $s->{same}//0, $s->{diff}//0, $s->{new}//0,
    ($s->{certified}//0)+($s->{converged}//0), $s->{infeasible}//0, ($s->{unknown}//0)+($s->{inaccurate}//0), $s->{error}//0, $s->{skipped}//0, $s->{sec}//0; }
print "\n## Rows that differ from the reference\n| library | id | tier | reference | 11/17 | difference |\n|---|---|---|---|---|---|\n";
for my $d (@diff) { my ($n,$r,$w) = @$d;
  printf "| %s | %s | %s | %s %s %s/%s %ss | %s %s %s/%s %ss | %s |\n", $n->{lib},$n->{id},$n->{set},$r->{st},$r->{res},$r->{ndv},$r->{m},$r->{sec},$n->{st},$n->{res},$n->{ndv},$n->{m},$n->{sec},$w; }
print "\n## Rows new in the 11/17 pass\n| library | id | tier | status | result | ndv | m | seconds | note |\n|---|---|---|---|---|---|---|---|---|\n";
for my $n (@new) { my $nt = $n->{note}; $nt =~ s/\|/\//g; $nt = substr($nt,0,140);
  printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", $n->{lib},$n->{id},$n->{set},$n->{st},$n->{res},$n->{ndv},$n->{m},$n->{sec},$nt; }
