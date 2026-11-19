#!/usr/bin/perl
# b5_list.pl BASE_CSV LOGDIR VARIANT: every 1-D row of the variant (ab_slack logs) against the
# base light rows (results_L1ctl.csv, the same tree), with the verdict, the value, the variables,
# and the B5 class: 'cand' = same verdict, value equal within 1e-4 relative or better in the
# row's direction (reach: higher; bound: lower), fewer variables. MMP, 11/18/2026.
use strict; use warnings;
my ($B,$L,$v) = @ARGV;
my %base;
open(my $fh,'<',$B) or die "$B"; <$fh>;
while (<$fh>) { chomp; s/\r$//; my @c = split /,/, $_, 10; next if @c<10;
  $base{"$c[0]:$c[1]"} = {st=>$c[3], res=>$c[4], ndv=>$c[8]}; }
close $fh;
my $cert = { certified=>1, converged=>1 };
for my $lib ('DDE','NDS','PDE1D') {
  open(my $g,'<',"$L/ab_slack_${v}_$lib.log") or next;
  while (<$g>) { chomp; s/\r$//;
    next unless /^\@\@AB $v\s+(\S+)\s+(\S+)\s+result\s+(\S+)\s+n=(\d+)\s+ndv=(\d+)/;
    my ($key,$st,$res,$n,$ndv) = ($1,$2,$3,$4,$5);
    my $b = $base{$key}; next unless $b;
    my $cls = 'other';
    my $bc = $cert->{$b->{st}} || 0; my $nc = $cert->{$st} || 0;
    my $num = ($b->{res} =~ /^[-\d.eE+]+$/ && $res =~ /^[-\d.eE+]+$/);
    my $stab = ($b->{st} eq 'certified' && $b->{res} eq '1' && $res eq 'NaN');
    if ($bc && $nc && $b->{st} eq $st) {
      my $dir = ($st eq 'converged') ? 1 : -1;                      # reach up, bound down
      my $rel = ($num) ? ($res - $b->{res}) / (abs($b->{res}) || 1) : 0;
      my $ok = $stab || !$num || $dir*$rel >= -1e-4;
      $cls = ($ok && $ndv < $b->{ndv}) ? 'cand' : ($ok ? 'same-more' : 'value-lost');
    } elsif ($nc && !$bc) { $cls = 'gained' } elsif ($bc && !$nc) { $cls = 'lost' }
    printf "%-24s %-10s %-11s %-14s -> %-11s %-14s ndv %7d -> %7d (%.2f)\n", $key, $cls, $b->{st}, $b->{res}, $st, $res, $b->{ndv}, $ndv, ($b->{ndv} =~ /^\d+$/ && $b->{ndv}>0) ? $ndv/$b->{ndv} : 0;
  }
  close $g;
}
