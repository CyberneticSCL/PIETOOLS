#!/usr/bin/perl
# cost_at_fixed_accuracy.pl HARNESS_DIR BASE_1028 DIR
# Complexity at a fixed accuracy, stock against present, without new runs.
# Per row of the 10/28 baseline, a common target: the best accuracy that both
# the stock tiers and the present tiers reach (stability rows: a certificate;
# bisection rows: the lower of the two best thresholds; value rows: the higher
# of the two best gains or norms; within 1e-3 relative). For each family the
# cost to reach the target is the fewest decision variables (cone size nx)
# over its measurements that reach it:
#   stock   - the six stock tiers of 10/28/2026 (stock_results.csv);
#   tiers   - the five present tiers in their default forms (current tree;
#             the 11/17 pass where a tier was not rerun);
#   tuned   - every present measurement: tiers, the variants of the 11/18
#             studies and of Phase B, and the final settings of the knob
#             search (DIR/phaseB_results.csv), which reach the present best.
# Harness gate on every side. Prints the ratios to stock by kind of row and
# number of variables. MMP, 11/19/2026.
use strict; use warnings;
my ($H, $B, $D) = @ARGV;  my $RT = 1e-3;
sub isnum { my $v = shift; return defined $v && $v =~ /^-?[0-9.]+(?:[eE][-+]?\d+)?$/ }
sub bis { $_[0] =~ /^(converged|max_solves|lo_fail|hi_pass)$/ }
sub rdall { my ($f, $tiers) = @_; my @r; open(my $h, '<', $f) or return (); <$h>;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = split /,/, $l, 10; next if @x < 10;
    next if $tiers && !$tiers->{$x[2]}; next if $x[3] eq 'skipped';
    push @r, { key => "$x[0]:$x[1]", tier => $x[2], st => $x[3], res => $x[4], ndv => $x[8], note => $x[9] }; }
  close $h; return @r; }
# measurement -> (certified, accuracy value or undef for verdict rows)
sub acc { my ($e) = @_;
  if (bis($e->{st})) { return ($e->{st} ne 'lo_fail' && isnum($e->{res})) ? (1, $e->{res}) : (0) }
  return (0) unless $e->{st} eq 'certified';
  return ($e->{note} =~ /Hinf_|H2_/ && isnum($e->{res})) ? (1, $e->{res}) : (1, undef); }
my %fam;
push @{$fam{stock}{$_->{key}}}, $_ for rdall("$B/stock_results.csv");
my @tierfiles = map { "$H/results_$_.csv" } qw(M1L H1ctl XSQ C1ext B4ctl B1c0 X2SQ E1str C1ext2 E3vh);
my %seen;
for my $f (@tierfiles) { for my $e (rdall($f)) { next if $f =~ /H2ctl/ && $e->{tier} eq 'light'; push @{$fam{tiers}{$e->{key}}}, $e; $seen{"$e->{key}|$e->{tier}"} = 1 } }
for my $e (rdall("$H/results_H2ctl.csv", { heavy => 1 })) { push @{$fam{tiers}{$e->{key}}}, $e; $seen{"$e->{key}|heavy"} = 1 }
for my $e (rdall("$H/results_all_1118b.csv")) { next if $seen{"$e->{key}|$e->{tier}"}; push @{$fam{tiers}{$e->{key}}}, $e }
# tuned: everything measured on the present side
for my $k (keys %{$fam{tiers}}) { push @{$fam{tuned}{$k}}, @{$fam{tiers}{$k}} }
for my $f (glob("$H/results_*.csv")) { next if $f =~ /results_all_1118b|_started/;
  for my $e (rdall($f)) { push @{$fam{tuned}{$e->{key}}}, $e } }
# the knob-search finals: certified at the present best (by construction)
my %final; { open(my $h, '<', "$D/phaseB_results.csv") or die; <$h>;
  while (my $l = <$h>) { my @x = ($l =~ /("(?:[^"]|"")*"|[^,]*)(?:,|$)/g); $final{$x[0]} = $x[7] if isnum($x[7]); } close $h; }
my (%R, %only);
for my $k (sort keys %{$fam{stock}}) {
  next unless $fam{tiers}{$k};
  my @s = grep { (acc($_))[0] } @{$fam{stock}{$k}};  my @p = grep { (acc($_))[0] } @{$fam{tiers}{$k}};
  if (!@s || !@p) { $only{@s ? 'stock' : @p ? 'tiers' : 'neither'}++; next }
  my ($kind, $A);
  my @sv = map { (acc($_))[1] } @s;  my @pv = map { (acc($_))[1] } @p;
  if (bis($s[0]{st}) || bis($p[0]{st})) { $kind = 'threshold';
    my ($sb) = sort { $b <=> $a } grep { defined } @sv; my ($pb) = sort { $b <=> $a } grep { defined } @pv;
    next unless defined $sb && defined $pb; $A = $sb < $pb ? $sb : $pb; }
  elsif (grep { defined } @sv) { $kind = 'gain/norm';
    my ($sb) = sort { $a <=> $b } grep { defined } @sv; my ($pb) = sort { $a <=> $b } grep { defined } @pv;
    next unless defined $sb && defined $pb; $A = $sb > $pb ? $sb : $pb; }
  else { $kind = 'stability' }
  my $reach = sub { my $e = shift; my ($c, $v) = acc($e); return 0 unless $c && isnum($e->{ndv});
    return 1 if $kind eq 'stability'; return 0 unless defined $v;
    return $kind eq 'threshold' ? $v >= $A - $RT*abs($A) : $v <= $A + $RT*abs($A); };
  my %cost;
  for my $fm ('stock', 'tiers', 'tuned') { my @c = map { $_->{ndv} } grep { $reach->($_) } @{$fam{$fm}{$k} // []};
    ($cost{$fm}) = sort { $a <=> $b } @c; }
  $cost{tuned} = $final{$k} if isnum($final{$k}) && (!defined $cost{tuned} || $final{$k} < $cost{tuned});
  next unless $cost{stock} && $cost{tiers};
  my ($id) = $k =~ /:(\d+)/;  my $dim = ($k =~ /^PDE/ && $id >= 33 && $id <= 48) ? '2-D' : '1-D';
  push @{$R{"$kind|$dim"}}, [$k, $cost{stock}, $cost{tiers}, $cost{tuned} // $cost{tiers}];
}
sub stats { my @v = sort { $a <=> $b } @_; my $n = @v; return ('', '', '', '') unless $n;
  my $g = 0; $g += log($_) for @v; return (sprintf("%.2f", $v[int($n/2)]), sprintf("%.2f-%.2f", $v[int($n/4)], $v[int(3*$n/4)]),
    sprintf("%.2f", exp($g/$n)), scalar(grep { $_ < 1 } @v)); }
printf "%-10s %-4s %4s | %-30s | %-30s\n", 'kind', 'dim', 'rows', 'present tiers / stock: median (IQR) geo, <1', 'tuned / stock: median (IQR) geo, <1';
for my $g (sort keys %R) { my ($kind, $dim) = split /\|/, $g; my @r = @{$R{$g}};
  my @a = stats(map { $_->[2] / $_->[1] } @r); my @b = stats(map { $_->[3] / $_->[1] } @r);
  printf "%-10s %-4s %4d | %5s (%s) %5s, %3d | %5s (%s) %5s, %3d\n", $kind, $dim, scalar(@r), @a, @b; }
open(my $o, '>', "$D/cost_at_fixed_accuracy.csv") or die;
print $o "key,kind,dim,stock,tiers,tuned\n";
for my $g (sort keys %R) { my ($kind, $dim) = split /\|/, $g;
  print $o join(',', $_->[0], $kind, $dim, @$_[1..3]), "\n" for @{$R{$g}}; }
close $o;
my @all = map { @$_ } values %R;
my @a = stats(map { $_->[2] / $_->[1] } @all); my @b = stats(map { $_->[3] / $_->[1] } @all);
printf "%-10s %-4s %4d | %5s (%s) %5s, %3d | %5s (%s) %5s, %3d\n", 'all', '', scalar(@all), @a, @b;
print "rows certified by one family only: ", join(', ', map { "$_ $only{$_}" } sort keys %only), "\n";
