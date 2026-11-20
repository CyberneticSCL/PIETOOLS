#!/usr/bin/perl
# phaseB_patterns.pl HARNESS_DIR DIR
# Which knobs the per-example search could lower, by kind of executive:
# from the build logs exs_<tag>.csv of ex_search.m, every kept change and
# every change that failed (certified no more, or lost accuracy), grouped by
# knob family (storage: d1, d2, d3, patterns, deg_Q, var.deg; slack:
# dD, psatz, cap_inc, wmin, raises, dw_feas, dw_P; cross terms) and by
# executive (stability, H-infinity gain, synthesis, H2 norm; from
# DIR/row_exec.txt), with the decision variables each family removed.
# MMP, 11/19/2026.
use strict; use warnings;
my ($H, $D) = @ARGV;
my %exe; { open(my $h, '<', "$D/row_exec.txt") or die; while (<$h>) { s/\r?\n$//; my ($k, $e) = split / /; $exe{$k} = $e } }
sub kind { my $e = shift // ''; return 'stability' if $e =~ /stability/; return 'gain' if $e =~ /Hinf_gain/;
  return 'synthesis' if $e =~ /control|estimator/; return 'H2 norm' if $e =~ /H2_norm/; return 'other'; }
sub fam { my $m = shift; return 'cross' if $m =~ /^x:/;
  # trade names: 't:' + raise + '+' + reduction; storage raises carry a '+' of
  # their own ('d2+1.1', 'degQ+1.2', 'vdeg+1.3'), slack raises are 'dD+',
  # 'psatz+k', 'wmin+k', 'raiseQ+k', 'dwfeas+k', 'dwP+k'
  if ($m =~ /^t:/) { return $m =~ /^t:(d[123]|degQ|vdeg)\+\d+\.\d+\+/ ? 'trade: storage up, slack down' : 'trade: slack up, storage down' }
  return 'storage' if $m =~ /^(d[123]|pat|degQ|vdeg)/; return 'slack' if $m =~ /^(dD0|psatz|capinc|wmin|raise|dwfeas|dwP)/;
  return 'other'; }
my (%kept, %fail, %gain, %rows, %cur);
for my $f (glob("$H/exs_*.csv")) { next if $f =~ /final/;
  open(my $h, '<', $f) or die; <$h>;
  while (my $l = <$h>) { my @x = split /,/, $l;  my ($k, $mv, $out, $ndv) = ($x[1], $x[3], $x[4], $x[10]);
    my $kd = kind($exe{$k}); $rows{$kd}{$k} = 1;
    if ($out eq 'start') { $cur{$k} = $ndv; next }
    next if $f =~ /TRADE/ && $mv !~ /^t:/;     # the trade pass: its trades only (no double count)
    my $fm = fam($mv);
    if ($out eq 'kept') { $kept{$kd}{$fm}++; $gain{$kd}{$fm} += ($cur{$k} // $ndv) - $ndv; $cur{$k} = $ndv; }
    elsif ($out eq 'fails') { $fail{$kd}{$fm}++; } } close $h; }
printf "%-10s %5s | %-31s %6s %6s %12s\n", 'executive', 'rows', 'family', 'kept', 'failed', 'vars removed';
for my $kd (sort keys %rows) { my $n = scalar keys %{$rows{$kd}};
  for my $fm ('storage', 'slack', 'cross', 'trade: slack up, storage down', 'trade: storage up, slack down') {
    next unless ($kept{$kd}{$fm} // 0) + ($fail{$kd}{$fm} // 0);
    printf "%-10s %5d | %-31s %6d %6d %12d\n", $kd, $n, $fm, $kept{$kd}{$fm} // 0, $fail{$kd}{$fm} // 0, $gain{$kd}{$fm} // 0; } }
# the slack moves kept, by name, per executive
print "\nkept moves by name (executive: move count)\n";
my %byname;
for my $f (glob("$H/exs_*.csv")) { next if $f =~ /final/; open(my $h, '<', $f) or die; <$h>;
  while (my $l = <$h>) { my @x = split /,/, $l; next unless $x[4] eq 'kept'; (my $m = $x[3]) =~ s/\d+(\.\d+)?$//;
    $byname{kind($exe{$x[1]})}{$m}++; } close $h; }
for my $kd (sort keys %byname) { print "  $kd: ", join(', ', map { "$_ $byname{$kd}{$_}" } sort { $byname{$kd}{$b} <=> $byname{$kd}{$a} } keys %{$byname{$kd}}), "\n"; }
