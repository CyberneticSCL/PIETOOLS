#!/usr/bin/perl
# phaseB_tierhints.pl DIR
# Patterns in the example-specific settings that may suggest changes to the
# default tiers (maintainer, 11/19/2026): from DIR/phaseB_results.csv (final
# settings of the knob search per row: form, tier, statements relative to
# lpisettings_sop(tier)) and DIR/row_exec.txt (the executive of each row), per
# kind of executive (stability, gain, synthesis, H2 norm) and dimension:
#   - the tier each tuned setting is built on, and the form;
#   - each change of a knob from its tier value, with the share of rows of
#     the group that carry it;
#   - the reduction factor of the group (start/final decision variables).
# Rows of PDE 33-48 and 56-72 count as two-variable, 53-55 as three.
# MMP, 11/19/2026.
use strict; use warnings;
my ($D) = @ARGV;
my %exe; { open(my $h, '<', "$D/row_exec.txt") or die; while (<$h>) { s/\r?\n$//; my ($k, $e) = split / /; $exe{$k} = $e } }
sub kind { my $e = shift // ''; return 'stability' if $e =~ /stability/; return 'gain' if $e =~ /Hinf_gain/;
  return 'synthesis' if $e =~ /control|estimator/; return 'H2 norm' if $e =~ /H2_norm/; return 'other'; }
sub dim { my $k = shift; return '1-D' unless $k =~ /^PDE:(\d+)/; my $n = $1;
  return '3-D' if $n >= 53 && $n <= 55; return '2-D' if ($n >= 33 && $n <= 48) || ($n >= 56 && $n <= 72); return '1-D'; }
my (%g, %tier, %form, %chg, %fac);
open(my $h, '<', "$D/phaseB_results.csv") or die; <$h>;
while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = ($l =~ /("(?:[^"]|"")*"|[^,]*)(?:,|$)/g); s/^"(.*)"$/$1/ for @x;
  my ($k, $tag, $fm, $t, $stm, $n0, $m0, $n1) = @x[0..7]; next unless $n1 =~ /^\d+$/;
  my $G = kind($exe{$k}) . ' ' . dim($k); $g{$G}++; $tier{$G}{$t}++; (my $F = $fm) =~ s/d$//; $form{$G}{$F || '-'}++;
  push @{$fac{$G}}, $n0 / $n1 if $n0 =~ /^\d+$/ && $n1 > 0;
  for my $s (grep { $_ ne '' } map { s/^\s+|\s+$//gr } split /(?<=;)\s+(?=st\.)/, $stm) { $s =~ s/;$//; $chg{$G}{$s}++; } }
close $h;
for my $G (sort keys %g) { my $n = $g{$G};
  my @f = sort { $a <=> $b } @{$fac{$G} // []};  my $gm = 0; $gm += log($_) for @f;
  printf "\n== %s: %d rows; start/final variables median %.2f, geometric mean %.2f\n", $G, $n,
    (@f ? $f[int(@f/2)] : 0), (@f ? exp($gm/@f) : 0);
  print "   tier: ", join(', ', map { "$_ $tier{$G}{$_}" } sort { $tier{$G}{$b} <=> $tier{$G}{$a} } keys %{$tier{$G}}), "\n";
  print "   form: ", join(', ', map { "$_ $form{$G}{$_}" } sort { $form{$G}{$b} <=> $form{$G}{$a} } keys %{$form{$G}}), "\n";
  my @c = sort { $chg{$G}{$b} <=> $chg{$G}{$a} } keys %{$chg{$G}};
  for my $s (grep { $chg{$G}{$_} / $n >= 0.15 } @c) { printf "   %3.0f%%  %s\n", 100 * $chg{$G}{$s} / $n, $s; } }
