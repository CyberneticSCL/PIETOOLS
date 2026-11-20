#!/usr/bin/perl
# phaseB_collect.pl HARNESS_DIR OUT_PREFIX
# Merges the knob-search results exs_<tag>_final.csv of ex_search.m (tags
# EXT, EXS1, EXS2, EXU) into one table: per row the form, the tier, the st
# statements from lpisettings_sop(tier), the decision variables and m of
# the start and of the result, the builds and the reason the search ended.
# A row in several files keeps its last finished entry (EXU after EXS*).
# Writes OUT_PREFIX.csv and OUT_PREFIX.md. MMP, 11/19/2026.
use strict; use warnings;
my ($H, $OUT) = @ARGV;
my (%r, @order);
# EXV1: the uncertified rows of 02:21-02:55 ran under that name (st_run lost
# its 'EXU' case); TRADE*: the trade pass from the final settings. A later
# file gives the final settings of a row; the first gives its start size.
for my $tag ('EXT', 'EXS1', 'EXS2', 'EXU', 'EXV1', 'TRADE1', 'TRADE2', 'EXX', 'EXU3D', 'EXU3D2', 'EXU3D3', 'EXU2DU', 'EXS2X', 'EXS2Y', 'EXS53', 'EXS54') {
  my $f = "$H/exs_${tag}_final.csv";  next unless -f $f;
  open(my $h, '<', $f) or die; <$h>;
  while (my $l = <$h>) { $l =~ s/\r?\n$//;
    # key,form,tier,"statements",ndv0,m0,ndv,m,builds,"reason"
    my @x = ($l =~ /("(?:[^"]|"")*"|[^,]*)(?:,|$)/g);
    s/^"(.*)"$/$1/ for @x;
    my ($key, $form, $tier, $stm, $n0, $m0, $n1, $m1, $nb, $why) = @x[0..9];
    if (exists $r{$key} && $n1 !~ /^\d+$/) { next }      # a later pass without a result
    my ($s0, $sm0, $stg) = exists $r{$key} && $r{$key}{n0} =~ /^\d+$/ ? ($r{$key}{n0}, $r{$key}{m0}, "$r{$key}{tag}+$tag") : ($n0, $m0, $tag);
    push @order, $key unless exists $r{$key};
    $r{$key} = { tag=>$stg, form=>$form, tier=>$tier, stm=>$stm, n0=>$s0, m0=>$sm0,
                 n1=>$n1, m1=>$m1, nb=>$nb, why=>$why }; }
  close $h; }
open(my $csv, '>', "$OUT.csv") or die;  open(my $md, '>', "$OUT.md") or die;
print $csv "key,search,form,tier,statements,ndv_start,m_start,ndv,m,factor,builds,reason\n";
print $md "| row | form | tier | st (from lpisettings_sop(tier)) | variables start -> final | factor | reason |\n|---|---|---|---|---|---|---|\n";
my ($nr, $nred, $sum0, $sum1) = (0, 0, 0, 0);
for my $k (@order) { my $e = $r{$k};
  my $fac = ($e->{n0} =~ /^\d+$/ && $e->{n1} =~ /^\d+$/ && $e->{n1} > 0) ? sprintf("%.2f", $e->{n0}/$e->{n1}) : '';
  if ($fac ne '') { $nr++; $nred++ if $e->{n1} < $e->{n0}; $sum0 += $e->{n0}; $sum1 += $e->{n1}; }
  (my $s = $e->{stm}) =~ s/"/""/g;
  print $csv join(',', $k, $e->{tag}, $e->{form}, $e->{tier}, "\"$s\"", $e->{n0}, $e->{m0}, $e->{n1}, $e->{m1},
    $fac, $e->{nb}, "\"$e->{why}\""), "\n";
  print $md "| $k | $e->{form} | $e->{tier} | `$e->{stm}` | $e->{n0} -> $e->{n1} | $fac | $e->{why} |\n"; }
close $csv;
printf $md "\nRows with a result: %d; reduced: %d; summed variables %d -> %d.\n", $nr, $nred, $sum0, $sum1;
close $md;
printf STDERR "rows %d, reduced %d, summed variables %d -> %d\n", $nr, $nred, $sum0, $sum1;
