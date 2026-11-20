#!/usr/bin/perl
# phaseC_plan.pl DIR
# Per example (PDE file and index, or DDE/NDS/DDF library block), the
# settings st and the forms to write, from the knob-search results of
# Phase B (DIR/phaseB_results.csv), the choices (DIR/phaseB_choices.csv), the
# executive of each row (DIR/row_exec.txt) and the PDE case-to-file map
# (DIR/case_file.txt). One st serves every executive of an example (the
# auto-executor passes one settings struct): when the rows of an example
# end at different settings, the row with the most decision variables
# gives st and the others are flagged for verification. Writes
# DIR/phaseC_plan.csv. MMP, 11/19/2026.
use strict; use warnings;
my ($D) = @ARGV;
sub csv { my ($f) = @_; my @r; open(my $h, '<', $f) or return (); my $hd = <$h>; $hd =~ s/\r?\n$//;
  my @c = split /,/, $hd;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = ($l =~ /("(?:[^"]|"")*"|[^,]*)(?:,|$)/g);
    s/^"(.*)"$/$1/s for @x; my %e; @e{@c} = @x; push @r, \%e; } return @r; }
my %res = map { $_->{key} => $_ } csv("$D/phaseB_results.csv");
my %cho = map { ("$_->{library}:$_->{id}") => $_ } csv("$D/phaseB_choices.csv");
my %exe; { open(my $h, '<', "$D/row_exec.txt") or die; while (<$h>) { s/\r?\n$//; my ($k, $e) = split / /; $exe{$k} = $e; } }
my @cases; { open(my $h, '<', "$D/case_file.txt") or die;
  while (<$h>) { s/\r?\n$//; my ($c, $f, $ix) = split / \| /; $c =~ s/[{}]//g;
    push @cases, { nums => [split /,/, $c], file => $f, idx => ($ix // '') }; } }
my %xfile = (X1 => 'PIETOOLS_PDE_Ex_Reaction_Diffusion_Robin', X2 => 'PIETOOLS_PDE_Example_Diagne',
             X3 => 'PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D');     # outside the library (11/19)
sub pde_group { my ($id) = @_; my $n = $id; $n =~ s{/.*}{};
  return ($xfile{$n}, '') if exists $xfile{$n};
  for my $c (@cases) { for my $v (@{$c->{nums}}) { if (abs($v - $n) < 1e-9) {
    my $iv = $c->{idx} ? ($n == int($n) ? sprintf("%.1f", $n) : $n) : '';
    return ($c->{file}, $iv); } } }
  return ('?', ''); }
my %grp;
for my $k (sort keys %exe) {
  my ($lib, $id) = split /:/, $k, 2;
  my ($g, $file, $iv);
  if ($lib eq 'PDE') { ($file, $iv) = pde_group($id); $g = $iv ne '' ? "$file\@$iv" : $file; }
  else { (my $b = $id) =~ s{/.*}{}; $g = "$lib:$b"; $file = $lib eq 'DDE' ? 'examples_DDE_library_PIETOOLS' : 'examples_NDSDDF_library_PIETOOLS'; $iv = $b; }
  push @{$grp{$g}{rows}}, $k; $grp{$g}{file} = $file; $grp{$g}{iv} = $iv; }
# verified overrides (DIR/phaseC_overrides.csv): one measured setting that
# reaches every row of the example, where the merged st failed verification
my %ovr = map { $_->{group} => $_ } csv("$D/phaseC_overrides.csv");
open(my $o, '>', "$D/phaseC_plan.csv") or die;
print $o "group,file,index,rows,executives,stability_form,Hinf_gain_form,tier,statements,ndv,status,notes\n";
my %cnt;
for my $g (sort keys %grp) {
  my @rows = @{$grp{$g}{rows}}; my (@ex, %sf, %gf, @have, @miss);
  for my $k (@rows) { my $e = $exe{$k}; push @ex, $e;
    my $r = $res{$k}; my $f = $r ? $r->{form} : ($cho{$k} ? $cho{$k}{form} : '');
    (my $L = $f) =~ s/d$//;
    if ($e =~ /stability/ && $L =~ /^[PQ]$/) { $sf{$L}++ }
    if ($e =~ /Hinf_gain/ && $L =~ /^[PQ]$/) { $gf{$L}++ }
    if ($r && $r->{ndv} =~ /^\d+$/) { push @have, $r } else { push @miss, $k } }
  my @notes; my ($tier, $stm, $ndv, $status) = ('', '', '', '');
  if (@have) {
    my %u = map { ("$_->{tier}|$_->{statements}") => 1 } @have;
    my ($top) = sort { $b->{ndv} <=> $a->{ndv} } @have;
    ($tier, $stm, $ndv) = ($top->{tier}, $top->{statements}, $top->{ndv});
    if (keys(%u) > 1) { push @notes, 'merged: st of ' . $top->{key} . '; verify ' . join(' ', map { $_->{key} } grep { $_ != $top } @have) }
    $status = @miss ? 'partial' : 'ready'; }
  else { $status = 'none' }
  push @notes, 'no result: ' . join(' ', @miss) if @miss;
  push @notes, 'stability forms differ: ' . join('/', sort keys %sf) if keys(%sf) > 1;
  push @notes, 'gain forms differ: ' . join('/', sort keys %gf) if keys(%gf) > 1;
  # the form of each kind: that of the largest row of the kind with a result
  # (the merged st is the st of the largest row, measured in its form); a
  # vote broke ties by hash order before 11/19/2026 04:30, so two runs of
  # this script could name different forms for one example
  my ($s1, $g1);
  for my $kind (['stability', \$s1, \%sf], ['Hinf_gain', \$g1, \%gf]) {
    my ($re, $ref, $h) = @$kind; next unless keys %$h;
    my ($t) = sort { $b->{ndv} <=> $a->{ndv} || $a->{key} cmp $b->{key} } grep { $exe{$_->{key}} =~ /$re/ } @have;
    if ($t) { (my $L = $t->{form}) =~ s/d$//; $$ref = $L if $L =~ /^[PQ]$/ }
    ($$ref) = sort { $h->{$b} <=> $h->{$a} || $a cmp $b } keys %$h unless defined $$ref; }
  if (my $v = $ovr{$g}) { ($tier, $stm, $s1, $g1) = ($v->{tier}, $v->{statements}, $v->{stability_form}, $v->{Hinf_gain_form});
    $status = 'ready'; @notes = ("override: $v->{note}"); $ndv = ''; }
  $cnt{$status}++;
  (my $q = $stm) =~ s/"/""/g;
  print $o join(',', $g, $grp{$g}{file}, $grp{$g}{iv}, join(' ', @rows), join(' ', map { s/^PIETOOLS_//r } @ex),
    $s1 // '', $g1 // '', $tier, "\"$q\"", $ndv, $status, '"' . join('; ', @notes) . '"'), "\n"; }
close $o;
printf STDERR "examples %d: %s\n", scalar(keys %grp), join(', ', map { "$_ $cnt{$_}" } sort keys %cnt);
