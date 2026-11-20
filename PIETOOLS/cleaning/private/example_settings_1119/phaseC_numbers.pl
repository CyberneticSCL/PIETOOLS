#!/usr/bin/perl
# phaseC_numbers.pl HARNESS_DIR DIR
# The numbers each example note cites, from the verification runs of the
# planned settings (results_EXV1B.csv, results_EXV2.csv; results_EXV1.csv
# for rows not rerun): per row the verified result in words (a certificate;
# the threshold of the bisected parameter, named from lb_meta.m; gamma, the
# H2 norm or the closed-loop gain), the decision variables and rows m, and
# the size of the same row at the light tier (results_M1L.csv, B1c0.csv) for
# comparison. Writes DIR/phaseC_numbers.csv: key, result, ndv, m, ndv0,
# tier0, verified. MMP, 11/19/2026.
use strict; use warnings;
my ($H, $D) = @ARGV;
sub isnum { my $v = shift; return defined $v && $v =~ /^-?[0-9.]+(?:[eE][-+]?\d+)?$/ }
my %par; { open(my $h, '<', "$H/lb_meta.m") or die;
  # entries: 'id', name or [P,'name'], then unquoted fields (DDE: [lines], a
  # number; PDE: one flag), then 'parameter'
  while (<$h>) { next unless /^\s*'([^']+)',\s*(?:\[[^\]]*\]|'[^']*')\s*,(?:\s*[^,']*\s*,)*?\s*'([^']*)'/; $par{$1} = $2 if $2 ne ''; } close $h; }
my %ver;
for my $f ('EXV1', 'EXV1B', 'EXV2', 'EXVX') { open(my $h, '<', "$H/results_$f.csv") or next; <$h>;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = split /,/, $l, 10; next if @x < 10;
    $ver{"$x[0]:$x[1]"} = { st => $x[3], res => $x[4], ndv => $x[8], note => $x[9], src => $f }; } close $h; }
my %light;
for my $f ('M1L', 'B1c0') { open(my $h, '<', "$H/results_$f.csv") or next; <$h>;
  while (my $l = <$h>) { my @x = split /,/, $l, 10; $light{"$x[0]:$x[1]"} = $x[8] if isnum($x[8]); } close $h; }
open(my $o, '>', "$D/phaseC_numbers.csv") or die;
print $o "key,result,ndv,m,ndv0,tier0,verified\n";
for my $k (sort keys %ver) { my $v = $ver{$k}; my ($lib, $id) = split /:/, $k, 2; (my $base = $id) =~ s{/.*}{};
  my ($m) = $v->{note} =~ /\bm=(\d+)/; my ($ex) = $v->{note} =~ /^sop:PIETOOLS_(\w+)/; $ex //= '';
  my $txt;
  if ($v->{st} =~ /^(converged|max_solves|hi_pass)$/ && isnum($v->{res})) {
    my $p = $par{$base} // q(the parameter); $p =~ s/^(\w+)=(.+)$/$1 ($2)/; $txt = "certified up to $p = " . sprintf('%.6g', $v->{res}); }
  elsif ($v->{st} eq 'certified' && $ex =~ /Hinf_gain/ && isnum($v->{res})) { $txt = 'gamma ' . sprintf('%.6g', $v->{res}) }
  elsif ($v->{st} eq 'certified' && $ex =~ /H2_norm/ && isnum($v->{res})) { $txt = 'H2 norm ' . sprintf('%.6g', $v->{res}) }
  elsif ($v->{st} eq 'certified' && $ex =~ /control|estimator/ && isnum($v->{res})) { $txt = 'closed-loop gain ' . sprintf('%.6g', $v->{res}) }
  elsif ($v->{st} eq 'certified') { $txt = 'certified' }
  else { $txt = '' }
  print $o join(',', $k, "\"$txt\"", $v->{ndv}, $m // '', $light{$k} // '', 'light', ($txt ne '' ? 1 : 0)), "\n"; }
close $o;
