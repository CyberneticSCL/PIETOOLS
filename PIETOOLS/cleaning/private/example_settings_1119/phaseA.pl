#!/usr/bin/perl
# phaseA.pl HARNESS_DIR OUT_PREFIX
# Phase A of the per-example settings (maintainer, 11/19/2026: "all examples
# should have optimized settings files. There is no default."): from the
# results files of the 11/18 studies (current tree; the 11/17 pass only where
# no study run covers a row and tier), for each library row
#   1. the best certified result over every measured (form, setting):
#      bisection rows the largest certified parameter, value rows (gains,
#      norms) the smallest certified value, stability rows the verdict;
#   2. the form and setting together: the cheapest candidate reaching the
#      best result in either form, 'P' on an exact tie (maintainer,
#      11/19/2026: "P form should be used when it is cheaper and just as
#      good for accuracy. Not blindly!"; supersedes 'P' wherever it
#      certifies);
#   3. cheapest: decision variables, then m, then
#      solves, then seconds.
# Certified: single solves need both gates (harness status 'certified' and
# the executive status exst= optimal/inaccurate when the note has one);
# bisections are decided by the harness gate at each step. Reaching the best:
# within RTOL = 1e-3 relative (the bisection resolution).
# Writes OUT_PREFIX.csv (one line per row) and OUT_PREFIX.md. MMP, 11/19/2026.
use strict; use warnings;
my ($H, $OUT) = @ARGV;  my $RTOL = 1e-3;
# 'P' within PBAND of the cheapest decision-variable count counts as a tie
# ('P' is the stronger form; the correction of 11/19 rejects large costs)
my $PBAND = 1.10;
# tag => [setting label, st description]; tier comes from the settings column
my %TAG = (
  M1L    => ['light', ''],
  B2dD0  => ['light, dD 0', 'st.lpi.ineq.dD = 0'],
  H1ctl  => ['heavy', ''],
  B3q522 => ['heavy, deg_Q{1} [5 2 2]', 'st.sop.deg_Q{1} = [5 2 2]'],
  B3q622 => ['heavy, deg_Q{1} [6 2 2]', 'st.sop.deg_Q{1} = [6 2 2]'],
  B4ctl  => ['veryheavy', ''],
  B4q866 => ['veryheavy, deg_Q{1} [8 6 6]', 'st.sop.deg_Q{1} = [8 6 6]'],
  C1ext  => ['extreme', ''],
  C2ext  => ['extreme, stripped 1-D entries', 'stripped1d'],
  C3ext  => ['extreme, stripped 1-D entries, dD 0', 'stripped1d; st.lpi.ineq.dD = 0'],
  C1ext2 => ['extreme', ''],
  C2ext2 => ['extreme, stripped 1-D entries', 'stripped1d'],
  # Phase B (11/19): stripped and light 'P' or as is, with dD 0 ('d0')
  XSP    => ['stripped', ''],  XSPd0 => ['stripped, dD 0', 'st.lpi.ineq.dD = 0'],
  XLPd0  => ['light, dD 0', 'st.lpi.ineq.dD = 0'],
  XSQ    => ['stripped', ''],  XSQd0 => ['stripped, dD 0', 'st.lpi.ineq.dD = 0'],
  X2SP   => ['stripped', ''],  X2SPd0 => ['stripped, dD 0', 'st.lpi.ineq.dD = 0'],
  X2LPd0 => ['light, dD 0', 'st.lpi.ineq.dD = 0'],
  X2SQ   => ['stripped', ''],  X2SQd0 => ['stripped, dD 0', 'st.lpi.ineq.dD = 0'],
  X2HPd0 => ['heavy, dD 0', 'st.lpi.ineq.dD = 0'],
  XVP8   => ['veryheavy', ''], XVH => ['veryheavy', ''],
  D1P    => ['light', ''],  D1Ph => ['heavy', ''],
  D1Ppde => ['light', ''],  D1Ppdeh => ['heavy', ''],  D2Pl => ['light', ''],
  D2st   => ['light, DDE block st (eec4f07d)', 'blockst'],
  GallQ  => ['light', ''],  GallP => ['light', ''],  GallQh => ['heavy', ''],  GallPh => ['heavy', ''],
  B1c0   => ['light', ''],
  B1c1   => ['light, 2-D dD 0', 'st.lpi.ineq.dD = 0'],
  B1c2   => ['light, 2-D dD 0, deg_Q{2} 0 raised by degree (T2)', "st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0; st.sop.raise_Q = 'deg'; st.sop.max_raise_Q(2) = 2"],
  B1c3   => ['light, 2-D deg_Q{2} 0 raised by degree (T3)', "st.sop.deg_Q{2} = 0; st.sop.raise_Q = 'deg'; st.sop.max_raise_Q(2) = 2"],
  H2ctl  => ['heavy', ''],
  H2T3   => ['heavy, 2-D T3', "st.lpi.ineq.dD = 1; st.sop.deg_Q{2} = 0; st.sop.raise_Q = 'deg'; st.sop.max_raise_Q(2) = 2"],
  E3vh   => ['veryheavy', ''],  E1str => ['stripped', ''],  E2l => ['light', ''],  E2h => ['heavy', ''],
);
my %rows; my %reach;      # rows{key} = [candidates]; reach{key} = [T-bisections]
sub note_f { my ($n, $re) = @_; my @m = ($n =~ /$re/g); return @m ? $m[-1] : '' }
sub add { my ($key, $c) = @_; push @{$rows{$key}}, $c }
sub rd { my ($file, $tag, $lab, $desc, $only_tiers) = @_;
  open(my $h, '<', $file) or return; <$h>;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = split /,/, $l, 10; next if @x < 10;
    my ($lib, $id, $tier, $st, $res, $doc, $n, $sec, $ndv, $note) = @x;
    next if $only_tiers && !$only_tiers->{$tier};
    next if $st eq 'skipped';
    my $c = { tag=>$tag, tier=>$tier, lab=>$lab, desc=>$desc, st=>$st, res=>$res, doc=>$doc, n=>$n, sec=>$sec,
      ndv=>$ndv, form=>note_f($note, qr{_sop/(\w+)}), m=>note_f($note, qr{\bm=(\d+)}),
      exst=>note_f($note, qr{exst=(\w+)}), eta=>note_f($note, qr{eta=([-0-9.eE+]+)}),
      fr=>note_f($note, qr{\bfr=([-0-9.eE+]+)}), exec=>note_f($note, qr{^sop:(\w+)}),
      passes=>note_f($note, qr{passes=(\d+)}) };    # executive raise passes of a single solve
    # light rows seeded into H2ctl/H2T3 from B1c0 are B1c0 measurements
    next if ($tag eq 'H2ctl' || $tag eq 'H2T3') && $tier eq 'light';
    if ($id =~ m{^(\d+)/reach$}) { push @{$reach{"$lib:$1"}}, $c; next }
    add("$lib:$id", $c); }
  close $h; }
for my $t (sort keys %TAG) { my $lab = $TAG{$t}[0];
  # tags whose label omits the tier take it from the settings column
  rd("$H/results_$t.csv", $t, $lab, $TAG{$t}[1]); }
# the B1 tightness bisections (T<k>_<system>): reach of light, T1, T2, T3
for my $f (glob("$H/results_T[0-3]_*.csv")) { my ($k) = $f =~ /results_T(\d)_/;
  my @lab = ('light', 'light, 2-D dD 0 (T1)', 'light, 2-D dD 0, deg_Q{2} 0 raised (T2)', 'light, 2-D deg_Q{2} 0 raised (T3)');
  rd($f, "T$k", $lab[$k], ''); }
# the example-file settings of 11/18 (EXAMPLE_SETTINGS.md Sections 3.3, 6, 7,
# 11.3; one solve at the library defaults, form 'Q', scratch exst_* scripts):
# [row, label, st, variables, m, harness status, executive status]
for my $e (
  ['PDE:33', 'light, dD 0, deg_Q{2} 0', 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0', 409, 111, 'certified', 'optimal'],
  ['PDE:35', 'light, dD 0, deg_Q{2} 0', 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0', 8009, 815, 'certified', 'optimal'],
  ['PDE:37', 'light, dD 0, deg_Q{2} [1 0 0]', 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [1 0 0]', 16656, 1117, 'certified', 'optimal'],
  ['PDE:40', 'light, dD 0, deg_Q{2} 0', 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0', 8518, 844, 'certified', 'optimal'],
  ['PDE:38', "light, deg_Q{2} 2, psatz{2} 'none', max_raise_Q(2) 0", "st.sop.deg_Q{2} = 2; st.lpi.ineq.psatz{2} = 'none'; st.sop.max_raise_Q(2) = 0", 595684, 30999, 'certified', 'optimal'],
  ['PDE:43', 'light, dD 0, deg_Q{2} 2', 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 2', 196513, 6376, 'certified', 'optimal'],
  ['PDE:57', 'light, dD 0, deg_Q{2} [2 1 1]', 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [2 1 1]', 23121, 2344, 'certified', 'optimal'],
  ['PDE:58', 'light, dD 0, deg_Q{2} 2', 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 2', 225705, 6608, 'certified', 'optimal'],
  ['PDE:66', 'light, dD 0', 'st.lpi.ineq.dD = 0', 16740, 1301, 'certified', 'optimal'],
  ['PDE:72', 'light, dD 0, deg_Q{2} 2', 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 2', 241648, 6257, 'unknown', 'optimal'],
  ['PDE:61', 'light, dD 0', 'st.lpi.ineq.dD = 0', 197415, 6969, 'certified', 'optimal'] ) {
  my ($k, $lab, $desc, $ndv, $m, $hs, $es) = @$e;
  add($k, { tag=>'EXST', tier=>'light', lab=>"$lab (example file)", desc=>$desc, st=>$hs, res=>1, doc=>'', n=>1, sec=>'',
    ndv=>$ndv, form=>'Q', m=>$m, exst=>$es, eta=>'', fr=>'', exec=>'PIETOOLS_PIE2PDEstability', passes=>1 }); }
# the 11/17 pass, only for (row, tier) pairs no study run measured
my %have; for my $k (keys %rows) { for my $c (@{$rows{$k}}) { $have{"$k|$c->{tier}"} = 1 if $c->{lab} eq $c->{tier} } }
{ open(my $h, '<', "$H/results_all_1118b.csv") or die; <$h>;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = split /,/, $l, 10; next if @x < 10;
    my ($lib, $id, $tier) = @x[0..2]; next if $have{"$lib:$id|$tier"};
    my $note = $x[9];
    my $c = { tag=>'P1117', tier=>$tier, lab=>"$tier (11/17 tree)", desc=>'', st=>$x[3], res=>$x[4], doc=>$x[5], n=>$x[6],
      sec=>$x[7], ndv=>$x[8], form=>note_f($note, qr{_sop/(\w+)}), m=>note_f($note, qr{\bm=(\d+)}),
      exst=>note_f($note, qr{exst=(\w+)}), eta=>note_f($note, qr{eta=([-0-9.eE+]+)}),
      fr=>note_f($note, qr{\bfr=([-0-9.eE+]+)}), exec=>note_f($note, qr{^sop:(\w+)}),
      passes=>note_f($note, qr{passes=(\d+)}) };
    next if $x[3] eq 'skipped'; add("$lib:$id", $c); }
  close $h; }
sub isnum { my $v = shift; return defined $v && $v =~ /^-?[0-9.]+(?:[eE][-+]?\d+)?$/ }
sub bis { my $s = shift; return $s =~ /^(converged|max_solves|lo_fail|hi_pass)$/ ? 1 : 0 }
sub kind { my ($cs) = @_;
  for my $c (@$cs) { return 'bisect' if bis($c->{st}) }
  for my $c (@$cs) { return 'value' if $c->{exec} =~ /Hinf_|H2_/ }
  return 'stab'; }
sub cert { my ($c, $kind) = @_;            # returns (certified, gates-disagree)
  if ($kind eq 'bisect') { return ($c->{st} =~ /^(converged|max_solves|hi_pass)$/ && isnum($c->{res})) ? (1,0) : (0,0) }
  my $hc = $c->{st} eq 'certified' ? 1 : 0;
  my $ec = $c->{exst} eq '' ? $hc : ($c->{exst} =~ /^(optimal|inaccurate)$/ ? 1 : 0);
  return ($hc && $ec && ($kind ne 'value' || isnum($c->{res}))) ? (1,0) : (0, ($hc != $ec) ? 1 : 0); }
sub isP { my $f = shift; return $f =~ /^Pd?$/ ? 1 : 0 }
sub cheaper { my ($a, $b) = @_;
  for my $f ('ndv', 'm', 'passes', 'n', 'sec') { my ($x, $y) = ($a->{$f}, $b->{$f});
    next unless isnum($x) && isnum($y);               # a measure missing on either side ties
    return $x <=> $y if $x != $y }
  my $t = (isP($b->{form}) <=> isP($a->{form}));             # exact tie: 'P', the stronger form
  return $t if $t;
  return ($b->{tag} eq 'EXST') <=> ($a->{tag} eq 'EXST'); }   # then keep the example-file st
open(my $csv, '>', "$OUT.csv") or die; open(my $md, '>', "$OUT.md") or die;
# the starts of the per-row knob search (ex_search.m): tier, changes, form, target
open(my $sts, '>', "${OUT}_starts.csv") or die;
print $sts "key,kind,form,tier,mods,best\n";
print $csv "library,id,kind,form,choice,st,result,best,best_by,variables,m,solves,seconds,candidates,flags\n";
my (%cnt, @md);
for my $k (sort { my ($la,$ia) = split /:/, $a, 2; my ($lb,$ib) = split /:/, $b, 2;
                  $la cmp $lb or (($ia =~ /^(\d+)/)[0] // 0) <=> (($ib =~ /^(\d+)/)[0] // 0) or $ia cmp $ib } keys %rows) {
  my $cs = $rows{$k}; my $kind = kind($cs); my ($lib, $id) = split /:/, $k, 2;
  my (@ok, $split);
  for my $c (@$cs) { my ($o, $s) = cert($c, $kind); push @ok, $c if $o; $split++ if $s; }
  my @flags; push @flags, "gates disagree on $split" if $split;
  my $ncand = scalar(@$cs);
  if (!@ok) { $cnt{none}++; push @flags, 'nothing certified';
    print $csv join(',', $lib, $id, $kind, '', '', '', '', '', '', '', '', '', '', $ncand, join('; ', @flags)), "\n";
    push @md, "| $k | $kind | – | nothing certified ($ncand runs) | | | | | | ".join('; ', @flags)." |"; next; }
  my $best; my $bestc;
  for my $c (@ok) { my $v = $kind eq 'stab' ? 1 : $c->{res};
    if (!defined $best || ($kind eq 'bisect' && $v > $best) || ($kind eq 'value' && $v < $best)) { $best = $v; $bestc = $c } }
  my $reaches = sub { my $c = shift; return 1 if $kind eq 'stab';
    return $kind eq 'bisect' ? ($c->{res} >= $best - $RTOL*abs($best)) : ($c->{res} <= $best + $RTOL*abs($best)); };
  my @r = grep { $reaches->($_) } @ok;
  my @rP = grep { isP($_->{form}) } @r;
  my @okP = grep { isP($_->{form}) } @ok;
  if (@okP && !@rP) { my $pb; for my $c (@okP) { $pb = $c->{res} if !defined $pb || ($kind eq 'bisect' ? $c->{res} > $pb : $c->{res} < $pb) }
    push @flags, sprintf("'P' certifies %s against %s", $kind eq 'stab' ? 'yes' : $pb, $best); $cnt{Pweaker}++; }
  my ($ch) = sort { cheaper($a, $b) } @r;
  if (@rP && !isP($ch->{form})) { my ($cp) = sort { cheaper($a, $b) } @rP;
    if (isnum($cp->{ndv}) && isnum($ch->{ndv}) && $cp->{ndv} <= $PBAND*$ch->{ndv}) {
      push @flags, sprintf("'P' within %d%% of the cheapest (%s against %s)", 100*($PBAND-1)+0.5, $cp->{ndv}, $ch->{ndv});
      $ch = $cp; $cnt{Pband}++; }
    else { push @flags, sprintf("'P' reaches it at %s variables", $cp->{ndv}); $cnt{Pcostlier}++; } }
  $cnt{$ch->{form} =~ /^P/ ? 'P' : 'other'}++; $cnt{tag}{$ch->{tag}}++;
  push @flags, '11/17 tree only' if $ch->{tag} eq 'P1117';
  my $res = $kind eq 'stab' ? 'certified' : $ch->{res};
  print $sts join(',', $k, $kind, $ch->{form}, $ch->{tier}, "\"$ch->{desc}\"", ($kind eq 'stab' ? '' : $best)), "\n";
  print $csv join(',', $lib, $id, $kind, $ch->{form}, "\"$ch->{lab}\"", "\"$ch->{desc}\"", $res, $best, $bestc->{tag},
    $ch->{ndv}, $ch->{m}, $ch->{n}, $ch->{sec}, $ncand, '"'.join('; ', @flags).'"'), "\n";
  push @md, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |", $k, $kind, $ch->{form} || '–', $ch->{lab},
    $res, ($kind eq 'stab' ? '' : "$best ($bestc->{tag})"), $ch->{ndv}, $ch->{m}, $ncand, join('; ', @flags));
}
print $md "| row | kind | form | cheapest setting reaching the best | result | best (source) | variables | m | runs | flags |\n";
print $md "|---|---|---|---|---|---|---|---|---|---|\n";
print $md "$_\n" for @md;
close $csv; close $md;
printf STDERR "rows %d: form P %d, other %d, nothing certified %d, 'P' weaker than best %d\n",
  scalar(keys %rows), $cnt{P}//0, $cnt{other}//0, $cnt{none}//0, $cnt{Pweaker}//0;
printf STDERR "chosen from: %s\n", join(' ', map { "$_:$cnt{tag}{$_}" } sort keys %{$cnt{tag}});
