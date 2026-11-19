#!/usr/bin/perl
# cmp_runs.pl CONTROL.csv CANDIDATE.csv [LABEL]
# Compares two lb_run results files of the 11/18 studies row by row (key:
# library, id) and prints a markdown summary: verdicts by both gates (the
# harness status and the executive status 'exst=' of the note), bisection
# thresholds, decision variables, solves and seconds. A change of verdict
# counts only when both gates agree (maintainer, 11/18/2026); rows on which
# the gates disagree are listed apart. Bisection rows (status converged,
# lo_fail, hi_pass) are compared by their threshold, which the harness
# bisection decides with the harness gate. MMP, 11/18/2026.
use strict; use warnings;
my ($fc, $fn, $label) = @ARGV;  $label //= "candidate";
sub rd { my ($f) = @_; my %r; open(my $h, '<', $f) or die "$f: $!";
  my $hd = <$h>;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = split /,/, $l, 10;
    next if @x < 10; $r{"$x[0]:$x[1]"} = { lib=>$x[0], id=>$x[1], tier=>$x[2], st=>$x[3], res=>$x[4],
      doc=>$x[5], n=>$x[6], sec=>$x[7], ndv=>$x[8], note=>$x[9] }; }
  close $h; return \%r; }
my $C = rd($fc); my $N = rd($fn);
sub hcert { my $s = shift; return ($s eq 'certified') ? 1 : 0 }
sub ecert { my $n = shift; my @m = ($n =~ /exst=(\w+)/g); return undef unless @m;
  return ($m[-1] eq 'optimal' || $m[-1] eq 'inaccurate') ? 1 : 0 }
sub eta { my $n = shift; my @m = ($n =~ /eta=([-0-9.eE+naNIf]+)/g); return @m ? $m[-1] : '' }
sub bis { my $s = shift; return ($s =~ /^(converged|lo_fail|hi_pass)$/) ? 1 : 0 }
my (@gain, @loss, @split, @thr, @other); my ($sc, $sn, $vc, $vn, $nc, $nn, $both) = (0,0,0,0,0,0,0);
for my $k (sort keys %$N) {
  next unless exists $C->{$k}; my ($c, $n) = ($C->{$k}, $N->{$k}); $both++;
  $sc += $c->{sec} =~ /^[0-9.]+$/ ? $c->{sec} : 0;  $sn += $n->{sec} =~ /^[0-9.]+$/ ? $n->{sec} : 0;
  $vc += $c->{ndv} =~ /^[0-9]+$/ ? $c->{ndv} : 0;    $vn += $n->{ndv} =~ /^[0-9]+$/ ? $n->{ndv} : 0;
  $nc += $c->{n} =~ /^[0-9]+$/ ? $c->{n} : 0;        $nn += $n->{n} =~ /^[0-9]+$/ ? $n->{n} : 0;
  my $row = sprintf("| %s | %s / %s | %s / %s | %s / %s | %s / %s | %s / %s |", $k,
     $c->{st}, $n->{st}, $c->{res}, $n->{res}, $c->{ndv}, $n->{ndv}, $c->{n}, $n->{n}, $c->{sec}, $n->{sec});
  if (bis($c->{st}) || bis($n->{st})) {
    if ($c->{st} ne $n->{st} || abs(($n->{res}||0) - ($c->{res}||0)) > 1e-12) {
      my $rel = ($c->{res} =~ /^-?[0-9.eE+-]+$/ && $c->{res} != 0) ? sprintf("%+.3f%%", 100*($n->{res}-$c->{res})/abs($c->{res})) : '';
      push @thr, "$row $rel |"; }
    next; }
  my ($hc, $hn) = (hcert($c->{st}), hcert($n->{st}));
  my ($ec, $en) = (ecert($c->{note}), ecert($n->{note}));
  $ec //= $hc; $en //= $hn;
  my $cc = ($hc && $ec) ? 1 : (!$hc && !$ec) ? 0 : -1;   # -1: gates disagree
  my $cn = ($hn && $en) ? 1 : (!$hn && !$en) ? 0 : -1;
  my $g = sprintf(" H %d/%d E %d/%d eta %s/%s |", $hc, $hn, $ec, $en, eta($c->{note}), eta($n->{note}));
  if ($cc == 0 && $cn == 1) { push @gain, $row.$g }
  elsif ($cc == 1 && $cn == 0) { push @loss, $row.$g }
  elsif ($cc == -1 || $cn == -1) { push @split, $row.$g if ($hc != $hn || $ec != $en) }
  elsif ($c->{st} ne $n->{st} || $c->{res} ne $n->{res}) { push @other, $row.$g }
}
my $hdr = "| row | status ctl / $label | result | variables | solves | seconds |";
my $sep = "|---|---|---|---|---|---|";
printf "## %s against control\n\nRows in both files: %d. Summed over them: variables %d -> %d (%.2fx), solves %d -> %d, seconds %.0f -> %.0f.\n\n",
  $label, $both, $vc, $vn, ($vc ? $vn/$vc : 0), $nc, $nn, $sc, $sn;
for my $s (['Gained (both gates)', \@gain, 1], ['Lost (both gates)', \@loss, 1], ['Gates disagree', \@split, 1],
           ['Bisection thresholds that differ', \@thr, 2], ['Other differences', \@other, 1]) {
  my ($t, $a, $kind) = @$s;  printf "### %s: %d\n\n", $t, scalar(@$a);  next unless @$a;
  if ($kind == 2) { print "$hdr change |\n|---|---|---|---|---|---|---|\n" } else { print "$hdr gates |\n|---|---|---|---|---|---|---|\n" }
  print "$_\n" for @$a;  print "\n"; }
