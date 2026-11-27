#!/usr/bin/perl
# p3_cmp.pl OLD_OUT NEW_OUT - 11/27/2026: Phase 3 identity of the opvar rename. Per row of the two
# runs (d3run.m): the status in each, and the sequence of SDP records (Sedumi2Mosek inputs A, b, c,
# K; MosekSol2SedumiSol outputs x, y; sedumi inputs and outputs), compared digest by digest. A row
# is IDENTICAL when both sequences have the same length and every record has the same kind and
# digests. Prints one line per row and a summary.
use strict; use warnings;
my ($O,$N) = @ARGV;
sub res { my $f = shift; my (%s,@o); open(my $h,'<',$f) or die "$f: $!";
  while (<$h>) { s/\r?\n$//; my @x = split /\t/, $_, -1; next unless @x >= 4;
    push @o, $x[1] unless exists $s{$x[1]}; $s{$x[1]} = [$x[2],$x[3],$x[4] // '']; } return (\%s,\@o) }
sub sdp { my $f = shift; my %r; return \%r unless -f $f; open(my $h,'<',$f) or die "$f: $!";
  while (<$h>) { s/\r?\n$//; my @x = split /\t/, $_, -1; push @{$r{$x[0]}}, [$x[2],$x[3],$x[4] // '']; } return \%r }
my ($so,$oo) = res("$O/results.tsv"); my ($sn,$on) = res("$N/results.tsv");
my $do = sdp("$O/sdp.tsv"); my $dn = sdp("$N/sdp.tsv");
my %seen; my @ids = grep { !$seen{$_}++ } (@$on, @$oo);
my ($rows,$same,$stdiff,$nsdp,$recs) = (0,0,0,0,0);
printf "%-32s %-13s %-9s %s\n", 'row', 'status o/n', 'SDPs o/n', 'verdict';
for my $id (@ids) {
  my $a = $so->{$id}; my $b = $sn->{$id};
  my $st = sprintf "%s/%s", $a ? $a->[0] : '-', $b ? $b->[0] : '-';
  my @x = @{$do->{$id} // []}; my @y = @{$dn->{$id} // []};
  my $cx = grep { $_->[0] eq 'S2M' || $_->[0] eq 'SEDin' } @x;
  my $cy = grep { $_->[0] eq 'S2M' || $_->[0] eq 'SEDin' } @y;
  my $v;
  if (!$a || !$b) { $v = 'in one run only' }
  else {
    $rows++; $stdiff++ if $a->[0] ne $b->[0];
    if (@x != @y) { $v = "DIFFER: ".scalar(@x)." vs ".scalar(@y)." records" }
    else {
      my $k = 0; $k++ while $k < @x && $x[$k][0] eq $y[$k][0] && $x[$k][1] eq $y[$k][1];
      if ($k == @x) { $v = @x ? 'IDENTICAL' : 'no SDP'; $same++; $nsdp += $cx; $recs += @x }
      else { $v = "DIFFER at record ".($k+1)." ($x[$k][0]): old [$x[$k][2]] new [$y[$k][2]]" }
    }
  }
  printf "%-32s %-13s %-9s %s\n", $id, $st, "$cx/$cy", $v;
}
print "\nrows in both runs: $rows; identical SDP sequences (or none): $same; status differs: $stdiff\n";
print "SDPs identical bit for bit (A, b, c, K, x, y): $nsdp, in $recs records\n";
