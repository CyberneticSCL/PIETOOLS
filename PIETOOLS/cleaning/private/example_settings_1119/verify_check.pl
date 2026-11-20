#!/usr/bin/perl
# verify_check.pl RESULTS DIR [OUT]
# Checks a verification run (results_EXV<k>.csv of st_run EXV1/EXV2: every
# row of an example under the planned st and forms of the example) against
# the target of the row (DIR/phaseB_choices.csv): stability rows certified;
# bisection rows certified up to the best measured threshold less 1e-3 of it;
# value rows certified with a value within 1e-3 of the best. Both gates for
# single solves (harness status and executive status exst=). Prints the rows
# that fail and the counts; OUT (optional) gets one line per row: key,
# kind, status, result, best, ok, ndv, exst. MMP, 11/19/2026.
use strict; use warnings;
my ($F, $D, $OUT) = @ARGV;  my $RT = 1e-3;
sub isnum { my $v = shift; return defined $v && $v =~ /^-?[0-9.]+(?:[eE][-+]?\d+)?$/ }
my %ch; { open(my $h, '<', "$D/phaseB_choices.csv") or die; <$h>;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = ($l =~ /("(?:[^"]|"")*"|[^,]*)(?:,|$)/g); s/^"(.*)"$/$1/ for @x;
    $ch{"$x[0]:$x[1]"} = { kind => $x[2], best => $x[7] }; } close $h; }
my ($n, $ok, @bad) = (0, 0);
open(my $o, '>', $OUT // '/dev/null') or die;
print $o "key,kind,status,result,best,ok,ndv,exst\n";
open(my $h, '<', $F) or die; <$h>;
while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = split /,/, $l, 10; next if @x < 10;
  my $k = "$x[0]:$x[1]"; my $c = $ch{$k} or next; $n++;
  my ($st, $res, $ndv, $note) = @x[3, 4, 8, 9]; my ($ex) = $note =~ /exst=(\w+)/; $ex //= '';
  my $good;
  if ($c->{kind} eq 'bisect') { $good = $st =~ /^(converged|max_solves|hi_pass)$/ && isnum($res) && isnum($c->{best}) && $res >= $c->{best} - $RT * abs($c->{best}) }
  elsif ($c->{kind} eq 'value') { $good = $st eq 'certified' && ($ex eq '' || $ex =~ /optimal|inaccurate/) && isnum($res) && isnum($c->{best}) && $res <= $c->{best} + $RT * abs($c->{best}) }
  else { $good = $st eq 'certified' && ($ex eq '' || $ex =~ /optimal|inaccurate/) }
  $ok++ if $good;  push @bad, sprintf("%-28s %-6s %-10s result %-12s best %-12s ndv %-8s exst %s", $k, $c->{kind}, $st, $res, $c->{best}, $ndv, $ex) unless $good;
  print $o join(',', $k, $c->{kind}, $st, $res, $c->{best}, $good ? 1 : 0, $ndv, $ex), "\n"; }
close $h; close $o;
print "rows $n, verified $ok, not verified ", scalar(@bad), "\n", map { "  $_\n" } @bad;
