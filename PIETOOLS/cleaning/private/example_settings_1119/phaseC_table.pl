#!/usr/bin/perl
# phaseC_table.pl DIR HARNESS > PHASE_C_TABLE.md
# One row per example of the plan DIR/phaseC_plan.csv: the form of its
# tests, the tier and the changes of its st, the measured result and
# decision variables (DIR/phaseC_numbers.csv: the harness verification) with
# the variables of the light tier as is, and the verdict of the library path
# (HARNESS/libcheck_results.csv for the PDE files, HARNESS/blockcheck_results.csv
# for the DDE and NDS/DDF blocks; the last line of an id and test counts).
# MMP, 11/19/2026.
use strict; use warnings;
my ($D, $H) = @ARGV;
sub csv { my ($f) = @_; my @r; open(my $h, '<', $f) or return (); my $hd = <$h>; $hd =~ s/\r?\n$//;
  my @c = split /,/, $hd;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = ($l =~ /("(?:[^"]|"")*"|[^,]*)(?:,|$)/g);
    s/^"(.*)"$/$1/s for @x; s/""/"/g for @x; my %e; @e{@c} = @x; push @r, \%e; } return @r; }
sub fmt { my $n = shift; return '' unless defined $n && $n =~ /^\d+$/; 1 while $n =~ s/^(\d+)(\d{3})/$1,$2/; return $n; }
my %num = map { $_->{key} => $_ } csv("$D/phaseC_numbers.csv");
my %case; { open(my $h, '<', "$D/case_file.txt") or die;
  while (<$h>) { s/\r?\n$//; my ($c, $f) = map { s/^\s+|\s+$//gr } split /\|/; next unless $f;
    $case{$f} = $c; } }
$case{'PIETOOLS_PDDE_Ex'} = '50';
my (%lib, %libf, %ids);                    # file:id and file -> test -> status (last line wins);
                                           # by file only where the file has one id (a file with
                                           # several cases, e.g. Saba, showed its last case on every row)
for my $e (csv("$H/libcheck_results.csv")) { next unless $e->{executive}; my $v = $e->{value} =~ /^[\d.]+$/ ? sprintf(" %.6g", $e->{value}) : ""; my $n = $e->{ndv} =~ /^[\d.e+]+$/ ? " at " . fmt(sprintf("%.0f", $e->{ndv})) : ""; $lib{"$e->{file}:$e->{id}"}{$e->{executive}} = "$e->{status}$v$n"; $libf{$e->{file}}{$e->{executive}} = "$e->{status}$v$n"; $ids{$e->{file}}{$e->{id}} = 1; }
my %blk;                                   # lib:id -> test -> verdict
for my $e (csv("$H/blockcheck_results.csv")) { $blk{"$e->{lib}:$e->{id}"}{$e->{test}} = $e->{verdict} . ($e->{settings_match} eq 'true' ? '' : " (settings $e->{settings_match})"); }
print "| Example | Case | Form | Tier and changes | Measured (harness) | Variables | Light as is | Library path |\n";
print "|---|---|---|---|---|---|---|---|\n";
for my $p (csv("$D/phaseC_plan.csv")) {
  my @rows = split / /, $p->{rows};
  my $form = join(' ', grep { $_ ne '' } ($p->{stability_form} ne '' ? "stab $p->{stability_form}" : ''), ($p->{Hinf_gain_form} ne '' ? "gain $p->{Hinf_gain_form}" : ''));
  $form = 'oco' if $p->{file} eq '?' || $p->{file} eq 'PIETOOLS_PDDE_Ex';
  my $chg = $p->{statements} eq '' ? '' : join('; ', map { s/^st\.//r =~ s/;$//r } split /(?<=;)\s+(?=st\.)/, $p->{statements});
  my $tc = $p->{status} eq 'none' ? '–' : $p->{tier} . ($chg ne '' ? ": $chg" : '');
  my (@res, @nd, @n0);
  for my $k (@rows) { my $n = $num{$k} or next; my $t = $k =~ m{/(\w+)$} ? " ($1)" : '';
    next if $n->{verified} eq '0';         # an earlier plan (PDE 50: the form 'o'), not this st
    push @res, ($n->{result} ne '' ? $n->{result} : 'certified') . $t; push @nd, fmt($n->{ndv}) . $t;
    push @n0, fmt($n->{ndv0}) . $t if $n->{ndv0} ne ''; }
  my $ex = $p->{file} eq '?' ? 'PIETOOLS_PDDE_Ex' : $p->{file};
  my $cs = $p->{group} =~ /^(DDE|NDSDDF):(\S+)/ ? $2 : ($p->{index} ne '' ? $p->{index} : ($case{$ex} // ($ex =~ /Robin|Example_Diagne|Boundary_Delay_2D/ ? 'X' : '')));
  my $lp = '';
  if ($p->{group} =~ /^(DDE|NDSDDF):(\S+)/) { my $b = $blk{"$1:$2"}; $lp = $b ? join(', ', map { "$_ $b->{$_}" } sort keys %$b) : ''; }
  else { my $l = $lib{"$ex:$cs"} // $lib{"$ex:" . ($cs =~ s/\.0$//r)} // ((keys %{$ids{$ex} // {}}) == 1 ? $libf{$ex} : undef); $lp = $l ? join(', ', map { "$_ $l->{$_}" } sort keys %$l) : ''; }
  my $name = $ex =~ s/^PIETOOLS_PDE_Ex_//r =~ s/^examples_(\w+)_library_PIETOOLS$/$1 library/r;
  print "| $name | $cs | $form | $tc | ", join('; ', @res) || ($p->{status} eq 'none' ? 'no certified setting' : ''),
        " | ", join('; ', @nd), " | ", join('; ', @n0), " | $lp |\n";
}
