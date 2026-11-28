#!/usr/bin/perl
# SOP_APPLY.PL - Phase 1 of the removal of the _sop suffix (../PLAN.md), 11/28/2026; adapted from
# ../../rename_copvar_opvar_1125/phase1/ren_apply.pl (no article rule).
#   perl sop_apply.pl ROOT TABLE paths
#     prints "old<TAB>new" for each file outside the records whose path holds a mapped token.
#   perl sop_apply.pl ROOT TABLE check|apply DATE EXISTING REPORTDIR
#     maps every live .m file at its path after the moves; check writes the reports only, apply
#     also writes the files. DATE: mm/dd/yyyy of the D1 log line. EXISTING: the files of the 2025
#     release, by their path before the moves (D1); each that changes gets the D1 line.
# ROOT is the PIETOOLS folder; TABLE is ../rename_table.tsv. Records, never edited: cleaning/,
# orphaned_file/, SOSTOOLS400/; in a live file, the header-log lines and the lines holding "(was)".
# On every other line of a live .m file:
#   1. the text PIETOOLS_*_sop becomes PIETOOLS_*; in code, a name built as [... '_sop'] loses
#      the suffix;
#   2. each word token in %MAP (the renamed rows of TABLE, the folders, the labels of the
#      executives) becomes its new name, in lower, upper and capitalized form; any other token
#      holding sopvar names the class and stays; any other token ending in _sop or holding _sop_
#      must be in @KEEP, else it is UNKNOWN and stops apply;
#   3. a trailing comment that followed 2 or more blanks keeps its column, so that the change
#      markers stay at column 77.
use strict; use warnings; use File::Find; use File::Path qw(make_path);

my ($root,$table,$mode,$date,$exfile,$rep) = @ARGV;
die "usage: ROOT TABLE paths | ROOT TABLE check|apply DATE EXISTING REPORTDIR\n"
  unless $root && $table && $mode && ($mode eq 'paths' || ($mode =~ /^(check|apply)$/ && $date && $exfile && $rep));

# the renamed rows of the table (new name differs), the folders among them
my @PAIRS;
open(my $T,'<',$table) or die "$table: $!"; <$T>;
while (<$T>) { s/\r?\n$//; my @f = split /\t/; next unless @f >= 4 && $f[3] ne $f[0]; push @PAIRS, [$f[0],$f[3]] }
close $T;
# the labels the executives print (label = 'H2_control_sop', ['stability_sop ' form], ...), and two
# abbreviations in a comment of bl_sop.m
push @PAIRS, map { [$_.'_sop', $_] } qw(stability Hinf_gain H2_norm well_posedness H2_control
  H2_estimator Hinf_OFcontrol Hinf_control Hinf_estimator);
push @PAIRS, ['_Hinf_gain_sop','_Hinf_gain'], ['_H2_norm_sop','_H2_norm'];
# kept: names of this table that stay, variables, file flags, identifiers named after moved
# scripts, and retired names (in comments)
my @KEEP = qw(_sop exec_sop_settings lb_full_sop lb_sop_plan bl_sop rand_sop PAUSE_sop results_sop
  is_sop t_sop s_sop peak_sop_MB max_rel_opvar_vs_sop ntest_settings_prog_sop ntest_lpi_ineq_sop
  stability_1d_sop heat3d_stability_sop stability_2d_sop hinf_gain_1d_sop
  lpiprogram_sop get_lift_degs_sop getController_direct_sop op2copvar_sop eye_copvar_sop
  zeros_copvar_sop pi4_copvar_sop mat2copvar_sop poslpivar_sop PIETOOLS_stability_2D_sop
  PIETOOLS_H2_norm_c_sop);
my %KEEP = map { lc($_) => 1 } @KEEP;
my %MAP;
for my $p (@PAIRS) { my ($o,$n) = @$p;
  $MAP{$o} = $n; $MAP{uc $o} = uc $n; $MAP{ucfirst $o} = ucfirst $n }
my $CAND = qr/_sop(?:var)?(?:$|_)/i;     # a token that may name a renamed function

my (%stat,%unknown,@realign,@records,@d1,@changed,@nolog,@star);
sub maptok { my ($t,$where) = @_;
  return $t unless $t =~ $CAND;
  if (exists $MAP{$t}) { $stat{map}{"$t\t$MAP{$t}"}++; return $MAP{$t} }
  if ($t =~ /sopvar/i)   { $stat{class}{$t}++; return $t }   # the class (sv_test_sopvar_kron, ...)
  if ($KEEP{lc $t})      { $stat{keep}{$t}++;  return $t }
  $unknown{"$t\t$where"}++; return $t }
sub mappath { my $p = shift; join '/', map { (my $c=$_) =~ s/(\w+)/maptok($1,"path $p")/ge; $c } split m#/#, $p }

# code and comment of a line; the comment starts at % or ... outside strings
sub split_line { my $l=shift; my ($o,$i,$n,$p)=('',0,length $l,'');
  while ($i<$n) { my $c=substr($l,$i,1);
    return ($o, substr($l,$i)) if $c eq '%' || substr($l,$i,3) eq '...';
    if ($c eq '"' || ($c eq "'" && $p !~ /[\w\)\]\}\.']/)) { my $q=$c; my $j=$i+1;
      while ($j<$n) { my $d=substr($l,$j,1); if ($d eq $q) { if ($j+1<$n && substr($l,$j+1,1) eq $q) { $j+=2; next } last } $j++ }
      $o.=substr($l,$i,$j-$i+1); $i=$j+1; $p=$q; next }
    $o.=$c; $p = ($c eq ' ' || $c eq "\t") ? ' ' : $c; $i++ }
  return ($o, '') }

# header-log lines (the rule of ren_apply.pl)
my $ENTRY = qr/^%\s*(?:Initial coding|[A-Z]{2,3}(?:,? ?[A-Z]{2,3})*(?: \(CC\))?,? ?[ -] ?\d{1,2}\/\d{1,2}\/\d{2,4}|[A-Z]{2,3}(?:, [A-Z]{2,3})*, \d{1,2}\/\d{1,2}\/\d{2,4}|\d{1,2}\/\d{1,2}\/\d{2,4}\s*-)/;
sub header_log { my $L=shift; my ($first,$seen)=(0,0);
  for my $i (0..$#$L) { my $l=$L->[$i];
    next if $l =~ /^\s*(?:function|classdef)\b/;
    next if !$seen && $i && $L->[$i-1] =~ /\.\.\.\s*$/;
    if ($l =~ /^\s*%/) { $seen=1; next }
    next unless $l =~ /\S/;
    next if !$seen && $i < 3;
    $first=$i+1; last }
  my ($in,%log)=(0);
  for my $i (0..($first ? $first-2 : -1)) { my $l=$L->[$i];
    if ($l =~ $ENTRY) { $in=1; $log{$i+1}=1; next }
    if ($in && $l =~ /^%\s{2,}\S/) { $log{$i+1}=1; next }
    $in=0 }
  return \%log }

sub mapseg { my ($s,$where,$iscode) = @_;
  if ($s =~ s/\bPIETOOLS_\*_sop\b/PIETOOLS_*/g) { push @star, $where }
  # a name built at run time, ['PIETOOLS_' nm{1} '_sop'] (test_executives_sop.m): drop the suffix
  if ($iscode && $s =~ s/\s+'_sop'\]/]/g) { push @star, "$where\tbuilt" }
  $s =~ s/(\w+)/maptok($1,$where)/ge; return $s }

my @files;
find({ no_chdir=>1, wanted=>sub { my $f=$File::Find::name; return unless -f $f;
  (my $r=$f) =~ s#^\Q$root\E/##; return if $r =~ m#^(cleaning|orphaned_file|SOSTOOLS400)/#;
  push @files, $r }}, $root);
@files = sort @files;

if ($mode eq 'paths') {
  for my $r (@files) { my $n = mappath($r); print "$r\t$n\n" if $n ne $r }
  for my $k (sort keys %unknown) { print STDERR "UNKNOWN\t$k\n" }
  exit(%unknown ? 1 : 0) }

open(my $E,'<',$exfile) or die "$exfile: $!";
my %existing = map { s/\r?\n$//r => 1 } grep { /\S/ } <$E>; close $E;
my %existing_new = map { mappath($_) => $_ } keys %existing; %unknown = (); %stat = ();
my ($nfiles,$nlines) = (0,0);
for my $r (grep { /\.m$/ } @files) {
  open(my $h,'<:raw',"$root/$r") or die "$r: $!"; local $/; my $s=<$h>; close $h;
  my @L; while ($s =~ /\G([^\r\n]*)(\r\n|\n|\z)/gc) { push @L,[$1,$2]; last if $2 eq '' }
  pop @L if @L && $L[-1][0] eq '' && $L[-1][1] eq '';
  my @T = map { $_->[0] } @L;
  next unless grep { /_sop/i } @T;
  my $log = header_log(\@T);
  my %rec; for my $i (0..$#T) {
    if ($log->{$i+1}) { $rec{$i}='log' } elsif ($T[$i] =~ /\(was\)/) { $rec{$i}='was' } }
  my $nch = 0;
  for my $i (0..$#T) { my $t=$T[$i];
    if ($rec{$i}) { push @records, "$r:".($i+1)."\t$rec{$i}\t$t" if $t =~ /_sop/i; next }
    my ($code,$cmt) = split_line($t); my $where = "$r:".($i+1);
    my $nc = mapseg($code,$where,1); my $nm = mapseg($cmt,$where,0);
    my ($gap) = $code =~ /( *)$/;
    if ($cmt =~ /^%/ && $code =~ /\S/ && length($gap) >= 2 && $nc ne $code) {
      (my $core=$nc) =~ s/ +$//; my $ng = length($code) - length($core); $ng = 1 if $ng < 1;
      $nc = $core . (' ' x $ng); push @realign, "$where\t".length($code)."\t".(length($core)+$ng) }
    my $new = $nc.$nm;
    if ($new ne $L[$i][0]) { $nch++; $L[$i][0] = $new } }
  next unless $nch;
  $nfiles++; $nlines += $nch; push @changed, $r;
  if (exists $existing_new{$r}) {
    my @ln = sort { $a <=> $b } keys %$log;
    if (!@ln) { push @nolog, $r }
    else { my $k = $ln[-1]; my %pf;
      for my $j (@ln) { my $l=$L[$j-1][0]; $pf{$1}++ if $l !~ $ENTRY && $l =~ /^(%\s{2,})\S/ }
      my ($pre) = sort { $pf{$b} <=> $pf{$a} || $a cmp $b } keys %pf; $pre //= '%   ';
      my $eol = $L[$k-1][1] ne '' ? $L[$k-1][1] : "\n";
      splice @L, $k, 0, ["% MMP, $date: _sop dropped from the names of the executives, LPI", $eol],
                        ["${pre}and settings functions (2025 names where they exist); names only.", $eol];
      $L[$k-1][1] = $eol;
      push @d1, "$r\tafter line $k\tprefix '".$pre."'" } }
  if ($mode eq 'apply') { open(my $o,'>:raw',"$root/$r") or die "$r: $!";
    print $o join('', map { $_->[0].$_->[1] } @L); close $o } }

make_path($rep);
sub wr { my ($f,@x)=@_; open(my $o,'>',"$rep/$f") or die "$rep/$f: $!"; print $o map { "$_\n" } @x; close $o }
sub tot { my $h = shift; my $t = 0; $t += $_ for values %{$h||{}}; $t }
wr('tokens.tsv', map { "$_\t$stat{map}{$_}" } sort keys %{$stat{map}||{}});
wr('kept.tsv', (map { "keep\t$_\t$stat{keep}{$_}" } sort keys %{$stat{keep}||{}}),
               (map { "class\t$_\t$stat{class}{$_}" } sort keys %{$stat{class}||{}}));
wr('records.tsv', @records); wr('realign.tsv', @realign); wr('star.tsv', @star);
wr('d1.tsv', @d1); wr('changed.txt', @changed); wr('nolog.txt', @nolog);
wr('unknown.tsv', map { "$_\t$unknown{$_}" } sort keys %unknown);
my $sum = sprintf("mode %s: %d live .m files changed, %d lines; tokens mapped %d (%d kinds), kept %d, class %d; ".
  "PIETOOLS_*_sop %d; realigned %d; record lines with _sop %d; D1 lines %d; existing without a log %d; UNKNOWN %d",
  $mode, $nfiles, $nlines, tot($stat{map}), scalar(keys %{$stat{map}||{}}), tot($stat{keep}), tot($stat{class}),
  scalar(@star), scalar(@realign), scalar(@records), scalar(@d1), scalar(@nolog), scalar(keys %unknown));
wr('summary.txt', $sum); print "$sum\n";
exit(%unknown ? 1 : 0);
