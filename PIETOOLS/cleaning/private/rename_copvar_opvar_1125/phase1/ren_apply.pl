#!/usr/bin/perl
# REN_APPLY.PL - Phase 1 of the rename copvar -> opvar, cdopvar -> dopvar (PLAN.md, Phase 1),
# 11/26/2026.
#   perl ren_apply.pl ROOT paths
#     prints "old<TAB>new" for each file outside the records whose path holds a mapped token.
#   perl ren_apply.pl ROOT check|apply DATE EXISTING REPORTDIR
#     maps every live .m file at its path after the moves; check writes the reports only, apply
#     also writes the files. DATE: mm/dd/yyyy of the D1 log line. EXISTING: the files of the 2025
#     release, by their path before the rename (PLAN.md, D1); each that changes gets the D1 line.
# ROOT is the PIETOOLS folder. Records, never edited: cleaning/, orphaned_file/, SOSTOOLS400/;
# in a live file, the header-log lines, the lines holding "(was)", and the retired names of
# @RETIRED (any case, also inside a longer token). On every other line of a live .m file:
#   1. "a"/"A" before a token that the map turns from cop... into op... becomes "an"/"An",
#      also when the article ends the comment of the line before;
#   2. each word token in %MAP becomes its new name; any other token that names a container
#      is UNKNOWN and stops apply;
#   3. a trailing comment that followed 2 or more blanks keeps its column, so that the change
#      markers stay at column 77.
use strict; use warnings; use File::Find; use File::Path qw(make_path);

my ($root,$mode,$date,$exfile,$rep) = @ARGV;
die "usage: ROOT paths | ROOT check|apply DATE EXISTING REPORTDIR\n"
  unless $root && $mode && ($mode eq 'paths' || ($mode =~ /^(check|apply)$/ && $date && $exfile && $rep));

# rename_table.tsv in lower case; upper-case and capitalized forms are derived
my @PAIRS = (
  [copvar=>'opvar'], [cdopvar=>'dopvar'], [poscopvar=>'poslpivar'], [lpivar_cdopvar=>'lpivar'],
  [mat2copvar_sop=>'mat2opvar'], [rand_copvar=>'rand_opvar'], [eye_copvar_sop=>'eye_opvar'],
  [zeros_copvar_sop=>'zeros_opvar'], [pi4_copvar_sop=>'pi4_opvar'], [int2copvar=>'int2opvar'],
  [qp_copvar=>'qp_opvar'], [lpi_eq_cdopvar=>'lpi_eq_dopvar'], [poscopvar_direct=>'poslpivar_direct'],
  [poscopvar_lift=>'poslpivar_lift'], [posmult_cdopvar=>'posmult_dopvar'], [copquadvar=>'opquadvar'],
  [lift_copvar=>'lift_opvar'], [cat_copvar_grid=>'cat_opvar_grid'],
  [derive_copvar_meta=>'derive_opvar_meta'], [dpvar_op_copvar=>'dpvar_op_opvar'],
  [eq_copvar=>'eq_opvar'], [mat2copvar_grid=>'mat2opvar_grid'],
  [merge_copvar_registry=>'merge_opvar_registry'], [split_copvar=>'split_opvar'],
  [verify_copvar_meta=>'verify_opvar_meta'], [copvar_space_list=>'opvar_space_list'],
  [parse_copvar_spaces=>'parse_opvar_spaces'], [copvar2tensopmat=>'opvar2tensopmat'],
  [tensopmat2copvar=>'tensopmat2opvar'], [mtimes_copvar=>'mtimes_opvar'],
  [poslpivar_sop=>'poslpivar_legacy'], [is_id_cdopvar=>'is_id_dopvar'],
  [check_copvar=>'check_opvar'], [maxabs_copvar=>'maxabs_opvar'], [apply_cop=>'apply_opvar'],
  [rand_cdopvar=>'rand_dopvar'], [test_copvar=>'test_opvar'], [test_cdopvar=>'test_dopvar'],
  [test_qp_copvar_nd=>'test_qp_opvar_nd'], [test_copquadvar_faces=>'test_opquadvar_faces']);
my @RETIRED = qw(op2copvar_sop pie2copvar opvar2copvar copvar2opvar opvar2d2copvar cg_copvar_pd
                 test_cg_copvar_nd test_copvar_kyp test_poscopvar_vs_poslpivar_2d);
my $NAMES = qr/c(?:d)?opvar|copquadvar|apply_cop|poslpivar_sop/i;   # a token that names a container
my (%MAP,@ART);
for my $p (@PAIRS) { my ($o,$n) = @$p;
  $MAP{$o} = $n; $MAP{uc $o} = uc $n; $MAP{ucfirst $o} = ucfirst $n }
for my $o (keys %MAP) { push @ART, $o if $o =~ /^cop/i && $MAP{$o} =~ /^op/i }
my $ART = join '|', sort { length($b) <=> length($a) } @ART;

my (%stat,%unknown,@art,@realign,@records,@d1,@changed,@nolog);
sub maptok { my ($t,$where) = @_;
  return $t unless $t =~ $NAMES;
  my $lt = lc $t;
  for my $r (@RETIRED) { if (index($lt,$r) >= 0) { $stat{retired}{$t}++; return $t } }
  if (exists $MAP{$t}) { $stat{map}{"$t\t$MAP{$t}"}++; return $MAP{$t} }
  $unknown{"$t\t$where"}++; return $t }
sub mappath { my $p = shift; join '/', map { (my $c=$_) =~ s/(\w+)/maptok($1,"path $p")/ge; $c } split m#/#, $p }

# code and comment of a line; the comment starts at % or ... outside strings (h2scan.pl rule)
sub split_line { my $l=shift; my ($o,$i,$n,$p)=('',0,length $l,'');
  while ($i<$n) { my $c=substr($l,$i,1);
    return ($o, substr($l,$i)) if $c eq '%' || substr($l,$i,3) eq '...';
    if ($c eq '"' || ($c eq "'" && $p !~ /[\w\)\]\}\.']/)) { my $q=$c; my $j=$i+1;
      while ($j<$n) { my $d=substr($l,$j,1); if ($d eq $q) { if ($j+1<$n && substr($l,$j+1,1) eq $q) { $j+=2; next } last } $j++ }
      $o.=substr($l,$i,$j-$i+1); $i=$j+1; $p=$q; next }
    $o.=$c; $p = ($c eq ' ' || $c eq "\t") ? ' ' : $c; $i++ }
  return ($o, '') }

# header-log lines: before the first code line (the function or classdef line, its continuations
# and a code line ahead of the header in the first 3 lines do not count), an entry line and the
# continuation lines (2 or more blanks after %) that follow it
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

sub mapseg { my ($s,$where) = @_;
  if ($s =~ s/\b([Aa])(\s+['"`(]?)(?=(?:$ART)\b)/$1n$2/g) { push @art, "$where\tsame line" }
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
  next unless grep { /$NAMES/ } @T;
  my $log = header_log(\@T);
  my %rec; for my $i (0..$#T) {
    if ($log->{$i+1}) { $rec{$i}='log' } elsif ($T[$i] =~ /\(was\)/) { $rec{$i}='was' } }
  # 1. the article across a line break: "... a" + "% 'copvar' ..."
  for my $i (0..$#T-1) { next if $rec{$i} || $rec{$i+1};
    my (undef,$cmt) = split_line($T[$i]);
    if ($cmt =~ /^%/ && $cmt =~ /\b[Aa]\s*$/ && $T[$i+1] =~ /^\s*%+\s*['"`(]?(?:$ART)\b/) {
      $T[$i] =~ s/\b([Aa])(\s*)$/$1n$2/; push @art, "$r:".($i+1)."\tline break" } }
  my $nch = 0;
  for my $i (0..$#T) { my $t=$T[$i];
    if ($rec{$i}) { push @records, "$r:".($i+1)."\t$rec{$i}\t$t" if $t =~ $NAMES; next }
    my ($code,$cmt) = split_line($t); my $where = "$r:".($i+1);
    my $nc = mapseg($code,$where); my $nm = mapseg($cmt,$where);
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
      splice @L, $k, 0, ["% MMP, $date: copvar, cdopvar renamed opvar, dopvar, with the", $eol],
                        ["${pre}functions of the rename table; names only.", $eol];
      $L[$k-1][1] = $eol;
      push @d1, "$r\tafter line $k\tprefix '".$pre."'" } }
  if ($mode eq 'apply') { open(my $o,'>:raw',"$root/$r") or die "$r: $!";
    print $o join('', map { $_->[0].$_->[1] } @L); close $o } }

make_path($rep);
sub wr { my ($f,@x)=@_; open(my $o,'>',"$rep/$f") or die "$rep/$f: $!"; print $o map { "$_\n" } @x; close $o }
wr('tokens.tsv', map { "$_\t$stat{map}{$_}" } sort keys %{$stat{map}||{}});
wr('retired.tsv', map { "$_\t$stat{retired}{$_}" } sort keys %{$stat{retired}||{}});
wr('records.tsv', @records); wr('article.tsv', @art); wr('realign.tsv', @realign);
wr('d1.tsv', @d1); wr('changed.txt', @changed); wr('nolog.txt', @nolog);
wr('unknown.tsv', map { "$_\t$unknown{$_}" } sort keys %unknown);
my $sum = sprintf("mode %s: %d live .m files changed, %d lines; tokens mapped %d (%d kinds), retired kept %d; ".
  "article %d; realigned %d; record lines with a name %d; D1 lines %d; existing without a log %d; UNKNOWN %d",
  $mode, $nfiles, $nlines, eval { my $t=0; $t+=$_ for values %{$stat{map}||{}}; $t }, scalar(keys %{$stat{map}||{}}),
  eval { my $t=0; $t+=$_ for values %{$stat{retired}||{}}; $t }, scalar(@art), scalar(@realign),
  scalar(@records), scalar(@d1), scalar(@nolog), scalar(keys %unknown));
wr('summary.txt', $sum); print "$sum\n";
exit(%unknown ? 1 : 0);
