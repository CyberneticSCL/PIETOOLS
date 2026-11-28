#!/usr/bin/perl
# SOP_CMP_TREE.PL ROOT EXPORT IGNORED [D_EXPORT D_ROOT] - 11/28/2026; cmp_tree.pl of
# ../../rename_copvar_opvar_1125/phase1 with orphaned_file/ excluded as well and the D1 text of
# this pass. Compares every file of the PIETOOLS folder ROOT (the checkout) with its copy in
# EXPORT, outside cleaning/ and orphaned_file/. Text files (no NUL byte) compare after CRLF -> LF,
# the conversion git applies on commit (core.autocrlf); other files byte for byte. IGNORED lists
# the git-ignored files of ROOT, which no export holds. With D_EXPORT and D_ROOT, the date of the
# D1 log line reads D_ROOT in EXPORT.
use strict; use warnings; use File::Find;
my ($root,$exp,$ign,$d0,$d1) = @ARGV;
die "usage: ROOT EXPORT IGNORED [D_EXPORT D_ROOT]\n" unless $root && $exp && $ign;
my %ign; open(my $h,'<',$ign) or die "$ign: $!"; while (<$h>) { s/\r?\n$//; $ign{$_} = 1 } close $h;
sub list { my $b = shift; my @f;
  find({ no_chdir=>1, wanted=>sub { return unless -f $_; (my $r = $_) =~ s#^\Q$b\E/##;
    push @f, $r unless $r =~ m#^(cleaning|orphaned_file)/# }}, $b); return @f }
sub rd { my $p = shift; open(my $f,'<:raw',$p) or die "$p: $!"; local $/; my $s = <$f>; close $f;
  $s //= ''; $s =~ s/\r\n/\n/g unless $s =~ /\0/; return $s }
my %R = map { $_ => 1 } grep { !$ign{$_} } list($root);
my %E = map { $_ => 1 } list($exp);
my ($same,$diff,$extra,$miss) = (0,0,0,0);
for my $r (sort keys %R) {
  if (!$E{$r}) { print "ONLY-ROOT $r\n"; $extra++; next }
  my ($x,$y) = (rd("$root/$r"), rd("$exp/$r"));
  $y =~ s/^(% MMP, )\Q$d0\E(: _sop dropped from the names of the executives, LPI)$/$1$d1$2/mg if defined $d1;
  if ($x eq $y) { $same++ } else { print "DIFF $r\n"; $diff++ } }
for my $r (sort keys %E) { next if $R{$r}; print "ONLY-EXPORT $r\n"; $miss++ }
print "same $same, differ $diff, only in ROOT $extra, only in EXPORT $miss\n";
exit(($diff || $extra || $miss) ? 1 : 0);
