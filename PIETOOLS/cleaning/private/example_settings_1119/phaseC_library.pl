#!/usr/bin/perl
# phaseC_library.pl REPO MODE [OUTFILE]
# The call sites of examples_PDE_library_PIETOOLS.m collect the output st of
# every example file whose signature returns it (read from the file when
# this script runs): '[PDE_t,PDE_b] = F(...)' becomes '[PDE_t,PDE_b,st] =
# F(...)', '[PDE_t] = F(...)' and 'PDE_t = F(...)' become '[PDE_t,st] =
# F(...)'. The old line is kept as a (was) comment above the new one, both
# marked % MMP, 11/19/2026 at column 77 (after the code when longer). MODE
# 'dry' writes OUTFILE, 'apply' rewrites the library. MMP, 11/19/2026.
use strict; use warnings;
my ($REPO, $MODE, $OUT) = @ARGV;  my $MK = '% MMP, 11/19/2026';
my $lib = "$REPO/PIETOOLS_examples/examples_PDE_library_PIETOOLS.m";
my %files; { my @f = `find "$REPO/PIETOOLS_examples/Examples_Library" -name 'PIETOOLS_PD*E_Ex*.m'`;
  for (@f) { chomp; my ($n) = m{/([^/]+)\.m$}; $files{$n} = $_; } }
sub returns_st { my ($n) = @_; my $p = $files{$n} or return 0; open(my $h, '<', $p) or return 0;
  my $l = <$h>; close $h; return $l =~ /^\s*function\s*\[[^\]]*\bst\s*\]/ ? 1 : 0; }
sub mk { my ($code, $tag) = @_; my $m = $tag ? "$MK $tag" : $MK; (my $v = $code) =~ s/\t/    /g; my $l = length($v);  # a tab counts 4 (the editor default)
  return $l < 76 ? $code . (' ' x (76 - $l)) . $m : "$code  $m"; }
open(my $h, '<:raw', $lib) or die; my $raw = do { local $/; <$h> }; close $h;
my $eol = $raw =~ /\r\n/ ? "\r\n" : "\n";  my $fe = $raw =~ /\r?\n\z/ ? 1 : 0;
my @L = split /\r?\n/, $raw, -1;  pop @L if @L && $L[-1] eq '';
my (@out, $n, @kept);
for my $l (@L) {
  if ($l !~ /^\s*%/ && $l =~ /^(\s*)(\[?PDE_t(?:,\s*PDE_b)?\]?)\s*=\s*(PIETOOLS_PDD?E_Ex\w*)(\(.*?\);)\s*(%.*)?$/) {
    my ($ind, $lhs, $fn, $args) = ($1, $2, $3, $4);
    if (returns_st($fn)) {
      (my $outs = $lhs) =~ s/[\[\]\s]//g;
      my $new = "$ind\[$outs,st\] = $fn$args";
      (my $old = $l) =~ s/\s+$//;  $old =~ s/\s*%\s*MMP,.*$//;
      push @out, mk('%' . $old, '(was)'), mk($new); $n++; next; }
    else { push @kept, $fn; } }
  push @out, $l; }
my $txt = join($eol, @out) . ($fe ? $eol : '');
if ($MODE eq 'apply') { open(my $w, '>:raw', $lib) or die; print $w $txt; close $w; }
else { open(my $w, '>:raw', $OUT) or die; print $w $txt; close $w; }
print "call sites changed: ", ($n // 0), "\n";
print "without st (unchanged): ", scalar(@kept), "\n", map { "  $_\n" } @kept;
