#!/usr/bin/perl
# phaseC_blocks.pl DIR REPO HARNESS MODE [OUTDIR]
# Writes the example-specific settings into the commented blocks of
# examples_DDE_library_PIETOOLS.m and examples_NDSDDF_library_PIETOOLS.m,
# from the plan DIR/phaseC_plan.csv (groups 'DDE:<id>', 'NDSDDF:<id>' with
# status 'ready'). A block is found by its code: the lines of the harness
# copy HARNESS/blocks/<id>.m (one comment level removed) appear in the
# library with '%' in front, in order. In each block, after its last code
# line (and after any settings lines of 11/18/2026):
#   % % st: <form, result, size>                         (a note)
#   % st = lpisettings_sop('<tier>');                    (and its changes)
#   % stability_form = 'P';  or  % Hinf_gain_form = 'P'; (the forms named)
# all commented one level, as the block is, so they run when the user
# uncomments the block. Settings lines of 11/18/2026 in the block become
# (was) lines with two comment levels, so they stay comments then. Lines
# carry % MMP, 11/19/2026 at column 77. MODE 'dry' writes the libraries to
# OUTDIR, 'apply' in place. MMP, 11/19/2026.
# Not idempotent: apply once. It writes every block one comment level down; the
# block active in the DDE library (Example 2b) was set to level 0 by hand,
# 11/19/2026 04:50.
use strict; use warnings;
my ($D, $REPO, $HB, $MODE, $OUT) = @ARGV;  my $MK = '% MMP, 11/19/2026';
sub csv { my ($f) = @_; my @r; open(my $h, '<', $f) or return (); my $hd = <$h>; $hd =~ s/\r?\n$//;
  my @c = split /,/, $hd;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = ($l =~ /("(?:[^"]|"")*"|[^,]*)(?:,|$)/g);
    s/^"(.*)"$/$1/s for @x; s/""/"/g for @x; my %e; @e{@c} = @x; push @r, \%e; } return @r; }
sub mk { my ($code, $tag) = @_; my $m = $tag ? "$MK $tag" : $MK; my $l = length($code);
  return $l < 76 ? $code . (' ' x (76 - $l)) . $m : "$code  $m"; }
sub fmt { my $n = shift; return $n unless defined $n && $n =~ /^\d+$/; 1 while $n =~ s/^(\d+)(\d{3})/$1,$2/; return $n; }
my @plan = grep { $_->{group} =~ /^(DDE|NDSDDF):/ && $_->{status} eq 'ready' } csv("$D/phaseC_plan.csv");
my %num = map { $_->{key} => $_ } csv("$D/phaseC_numbers.csv");
my %lib = (DDE => "$REPO/PIETOOLS_examples/examples_DDE_library_PIETOOLS.m",
           NDSDDF => "$REPO/PIETOOLS_examples/examples_NDSDDF_library_PIETOOLS.m");
for my $which ('DDE', 'NDSDDF') {
  my $path = $lib{$which};
  open(my $h, '<:raw', $path) or die; my $raw = do { local $/; <$h> }; close $h;
  my $eol = $raw =~ /\r\n/ ? "\r\n" : "\n";  my $fe = $raw =~ /\r?\n\z/ ? 1 : 0;
  my @L = split /\r?\n/, $raw, -1;  pop @L if @L && $L[-1] eq '';
  my (@edits, @log);
  for my $p (grep { $_->{group} =~ /^\Q$which\E:/ } @plan) {
    my ($id) = $p->{group} =~ /:(\S+)$/;
    my $bf = "$HB/blocks/$id.m";  unless (-f $bf) { push @log, "$id: no harness copy"; next }
    open(my $b, '<', $bf) or die; my @B = <$b>; close $b; shift @B;
    my @code = grep { $_ !~ /^\s*%/ && $_ =~ /\S/ } map { s/\r?\n$//r } @B;
    unless (@code) { push @log, "$id: no code lines"; next }
    # the library lines of the block code, normalized as lb_blocks reads them
    # (the first '%' after the leading blanks removed), in order, with only
    # prose (a second comment level) or blank lines between them
    my @S = map { $_ =~ s/^(\s*)%/$1/r } @L;
    my @at;
    for my $i (0..$#L) {
      next unless $S[$i] eq $code[0];
      my @hit = ($i); my $j = $i + 1; my $ok = 1;
      for my $c (@code[1..$#code]) {
        while ($j <= $#L && $S[$j] ne $c && ($S[$j] =~ /^\s*%/ || $S[$j] =~ /^\s*$/)) { $j++ }
        if ($j > $#L || $S[$j] ne $c) { $ok = 0; last } push @hit, $j; $j++; }
      push @at, [@hit] if $ok; }
    if (!@at) { push @log, "$id: block not found"; next }
    # a block found more than once is the same code repeated (DDF02 repeats
    # DDF01, lb_meta): every copy gets the same settings
    for my $hit (@at) {
    my $last = $hit->[-1];
    # settings lines of 11/18/2026 right after the code: (was), two levels
    my $k = $last + 1; my @was;
    while ($k <= $#L && $L[$k] =~ /MMP, 11\/18\/2026/ && $L[$k] =~ /^\s*%\s*(?:%\s*)?st\b/) { push @was, $k; $k++ }
    # the new lines
    my @rows = split / /, $p->{rows};
    my @stm = grep { $_ ne '' } map { s/^\s+|\s+$//gr } split /(?<=;)\s+(?=st\.)/, $p->{statements};
    my ($sf, $gf) = ($p->{stability_form}, $p->{Hinf_gain_form});
    my @res; for my $r (@rows) { my $n = $num{$r}; next unless $n && $n->{result} ne '';
      push @res, $n->{result} . ($r =~ m{/(\w+)$} ? " ($1)" : ''); }
    my $n0 = $num{$rows[0]}; my $nd = $n0 ? $n0->{ndv} : $p->{ndv};
    my @fm; push @fm, "form '$sf'" if $sf ne ''; push @fm, "gain form '$gf'" if $gf ne '';
    my $note = "st: example-specific settings (11/19/2026)" . (@fm ? ', ' . join(', ', @fm) : '')
      . (@res ? '; ' . join(', ', @res) : '') . ', ' . fmt($nd) . ' decision variables';
    # the indentation of the block; the note wrapped at 70 characters
    my ($ind) = $L[$last] =~ /^(\s*)/;
    my (@nl, $cur); $cur = '';
    for my $w (split / /, $note) { if ($cur ne '' && length($cur) + length($w) + 1 > 66 - length($ind)) { push @nl, $cur; $cur = $w } else { $cur = $cur eq '' ? $w : "$cur $w" } }
    push @nl, $cur if $cur ne '';
    my @new = ((map { mk("$ind% % " . ($_ == 0 ? '' : '    ') . $nl[$_]) } 0..$#nl),
               mk("$ind% st = lpisettings_sop('$p->{tier}');"), map { mk("$ind% $_") } @stm);
    push @new, mk("$ind% stability_form = '$sf';") if $sf ne '';
    push @new, mk("$ind% Hinf_gain_form = '$gf';") if $gf ne '';
    push @edits, [$k, \@was, \@new, $id . (@at > 1 ? " (copy at line " . ($hit->[0] + 1) . ")" : '')];
    }
  }
  # apply from the bottom up
  for my $e (sort { $b->[0] <=> $a->[0] } @edits) {
    my ($k, $was, $new, $id) = @$e;
    for my $w (@$was) { (my $t = $L[$w]) =~ s/\s*%\s*MMP, 11\/18\/2026\s*$//; $t =~ s/^(\s*)%/$1% %/; $L[$w] = mk($t, '(was)'); }
    splice @L, $k, 0, @$new;
    push @log, "$id: " . scalar(@$new) . " lines" . (@$was ? ", " . scalar(@$was) . " (was)" : ''); }
  # header log entry after the last log entry of 11/18/2026 (or 11/xx) in the header
  if (@edits) {
    my ($last) = reverse grep { $L[$_] =~ m{^% MMP, \d\d/\d\d/2026:} } 0..120;
    ($last) = reverse grep { $L[$_] =~ m{^% (?:Initial coding|[A-Z]{2,3}[ ,-]+\d\d[/_]\d\d[/_]\d{2,4})} } 0..120 unless defined $last;
    my $nwas = 0;  $nwas += scalar(@{$_->[1]}) for @edits;
    if (defined $last) { my $j = $last + 1; $j++ while $j <= $#L && $L[$j] =~ /^%\s{3,}\S/;
      my @h = ('% MMP, 11/19/2026: Example-specific settings st and the forms of the tests',
        '%                   (stability_form, Hinf_gain_form) in every block with a',
        '%                   certified setting, from the per-example search of',
        '%                   11/19/2026 (maintainer: no default settings);');
      if ($nwas) { push @h, '%                   the st lines of 11/18/2026 are kept as (was) lines with',
        '%                   two comment levels, so they stay comments when a block',
        '%                   is uncommented;'; }
      splice @L, $j, 0, @h; }
  }
  my $txt = join($eol, @L) . ($fe ? $eol : '');
  my $name = $path =~ s{.*/}{}r;
  if ($MODE eq 'apply') { open(my $w, '>:raw', $path) or die; print $w $txt; close $w; }
  else { open(my $w, '>:raw', "$OUT/$name") or die; print $w $txt; close $w; }
  print "$name: ", scalar(@edits), " blocks\n", map { "  $_\n" } @log;
}
