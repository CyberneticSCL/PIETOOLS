#!/usr/bin/perl
# phaseC_write.pl DIR REPO MODE [OUTDIR]
# Writes the example-specific settings into the PDE example files of the
# plan DIR/phaseC_plan.csv (groups with status 'ready', one file, no index,
# no output st yet), with the annotation protocol of PIETOOLS/CLAUDE.md:
#   - output st added to the signature, the old line kept as (was) in a
#     header log entry;
#   - an OUTPUT line for st;
#   - a note in the description: the form, the settings, the measured result
#     (DIR/phaseC_numbers.csv when present: key, result, ndv, m, ndv0, tier0);
#   - the form of the tests named next to the executive flags
#     (evalin('base','stability_form = ...') or Hinf_gain_form);
#   - the st block before the GUI branch.
# Changed lines carry % MMP, 11/19/2026 at column 77 (after the code when it
# is longer). MODE 'dry' writes the edited files to OUTDIR; 'apply' writes
# them in place. Lists the files edited and skipped. MMP, 11/19/2026.
use strict; use warnings;
my ($D, $REPO, $MODE, $OUT) = @ARGV;  my $MK = '% MMP, 11/19/2026';
sub csv { my ($f) = @_; my @r; open(my $h, '<', $f) or return (); my $hd = <$h>; $hd =~ s/\r?\n$//;
  my @c = split /,/, $hd;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = ($l =~ /("(?:[^"]|"")*"|[^,]*)(?:,|$)/g);
    s/^"(.*)"$/$1/s for @x; s/""/"/g for @x; my %e; @e{@c} = @x; push @r, \%e; } return @r; }
my @plan = csv("$D/phaseC_plan.csv");
my %num = map { $_->{key} => $_ } csv("$D/phaseC_numbers.csv");
my %files; { my @f = `find "$REPO/PIETOOLS_examples/Examples_Library" -name 'PIETOOLS_PDE_Ex*.m'`;
  for (@f) { chomp; my ($n) = m{/([^/]+)\.m$}; $files{$n} = $_; } }
sub mk { my ($code) = @_; my $l = length($code);
  return $l < 76 ? $code . (' ' x (76 - $l)) . $MK : "$code  $MK"; }
sub wrapnote { my ($txt, $pre) = @_; $pre //= '% % '; my @w = split / /, $txt; my (@L, $cur); $cur = '';
  for my $w (@w) { if (length($cur) + length($w) + 1 > 70) { push @L, $cur; $cur = $w } else { $cur = $cur eq '' ? $w : "$cur $w" } }
  push @L, $cur if $cur ne ''; return map { mk("$pre$_") } @L; }
sub fmt { my $n = shift; return $n unless defined $n && $n =~ /^\d+$/; 1 while $n =~ s/^(\d+)(\d{3})/$1,$2/; return $n; }
my (@done, @skip);
for my $p (@plan) {
  my $file = $p->{file};
  if ($p->{status} ne 'ready')       { push @skip, "$p->{group}: status $p->{status}"; next }
  if ($p->{index} ne '')             { push @skip, "$p->{group}: indexed (by hand)"; next }
  if ($file !~ /^PIETOOLS_PDE_Ex/)  { push @skip, "$p->{group}: library block (by hand)"; next }
  my $path = $files{$file};          unless ($path) { push @skip, "$p->{group}: file not found"; next }
  open(my $h, '<:raw', $path) or die; my $raw = do { local $/; <$h> }; close $h;
  my $eol = $raw =~ /\r\n/ ? "\r\n" : "\n";
  my $final_eol = $raw =~ /\r?\n\z/ ? 1 : 0;      # keep the end of file as it is
  my @L = split /\r?\n/, $raw, -1; pop @L if @L && $L[-1] eq '';
  if ($L[0] =~ /\bst\s*\]/)          { push @skip, "$p->{group}: has an output st (by hand)"; next }
  # 1. signature
  my $old = $L[0]; (my $sig = $old) =~ s/\s+%.*$//;
  my $new;
  if ($sig =~ /^function\s*\[([^\]]*)\]\s*=\s*(.*)$/) { $new = "function [$1,st] = $2" }
  elsif ($sig =~ /^function\s+(\w+)\s*=\s*(.*)$/)       { $new = "function [$1,st] = $2" }
  else { push @skip, "$p->{group}: signature not recognized"; next }
  $L[0] = mk($new);
  # rows, forms, numbers
  my @rows = split / /, $p->{rows}; my $ex = $p->{executives};
  my ($sf, $gf) = ($p->{stability_form}, $p->{Hinf_gain_form});
  # statements end with ';' and may hold matrices with ';' ([2;1;1]): split
  # only where the next assignment to st begins
  my @stm = grep { $_ ne '' } map { s/^\s+|\s+$//gr } split /(?<=;)\s+(?=st\.)/, $p->{statements};
  my $chg = @stm ? join(', ', map { s/;$//r } @stm) : 'no change';
  my @res;
  for my $k (@rows) { my $n = $num{$k}; next unless $n;
    push @res, ($n->{result} ne '' ? "$n->{result}" : 'certified') . ' (' . ($k =~ m{/(.*)} ? $1 : 'primal') . ')'; }
  my $n0 = $num{$rows[0]}; my ($nd, $mm, $nd0, $t0) = $n0 ? ($n0->{ndv}, $n0->{m}, $n0->{ndv0}, $n0->{tier0}) : ($p->{ndv}, '', '', '');
  my @forms; push @forms, "form '$sf' for the stability tests" if $sf ne ''; push @forms, "form '$gf' for the H-infinity gain" if $gf ne '';
  push @forms, "form 'o' for the H2 norm" if grep { /^\s*evalin\('base','(H2_norm|H2_norm_dual)\s*=\s*1;'\)/ } @L;
  my $note = "Example-specific settings (output st; maintainer, 11/19/2026: each example names its settings): "
    . (@forms ? join(' and ', @forms) . '; ' : '')
    . "lpisettings_sop('$p->{tier}') with $chg. MOSEK, 11/19/2026: "
    . (@res ? join(', ', @res) : 'certified') . " at " . fmt($nd) . " decision variables"
    . ($mm ne '' ? ", " . fmt($mm) . " rows" : '')
    . ($nd0 =~ /^\d+$/ && $nd =~ /^\d+$/ && $nd0 > $nd ? " ($t0 as is: " . fmt($nd0) . ")" : '')
    . ". Measured in cleaning/private/example_settings_1119.";
  # 2. OUTPUT line for st, after the last output item of the help
  my ($io) = grep { $L[$_] =~ /^%\s*OUTPUT/ } 0..$#L;
  if (defined $io) { my $j = $io + 1; my $last = $j;
    while ($j <= $#L && $L[$j] =~ /^%/ && $L[$j] !~ /^%\s*$/ && $L[$j] !~ /^% %-/) { $last = $j if $L[$j] =~ /^%\s+-\s|^%\s{6,}/; $j++ }
    splice @L, $last + 1, 0, mk('% - st:         example-specific settings (lpisettings_sop structure)'),
                              mk('%               for the tests this file names (header)'); }
  # 3. note before the closing separator of the description
  #    (a description with one separator ends at the comment 'Determine the
  #    location' of the code: the note goes before it, in the plain style)
  my @sep = grep { $L[$_] =~ /^% %-{10,}/ } 0..$#L;
  my ($c0) = grep { $_ > 0 && $L[$_] !~ /^\s*%/ && $L[$_] !~ /^\s*$/ } 0..$#L;
  my ($dl) = grep { $_ < $c0 && $L[$_] =~ /^%\s*Determine the location/ } 0..$#L;
  if (@sep >= 2) { splice @L, $sep[1], 0, mk('% %'), wrapnote($note); }
  elsif (defined $dl) { splice @L, $dl, 0, mk('%'), wrapnote($note, '% '); }
  else { push @skip, "$p->{group}: no place for the note"; next }
  # 4. header log entry right after the description block (its closing
  #    separator, or the note when it has one), else before the first code
  #    line after the header
  @sep = grep { $L[$_] =~ /^% %-{10,}/ } 0..$#L;
  ($c0) = grep { $_ > 0 && $L[$_] !~ /^\s*%/ && $L[$_] !~ /^\s*$/ } 0..$#L;
  ($dl) = grep { $_ < $c0 && $L[$_] =~ /^%\s*Determine the location/ } 0..$#L;
  my $at = @sep >= 2 ? $sep[1] + 1 : (defined $dl ? $dl : $c0);
  # an existing header log (entries '% XY, mm/dd/yyyy: ...' before the code):
  # the new entry goes after its last entry and that entry's continuation
  my @lg = grep { $_ < $c0 && $L[$_] =~ m{^%\s*[A-Z]{2,4}(?:, [A-Z]{2,4})*,? \d{1,2}[/_]\d{1,2}[/_]\d{2,4}} } 0..$#L;
  if (@lg) { $at = $lg[-1] + 1; $at++ while $at < $c0 && $L[$at] =~ /^%\s{2,}\S/; }
  my $lead = ($at > 0 && ($L[$at-1] =~ /^% %-{10,}/ || (@sep < 2 && $L[$at-1] =~ /\Q$MK\E$/))) ? 1 : 0;
  my @log = ('% MMP, 11/19/2026: Output st, the example-specific settings, and the form',
             '%                   of its tests named in the workspace (stability_form,',
             '%                   Hinf_gain_form, H2_norm_form), from the search of',
             '%                   11/19/2026 (maintainer: no default settings);',
             mk('% ' . $sig =~ s/\s+$//r) =~ s/\Q$MK\E$/$MK (was)/r);
  unshift @log, '%' if $lead;              # a blank comment line after the separator
  splice @L, $at, 0, @log;
  # 5. forms after the last active flag line of each kind
  my @fl = grep { $L[$_] =~ /^\s*evalin\('base','(stability|stability_dual)\s*=\s*1;'\)/ } 0..$#L;
  my @gl = grep { $L[$_] =~ /^\s*evalin\('base','(Hinf_gain|Hinf_gain_dual)\s*=\s*1;'\)/ } 0..$#L;
  my @hl = grep { $L[$_] =~ /^\s*evalin\('base','(H2_norm|H2_norm_dual)\s*=\s*1;'\)/ } 0..$#L;
  my @ins;
  push @ins, [$hl[-1] + 1, mk("evalin('base','H2_norm_form = ''o'';');  % form of the H2 norm test")] if @hl;
  push @ins, [$gl[-1] + 1, mk("evalin('base','Hinf_gain_form = ''$gf'';');  % form of the gain test")] if @gl && $gf ne '';
  push @ins, [$fl[-1] + 1, mk("evalin('base','stability_form = ''$sf'';');  % form of the stability test")] if @fl && $sf ne '';
  for my $i (sort { $b->[0] <=> $a->[0] } @ins) { splice @L, $i->[0], 0, $i->[1]; }
  # 6. st block before the last GUI branch, else before the final end
  my @gui = grep { $L[$_] =~ /^\s*if\s+GUI\s*~=\s*0/ } 0..$#L;
  my $bt = @gui ? $gui[-1] : (grep { $L[$_] =~ /^end\s*$/ } 0..$#L)[-1];
  my @blk = (mk('% Example-specific settings (header)'), mk("st = lpisettings_sop('$p->{tier}');"), map { mk($_ =~ /;$/ ? $_ : "$_;") } @stm);
  splice @L, $bt, 0, @blk, '';
  my $txt = join($eol, @L) . ($final_eol ? $eol : '');
  if ($MODE eq 'apply') { open(my $w, '>:raw', $path) or die; print $w $txt; close $w; }
  else { open(my $w, '>:raw', "$OUT/$file.m") or die; print $w $txt; close $w; }
  push @done, "$p->{group} ($path)";
}
print "edited ", scalar(@done), "\n", map { "  $_\n" } @done;
print "skipped ", scalar(@skip), "\n", map { "  $_\n" } @skip;
