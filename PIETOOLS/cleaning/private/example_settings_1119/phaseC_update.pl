#!/usr/bin/perl
# phaseC_update.pl DIR REPO MODE [OUTDIR]
# The companion of phaseC_write.pl for the PDE example files that already
# return st (the example settings of 11/17-18/2026): updates their st block
# to the plan DIR/phaseC_plan.csv, with the annotation protocol of
# PIETOOLS/CLAUDE.md:
#   - in the block 'if nargout>1 ... end', a changed statement is kept as a
#     (was) line, its new form follows; a new statement is added before 'end';
#     statements the plan keeps are left as they are;
#   - a note in the description (form, settings, measured result);
#   - the form of the tests named next to the executive flags;
#   - a header log entry after the last entry of the header log.
# Files that already carry an 11/19/2026 entry are skipped. MODE 'dry' writes
# to OUTDIR, 'apply' in place. MMP, 11/19/2026.
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
sub mk { my ($code, $tag) = @_; $tag //= $MK; my $l = length($code);
  return $l < 76 ? $code . (' ' x (76 - $l)) . $tag : "$code  $tag"; }
sub wrapnote { my ($txt) = @_; my @w = split / /, $txt; my (@L, $cur); $cur = '';
  for my $w (@w) { if (length($cur) + length($w) + 1 > 70) { push @L, $cur; $cur = $w } else { $cur = $cur eq '' ? $w : "$cur $w" } }
  push @L, $cur if $cur ne ''; return map { mk("% % $_") } @L; }
sub fmt { my $n = shift; return $n unless defined $n && $n =~ /^\d+$/; 1 while $n =~ s/^(\d+)(\d{3})/$1,$2/; return $n; }
sub lhs { my $s = shift; $s =~ s/^\s+//; $s =~ s/\s*=.*$//s; return $s; }   # 'st.lpi.ineq.dD'
my (@done, @skip);
for my $p (@plan) {
  my $file = $p->{file};
  next unless $p->{status} eq 'ready' && $p->{index} eq '' && $file =~ /^PIETOOLS_PDE_Ex/;
  my $path = $files{$file};  next unless $path;
  open(my $h, '<:raw', $path) or die; my $raw = do { local $/; <$h> }; close $h;
  next if $raw =~ /MMP, 11\/19\/2026/;                       # done (generator or this script)
  my $eol = $raw =~ /\r\n/ ? "\r\n" : "\n";  my $final_eol = $raw =~ /\r?\n\z/ ? 1 : 0;
  my @L = split /\r?\n/, $raw, -1; pop @L if @L && $L[-1] eq '';
  next unless $L[0] =~ /\bst\s*\]/;
  # the st block: 'st = lpisettings_sop(...)' inside 'if nargout>1 ... end'
  my ($b0) = grep { $L[$_] =~ /^\s*st\s*=\s*lpisettings_sop\('(\w+)'\)\s*;/ } 0..$#L;
  unless (defined $b0) { push @skip, "$p->{group}: no st block"; next }
  my ($b1) = grep { $_ > $b0 && $L[$_] =~ /^\s*end\b/ } 0..$#L;
  my ($ind) = $L[$b0] =~ /^(\s*)/;
  my %old;                                                   # lhs -> line index
  for my $i ($b0 + 1 .. $b1 - 1) { next if $L[$i] =~ /^\s*%/; $old{lhs($L[$i])} = $i if $L[$i] =~ /^\s*st\./; }
  my ($tier_old) = $L[$b0] =~ /lpisettings_sop\('(\w+)'\)/;
  my @stm = grep { $_ ne '' } map { s/^\s+|\s+$//gr } split /(?<=;)\s+(?=st\.)/, $p->{statements};
  my %new = map { lhs($_) => $_ } @stm;
  # statements of the old block the plan drops: the tier default returns
  my @drop = grep { !exists $new{$_} } sort keys %old;
  my (@chg, @add, %rep);
  for my $s (@stm) { my $k = lhs($s); (my $code = $s) =~ s/;?$/;/;
    if (exists $old{$k}) { (my $o = $L[$old{$k}]) =~ s/\s+%.*$//; $o =~ s/^\s+//;
      $rep{$old{$k}} = $code if $o ne $code; }
    else { push @add, $code } }
  my $tchg = $tier_old ne $p->{tier};
  # the block, bottom up: added lines before 'end', replaced and dropped lines as (was)
  splice @L, $b1, 0, map { mk("$ind$_") } @add;
  for my $i (sort { $b <=> $a } (keys %rep, map { $old{$_} } @drop)) {
    (my $o = $L[$i]) =~ s/\s+%\s*MMP, (\d\d\/\d\d\/\d{4})\s*$//;  my $od = $1 // '';
    my $was = mk('%' . substr($o, 1), "$MK (was)" . ($od ne '' ? "; MMP, $od" : ''));
    splice @L, $i, 1, $was, (exists $rep{$i} ? mk("$ind$rep{$i}") : ());
  }
  if ($tchg) { (my $o = $L[$b0]) =~ s/\s+%\s*MMP, (\d\d\/\d\d\/\d{4})\s*$//;  my $od = $1 // '';
    splice @L, $b0, 1, mk('%' . substr($o, 1), "$MK (was)" . ($od ne '' ? "; MMP, $od" : '')),
                       mk("${ind}st = lpisettings_sop('$p->{tier}');"); }
  my @what; push @what, "tier '$p->{tier}'" if $tchg; push @what, map { s/;$//r } (values %rep, @add);
  push @what, map { "$_ back to the tier value" } @drop;
  my $what = @what ? join(', ', @what) : 'settings unchanged';
  # the header comment of the block
  my ($hc) = grep { $_ < $b0 && $_ > $b0 - 8 && $L[$_] =~ /^% Example-specific settings \(header\)/ } 0..$#L;
  ($hc) = grep { $_ < $b0 && $L[$_] =~ /^\s*if\s+nargout\s*>\s*1/ } reverse 0..$#L unless defined $hc;
  if (defined $hc) { my $j = $hc; $j++ while $j + 1 < $#L && $L[$j+1] =~ /^%/;
    splice @L, $j + 1, 0, mk("% 11/19/2026: $what (per-example search)") if @what; }
  # note before the closing separator of the description
  my ($sf, $gf) = ($p->{stability_form}, $p->{Hinf_gain_form});
  my @rows = split / /, $p->{rows};  my @res;
  for my $k (@rows) { my $n = $num{$k}; next unless $n;
    push @res, ($n->{result} ne '' ? "$n->{result}" : 'certified') . ' (' . ($k =~ m{/(.*)} ? $1 : 'primal') . ')'; }
  my $n0 = $num{$rows[0]}; my ($nd, $mm) = $n0 ? ($n0->{ndv}, $n0->{m}) : ($p->{ndv}, '');
  my @forms; push @forms, "form '$sf' for the stability tests" if $sf ne ''; push @forms, "form '$gf' for the H-infinity gain" if $gf ne '';
  my $chg = @stm ? join(', ', map { s/;$//r } @stm) : 'no change';
  my $note = "Example-specific settings of 11/19/2026 (output st, the per-example search; replaces the st above): "
    . (@forms ? join(' and ', @forms) . '; ' : '')
    . "lpisettings_sop('$p->{tier}') with $chg. MOSEK, 11/19/2026: "
    . (@res ? join(', ', @res) : 'certified') . " at " . fmt($nd) . " decision variables"
    . ($mm ne '' ? ", " . fmt($mm) . " rows" : '') . ". Measured in cleaning/private/example_settings_1119.";
  my @sep = grep { $L[$_] =~ /^% %-{10,}/ } 0..$#L;
  if (@sep >= 2) { splice @L, $sep[1], 0, mk('% %'), wrapnote($note); } else { push @skip, "$p->{group}: no description separator"; next }
  # header log entry after the last entry before the code
  my ($c0) = grep { $_ > 0 && $L[$_] !~ /^\s*%/ && $L[$_] !~ /^\s*$/ } 0..$#L;
  my @lg = grep { $_ < $c0 && $L[$_] =~ m{^%\s*[A-Z]{2,4}(?:, [A-Z]{2,4})*,? \d{1,2}[/_]\d{1,2}[/_]\d{2,4}} } 0..$#L;
  unless (@lg) { push @skip, "$p->{group}: no header log"; next }
  my $at = $lg[-1] + 1; $at++ while $at < $c0 && $L[$at] =~ /^%\s{2,}\S/;
  splice @L, $at, 0, '% MMP, 11/19/2026: st updated to the per-example search of 11/19/2026',
                     '%                   (maintainer: no default settings): ' . (length($what) < 40 ? $what : 'see the block') . ';',
                     '%                   the form of its tests named in the workspace.';
  # forms after the last active flag line of each kind
  my @fl = grep { $L[$_] =~ /^\s*evalin\('base','(stability|stability_dual)\s*=\s*1;'\)/ } 0..$#L;
  my @gl = grep { $L[$_] =~ /^\s*evalin\('base','(Hinf_gain|Hinf_gain_dual)\s*=\s*1;'\)/ } 0..$#L;
  my @ins;
  push @ins, [$gl[-1] + 1, mk("evalin('base','Hinf_gain_form = ''$gf'';');  % form of the gain test")] if @gl && $gf ne '';
  push @ins, [$fl[-1] + 1, mk("evalin('base','stability_form = ''$sf'';');  % form of the stability test")] if @fl && $sf ne '';
  for my $i (sort { $b->[0] <=> $a->[0] } @ins) { splice @L, $i->[0], 0, $i->[1]; }
  my $txt = join($eol, @L) . ($final_eol ? $eol : '');
  my $dst = $MODE eq 'apply' ? $path : "$OUT/$file.m";
  open(my $w, '>:raw', $dst) or die; print $w $txt; close $w;
  push @done, "$p->{group}: $what";
}
print "updated ", scalar(@done), "\n", map { "  $_\n" } @done;
print "skipped ", scalar(@skip), "\n", map { "  $_\n" } @skip;
