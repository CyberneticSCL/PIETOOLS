# Lines naming each candidate token (whole word; '_' counts as a word character, so
# settings_sopvar does not match settings_sopvar_light), by area.
use strict; use warnings;
my ($namesf, $root) = (shift, shift);
open my $n, '<', $namesf or die; my @names = map { s/\r?\n//r } <$n>; close $n;
my $re = join '|', map { quotemeta } sort { length $b <=> length $a } @names;
my (%lines, %files, %area_lines);
while (my $f = <STDIN>) { $f =~ s/\r?\n//; my $p = "$root/$f"; open my $fh, '<', $p or next;
  (my $rel = $f) =~ s#^PIETOOLS/##;
  my $area = $rel =~ m{^cleaning/} ? 'cleaning' : $rel =~ m{^orphaned_file/} ? 'orphaned' : $rel =~ m{^testing_codes/} ? 'tests' : $rel =~ /\.(md|tex|txt|bib)$/ ? 'toolbox-doc' : 'toolbox';
  my %seen;
  while (my $l = <$fh>) { my %hit; while ($l =~ /(?<![A-Za-z0-9_])($re)(?![A-Za-z0-9_])/g) { $hit{$1}++ }
    next unless %hit; $area_lines{$area}++; for my $h (keys %hit) { $lines{$h}{$area}++; $seen{$h}++ } }
  close $fh; for my $h (keys %seen) { $files{$h}{$area}++ } }
my @areas = qw(toolbox toolbox-doc tests cleaning orphaned);
printf "%-28s %s\n", 'name', join(' ', map { sprintf '%12s', $_ } @areas);
for my $h (@names) { printf "%-28s %s\n", $h, join(' ', map { sprintf '%5d/%-6d', $files{$h}{$_}//0, $lines{$h}{$_}//0 } @areas) }
printf "%-28s %s\n", 'LINES (any name)', join(' ', map { sprintf '%12d', $area_lines{$_}//0 } @areas);
