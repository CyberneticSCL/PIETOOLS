# As d1_release.pl, but compares lines with per-line markers and surrounding blanks removed,
# and splits the count into code lines and comment-only lines.
use strict; use warnings;
my ($R,$tsv,$tab,$show) = @ARGV; $show //= '^$';
my %idx; for (`git -C "$R" ls-tree -r --name-only origin/PIETOOLS-2025`) { chomp; my $b=$_; $b=~s#.*/##; push @{$idx{$b}},$_; }
my @map; open(my $T,'<',$tab) or die; <$T>;
while(<$T>){ chomp; my @c=split /\t/; next if $c[5]=~/retire/; push @map,[$c[0],$c[5]] if $c[0] ne 'copvar' && $c[0] ne 'cdopvar'; }
@map = sort { length($b->[0]) <=> length($a->[0]) } @map;
sub mapline { my $s=shift; for my $m (@map){ my($o,$n)=@$m; my $O=uc $o; my $N=uc $n; $s=~s/\b$o\b/$n/g; $s=~s/\b$O\b/$N/g; }
  $s=~s/\ba (['`]?)copvar/an $1opvar/g; $s=~s/\bA (['`]?)copvar/An $1opvar/g;
  $s=~s/cdopvar/dopvar/g; $s=~s/CDOPVAR/DOPVAR/g; $s=~s/Cdopvar/Dopvar/g;
  $s=~s/copvar/opvar/g; $s=~s/COPVAR/OPVAR/g; $s=~s/Copvar/Opvar/g; return $s; }
sub norm { my $s=shift; $s=~s/\s+%\s*[A-Z][A-Za-z ,()]*\d{1,2}\/\d{1,2}\/\d{2,4}.*$//; $s=~s/^\s+|\s+$//g; return $s; }
open(my $F,'<',$tsv) or die; my ($tc,$tm,$rc,$rm)=(0,0,0,0);
while(<$F>){ chomp; my ($k,$f)=split /\t/; next unless $k eq 'existing';
  my $b=$f; $b=~s#.*/##; my %old;
  for my $p (@{$idx{$b}||[]}){ for (`git -C "$R" show "origin/PIETOOLS-2025:$p"`){ s/\r?\n$//; $old{norm($_)}=1; } }
  for (`git -C "$R" show "HEAD:$f"`){ s/\r?\n$//; next unless /copvar|cdopvar/i;
    my $com = /^\s*%/ ? 1 : 0; my $m = norm(mapline($_));
    if ($com) { $tm++ } else { $tc++ }
    if ($old{$m}) { if ($com) { $rm++ } else { $rc++ } }
    elsif ($f =~ /$show/) { printf "  %s  %s\n", ($com?'comment':'code   '), substr($m,0,140); } } }
print "code lines $tc (identical to a 2025 line after the map: $rc); comment-only lines $tm (identical: $rm)\n";
