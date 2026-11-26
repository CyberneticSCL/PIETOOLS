# H1 rescan: code (comments and continuations stripped, string literals KEPT) naming the 2025
# classes opvar/dopvar as whole words. Quote rule: ' after [\w)\]}.'] is a transpose, else a
# char literal; "" and '' inside literals are escapes; % or ... outside literals ends the code.
use strict; use warnings;
my $root = shift;
sub code_of { my $l=shift; my ($o,$i,$n,$p)=('',0,length $l,'');
  while ($i<$n) { my $c=substr($l,$i,1);
    last if $c eq '%' || substr($l,$i,3) eq '...';
    if ($c eq '"' || ($c eq "'" && $p !~ /[\w\)\]\}\.']/)) { my $q=$c; my $j=$i+1;
      while ($j<$n) { my $d=substr($l,$j,1); if ($d eq $q) { if ($j+1<$n && substr($l,$j+1,1) eq $q) { $j+=2; next } last } $j++ }
      $o.=substr($l,$i,$j-$i+1); $i=$j+1; $p=$q; next }
    $o.=$c; $p=$c; $i++ }
  return $o }
my @files = split /\n/, `git -C "$root" ls-files -- "*.m"`;
for my $f (@files) { next if $f =~ m#^(cleaning|orphaned_file)/#;
  open(my $h,'<',"$root/$f") or next; my $k=0;
  while (my $l=<$h>) { $k++; $l=~s/\r?\n$//; my $c=code_of($l);
    print "$f:$k:$l\n" if $c =~ /\b(?:opvar|dopvar)\b/ } close $h }
