# Code (comments and strings excluded) in live .m files that uses a NEW name of the table as a
# word: such a call resolves nowhere today and binds to the renamed function after Phase 1.
use strict; use warnings;
my ($tokf, $root) = (shift, shift);
open my $t, '<', $tokf or die; my @tok = map { s/\r?\n//r } <$t>; close $t;
my $re = join '|', map { quotemeta } sort { length $b <=> length $a } @tok;
sub code_of { my $l = shift; my ($o,$i,$n,$p) = ('',0,length $l,'');
  while ($i < $n) { my $c = substr($l,$i,1);
    return $o if $c eq '%' || substr($l,$i,3) eq '...';
    if ($c eq '"' || ($c eq "'" && $p !~ /[\w\)\]\}\.']/)) { my $q=$c; my $j=$i+1;
      while ($j<$n) { my $d=substr($l,$j,1); if ($d eq $q) { if ($j+1<$n && substr($l,$j+1,1) eq $q) { $j+=2; next } last } $j++ }
      $o .= ' ' x ($j-$i+1); $i=$j+1; $p=$q; next }
    $o .= $c; $p = ($c eq ' ' || $c eq "\t") ? ' ' : $c; $i++ }
  return $o }
while (my $f = <STDIN>) { $f =~ s/\r?\n//; open my $h, '<', "$root/$f" or next; my $blk = 0; my $ln = 0;
  (my $base = $f) =~ s#.*/##; $base =~ s/\.m$//;
  while (my $l = <$h>) { $ln++; $l =~ s/\r?\n$//;
    if ($l =~ /^\s*%\{\s*$/) { $blk++; next } if ($blk && $l =~ /^\s*%\}\s*$/) { $blk--; next } next if $blk;
    my $c = code_of($l);
    while ($c =~ /(?<![A-Za-z0-9_.])($re)(?![A-Za-z0-9_])/g) { print "$f:$ln\t$1\t", ($l =~ s/^\s+//r), "\n" } }
  close $h }
