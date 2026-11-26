# For each name: files whose CODE (comments and block comments removed, strings kept) uses it as a word,
# split into release (toolbox) and testing_codes; the defining file (basename == name) excluded.
use strict; use warnings;
my @names = split ' ', shift @ARGV; my %rel; my %tst; my %relf; my %tstf;
for my $f (@ARGV) {
  open my $fh, '<', $f or next; my $blk = 0; my $code = '';
  while (my $l = <$fh>) { $l =~ s/\r?\n$//;
    if ($l =~ /^\s*%\{\s*$/) { $blk++; next } if ($blk && $l =~ /^\s*%\}\s*$/) { $blk--; next } next if $blk;
    my $inq=0; my $prev=''; my $cut=length $l;
    for (my $i=0;$i<length $l;$i++){ my $c=substr($l,$i,1);
      if ($c eq "'") { if ($inq) {$inq=0} elsif ($prev !~ /[\w\)\]\}\.'"]/) {$inq=1} }
      elsif ($c eq '"' && !$inq) { my $j=index($l,'"',$i+1); $i=$j if $j>0 }
      elsif ($c eq '%' && !$inq) { $cut=$i; last }
      $prev=$c if $c ne ' '; }
    $code .= substr($l,0,$cut) . "\n"; }
  close $fh; (my $base = $f) =~ s#.*/##; $base =~ s/\.m$//; my $istest = $f =~ m{/testing_codes/};
  for my $n (@names) { next if $n eq $base; my $k = () = $code =~ /(?<![A-Za-z0-9_])\Q$n\E(?![A-Za-z0-9_])/g; next unless $k;
    if ($istest) { $tst{$n} += $k; $tstf{$n}++ } else { $rel{$n} += $k; $relf{$n}++ } }
}
printf "%-24s %14s %14s\n", 'name', 'release f/uses', 'tests f/uses';
for my $n (@names) { printf "%-24s %8d/%-5d %8d/%-5d\n", $n, $relf{$n}//0, $rel{$n}//0, $tstf{$n}//0, $tst{$n}//0 }
