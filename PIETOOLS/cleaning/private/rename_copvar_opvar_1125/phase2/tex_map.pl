#!/usr/bin/perl
# tex_map.pl FILE... - 11/27/2026, Phase 2 of the opvar rename: the map of Phase 1 on LaTeX notes,
# in place. Names with an underscore appear as "\_" in \code{} and as "_" in \file{}; both are
# mapped. Kept as written: the retired test files test_poscopvar_lift and test_poscopvar_direct.
# Article: "a copvar" / "a \code{copvar}" -> "an ...". Prints each changed line, old and new.
use strict; use warnings;
for my $f (@ARGV) {
  open(my $h,'<:raw',$f) or die "$f: $!"; local $/; my $s = <$h>; close $h;
  my @L = split /(?<=\n)/, $s; my $n = 0;
  for my $i (0..$#L) {
    my $o = $L[$i]; my $x = $o;
    $x =~ s/test(\\?_)poscopvar(\\?_)(lift|direct)/TEST\x01$1\x02$2$3/g;   # protect retired names
    $x =~ s/poslpivar(\\?_)sop\b/poslpivar$1legacy/g;
    $x =~ s/\bcopquadvar/opquadvar/g;
    $x =~ s/poscopvar/poslpivar/g;
    $x =~ s/\b(a)(\s+(?:\\code\{)?)(c?d?opvar)/$1n$2$3/g if 0;                # article, below
    $x =~ s/\b([Aa])(\s+)(\\code\{)?copvar/$1n$2${\($3 \/\/ '')}opvar/g;
    $x =~ s/cdopvar/dopvar/g;
    $x =~ s/copvar/opvar/g;
    $x =~ s/TEST\x01(\\?_)\x02(\\?_)(lift|direct)/test$1poscopvar$2$3/g;
    if ($x ne $o) { $n++; print "$f:", $i+1, ":\n  - $o  + $x"; $L[$i] = $x; }
  }
  open($h,'>:raw',$f) or die "$f: $!"; print $h join('',@L); close $h;
  print STDERR "$f: $n lines changed\n";
}
