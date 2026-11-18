#!/usr/bin/perl
# summ_qd.pl LOG...: a markdown table of QD (direct test) lines: system,
# parameter, term set (the psatz_rows of the terms JSON, abbreviated), Q
# length, dw, status, ndv, t. SCRATCH, MMP, 11/16/2026.
use strict; use warnings;
my %abbr = (
  'none' => 'none', 'auto' => 'auto',
  '[[1,0],[2,0],[0,1],[0,2]]|[1,1,1,1]|false' => 'facesS',
  '[[1,0],[2,0],[0,1],[0,2]]|[0,0,0,0]|true'  => 'faces0',
  '[[1,0],[2,0],[0,1],[0,2]]|[1,1,1,1]|true'  => 'faces',
  '[[3,0],[0,3]]|[1,1]|true' => 'box',
  '[[1,0],[2,0],[0,1],[0,2],[1,1],[1,2],[2,1],[2,2]]|[1,1,1,1,1,1,1,1]|true' => 'corners',
  '[[1,0],[2,0]]|[1,1]|true' => 'faces_s1',
  '[[0,1],[0,2]]|[1,1]|true' => 'faces_s2',
  '[3,0]|1|true' => 'box_s1',
  '[0,3]|1|true' => 'box_s2',
  '[[3,0],[1,0],[2,0],[0,1],[0,2]]|[1,1,1,1,1]|true' => 'box_s1+faces',
);
print "| system | parameter | terms | Q length | dw | status | ndv | t (s) | eta |\n|---|---|---|---|---|---|---|---|---|\n";
for my $log (@ARGV) {
    open(my $fh, '<', $log) or next;
    while (my $l = <$fh>) {
        next unless $l =~ /^\@\@ QD (\S+) par=(\S+) qdeg=(\S+) terms=(\S+) dw=(\d+) (\S+) ndv=(\S+) t=(\S+) .*eta=(\S+)/;
        my ($s, $p, $q, $t, $dw, $st, $nd, $tt, $eta) = ($1, $2, $3, $4, $5, $6, $7, $8, $9);
        my $key;
        if ($t =~ /"psatz":"none"/) { $key = 'none' }
        elsif ($t eq '{}') { $key = 'auto' }
        elsif ($t =~ /"psatz_rows":(\[.*?\]),"psatz_rows_offset":(\[.*?\]|\d+),"psatz_rows_perdir":(true|false)/) { $key = "$1|$2|$3" }
        my $nm = $abbr{$key // ''} // ($key // $t);
        my ($ql) = $q =~ /"mult":(\d+)/;  $ql = ($ql // 0) + 1;
        $s =~ s/_.*//;
        print "| $s | $p | $nm | $ql | $dw | $st | $nd | $tt | $eta |\n";
    }
    close $fh;
}
