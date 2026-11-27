#!/usr/bin/perl
# mkq.pl Q_NEW.txt > q_old.tsv - 11/27/2026: the queue of the export before the rename from that
# of the export after it: the class tests take their old names, P0_ROOT names the old export, and
# the D6 test (new in Phase 1b) is dropped. Row ids stay the same.
use strict; use warnings;
while (my $l = <>) {
  $l =~ s/\r?\n$//;
  my @f = split /\|/, $l, -1;
  next if $f[1] eq 'lps_test_call2025_sop';
  $f[3] =~ s/^Test_opvar$/Test_copvar/;
  $f[3] =~ s/^Test_dopvar$/Test_cdopvar/;
  $f[3] =~ s/^Test_qp_opvar_nd$/Test_qp_copvar_nd/;
  $f[3] =~ s/wt1127\\new\\/wt1127\\old\\/;
  print join("\t", @f), "\n";
}
