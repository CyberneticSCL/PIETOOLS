# pc_plan_t3.pl PLAN TRADE3_FINAL > OUT: a plan (phaseC_plan.csv columns) with one row per TRADE3 row whose trade
# pass cut the program: every key of its plan row (the example has one st for the primal and the dual
# test), its tier and statements, the forms of the plan row, status ready.
use strict; use warnings; my %used;
sub csv { my $l = shift; $l =~ s/\r?\n$//; return ($l =~ /("[^"]*"|[^,]*)(?:,|$)/g); }
open(my $p,'<',$ARGV[0]) or die; my $hdr = <$p>; my @P;
while (<$p>) { my @f = csv($_); pop @f if @f > 12 && $f[-1] eq ''; push @P, [@f]; }
open(my $t,'<',$ARGV[1]) or die; <$t>;
print $hdr;
while (<$t>) {
    my @f = csv($_); my ($key,$form,$tier,$stm,$n0,$m0,$n1) = @f[0..6];
    next unless $n1 ne 'NaN' && $n0 ne 'NaN' && $n1 < $n0;
    my ($row) = grep { (" ".($_->[3])." ") =~ / \Q$key\E / } @P;
    die "no plan row for $key\n" unless $row;
    my @ks = split / /, $row->[3]; if (grep { $used{$_} } @ks) { print STDERR "skip $key: its plan row is verified with another cut
"; next; } $used{$_} = 1 for @ks;
    my @o = @$row; $o[7] = $tier; $o[8] = $stm; $o[9] = $n1; $o[10] = 'ready'; $o[11] = '"TRADE3 11/20"';
    print join(',', @o[0..11]), "\n";
}
