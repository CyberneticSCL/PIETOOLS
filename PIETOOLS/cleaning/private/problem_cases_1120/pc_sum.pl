# pc_sum.pl RESULTS.csv... - summary of diag_amt results of the problem-case run:
# per job one line (status, feasibility ratio, numerr, variables, m, seconds,
# mods); then per key the counts by status and the job closest to fr = 1
# among the unknown ones.
use strict; use warnings;
my (%cnt, %best);
for my $f (@ARGV) {
    open(my $fh,'<',$f) or next; <$fh>;
    (my $tag = $f) =~ s/.*results_(.*)\.csv$/$1/;
    while (my $l = <$fh>) {
        $l =~ s/\r?\n$//;
        my @F = ($l =~ /("[^"]*"|[^,]*)(?:,|$)/g);
        my ($key,$form,$var,$tier,$mods,$st,$fr,$pinf,$dinf,$ne,$ndv,$m,$sec,$note,$exst,$eta) = @F; $exst //= q{}; $eta //= q{};
        $mods =~ s/"//g; $note //= ''; $note =~ s/"//g;
        printf "%-6s %-8s %-3s %-9s %-10s fr=%-8s ne=%-2s ex=%-10s eta=%-8s ndv=%-9s m=%-6s s=%-6s %s%s\n",
            $tag,$key,$form,$tier,$st,$fr,$ne,$exst,$eta,$ndv,$m,$sec,$mods,($note ne '' ? "  [".substr($note,0,60)."]" : '');
        $cnt{$key}{$st}++;
        if ($st eq 'unknown' && $fr ne 'NaN') {
            my $d = abs($fr - 1);
            $best{$key} = [$d,"$form $tier fr=$fr ($mods)"] if !$best{$key} || $d < $best{$key}[0];
        }
    }
}
print "\n";
for my $k (sort keys %cnt) {
    printf "%-8s %s%s\n", $k, join(' ', map {"$_=$cnt{$k}{$_}"} sort keys %{$cnt{$k}}),
        ($best{$k} ? "  | closest unknown: $best{$k}[1]" : '');
}
