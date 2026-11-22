# pc_bycase.pl DIR KEY...
# Prints, for each KEY (e.g. PDE:69), every diag_amt row of results_<TAG>.csv in DIR whose key matches:
# tag, form, tier, status, measure, variables, seconds, mods. For drafting the tried / not tried sections.
use strict; use warnings;
sub csv { my $l = shift; $l =~ s/\r?\n$//; my @f = ($l =~ /("[^"]*"|[^,]*)(?:,|$)/g); s/^"|"$//g for @f; return @f; }
my ($dir, @keys) = @ARGV;
opendir(my $d, $dir) or die; my @files = sort grep { /^results_\w+\.csv$/ } readdir $d; closedir $d;
for my $k (@keys) {
    print "=== $k\n";
    for my $f (@files) {
        my ($tag) = $f =~ /^results_(\w+)\.csv$/;
        open(my $in, '<', "$dir/$f") or next; my $h = <$in>; next unless defined $h;
        my @H = csv($h); my %ix; @ix{@H} = (0..$#H);
        next unless exists $ix{key} && exists $ix{status};
        while (<$in>) {
            my @r = csv($_); next unless ($r[$ix{key}] // '') eq $k;
            my $g = sub { my $c = shift; exists $ix{$c} ? ($r[$ix{$c}] // '') : '' };
            printf "%-7s %-6s %-2s %-9s %-11s fr %-8s n %-9s %-7s %s%s\n", $tag, $g->("variant").":".$g->("param"), $g->("form"), $g->("tier"),
                $g->('status'), $g->('fr'), $g->('ndec'), $g->('sec'), $g->('mods'),
                ($g->('note') ne '' ? "  [" . substr($g->('note'), 0, 60) . "]" : '');
        }
    }
}
