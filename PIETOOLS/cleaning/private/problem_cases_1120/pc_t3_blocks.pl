# pc_t3_blocks.pl PLAN VERIFY BLOCKSDIR DDELIB NDSLIB
# For the DDE and NDS/DDF library blocks among the TRADE3 cuts that every test verified, adds after the block's
# 11/19 settings lines (worktree copy of the library only) the cut setting at the block's comment level, marked
# 11/20/2026 (blockcheck of the worktree takes lines with either marker into the block). The block is found by its
# last two code lines (BLOCKSDIR/<id>.m). Keeps the line ending of the library file.
use strict; use warnings;
sub csv { my $l = shift; $l =~ s/\r?\n$//; return ($l =~ /("[^"]*"|[^,]*)(?:,|$)/g); }
sub strip1 { my $s = shift; $s =~ s/\r?\n$//; $s =~ s/^(\s*)%/$1/; $s =~ s/^\s+|\s+$//g; return $s; }
my ($plan, $verify, $bdir, $ddelib, $ndslib) = @ARGV;
my %ok; open(my $v,'<',$verify) or die; <$v>; while (<$v>) { my @f = csv($_); $ok{$f[0]} = $f[5]; }
my %todo;   # library file -> list of [id, tier, statements]
open(my $p,'<',$plan) or die; <$p>;
while (<$p>) {
    my @f = csv($_); my ($file, $idx, $rows, $tier, $stm) = @f[1,2,3,7,8]; $stm =~ s/"//g;
    next unless $file =~ /^examples_(DDE|NDSDDF)_library/;
    my @ks = split / /, $rows; next if grep { ($ok{$_} // '') ne '1' } @ks;
    my $lib = $file =~ /DDE_library/ ? $ddelib : $ndslib;
    push @{$todo{$lib}}, [$idx, $tier, $stm];
}
sub mark { my ($code) = @_; my $m = '% MMP, 11/20/2026 (worktree test)'; return length($code) < 76 ? sprintf("%-76s%s", $code, $m) : "$code  $m"; }
for my $lib (sort keys %todo) {
    open(my $in,'<',$lib) or die; my @L = <$in>; close $in;
    my $eol = $L[0] =~ /\r\n$/ ? "\r\n" : "\n";
    for my $t (@{$todo{$lib}}) {
        my ($id, $tier, $stm) = @$t;
        open(my $b,'<',"$bdir/$id.m") or die "no block $id\n"; my @B = <$b>; close $b; shift @B;
        my @code = grep { $_ ne '' && $_ !~ /^%/ } map { my $s = $_; $s =~ s/\r?\n$//; $s =~ s/^\s+|\s+$//g; $s } @B;
        my ($c1, $c2) = @code[-2, -1];
        my @hit;
        for my $i (1..$#L) {
            next unless strip1($L[$i]) eq $c2;
            my $j = $i - 1; $j-- while $j > 0 && strip1($L[$j]) ne $c1 && ($L[$j] =~ /^\s*%\s*%/ || strip1($L[$j]) eq '');
            push @hit, $i if strip1($L[$j]) eq $c1;
        }
        die "block $id: ".scalar(@hit)." matches\n" unless @hit == 1;
        my $k = $hit[0]; $k++ while $k < $#L && $L[$k+1] =~ /MMP, 11\/(19|20)\/2026/;
        my ($lvl) = $L[$k] =~ /^(\s*%\s?)/; $lvl //= '';
        splice(@L, $k+1, 0, mark("${lvl}st = lpisettings_sop('$tier');").$eol, mark("$lvl$stm").$eol);
        print STDERR "block $id (".($lib =~ /NDSDDF/ ? 'NDSDDF' : 'DDE')."): after line ".($k+1).": $tier; $stm\n";
    }
    open(my $w,'>',$lib) or die; print $w @L; close $w;
}
