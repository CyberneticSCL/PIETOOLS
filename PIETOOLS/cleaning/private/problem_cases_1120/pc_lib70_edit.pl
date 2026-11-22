# pc_lib70_edit.pl FILE
# Worktree copy of the PDE 70 example (PIETOOLS_PDE_Ex_2D_Wave_Acoustic_Boundary_Beale.m) only: replaces the state
# declarations and PDE_t by the mixed position-strain form of DIAG2M ('70mix': ph of order (2,0), e2 = ph_s2 of order
# (0,1), the acoustic coupling as e2(s1,0) = -et), the old lines kept commented as (was), and adds st.lpi.ineq.dD = 0
# after st = lpisettings_sop('light'); so libcheck runs the file with the form and the
# setting certified in the harness. Keeps the line endings of the file.
use strict; use warnings;
my ($f) = @ARGV;
open(my $in, '<', $f) or die "no $f\n"; my @L = <$in>; close $in;
my $eol = $L[0] =~ /\r\n$/ ? "\r\n" : "\n";
sub mark { my ($code, $sfx) = @_; my $m = "% MMP, 11/20/2026 (worktree test$sfx)";
           return length($code) < 76 ? sprintf("%-76s%s", $code, $m) : "$code  $m"; }
my ($i0, $i1, $ist) = (-1, -1, -1);
for my $i (0..$#L) {
    $i0  = $i if $L[$i] =~ /^ph = pde_var\(\[s1;s2\],\[0,1;0,1\],\[2;2\]\);/;
    $i1  = $i if $i0 >= 0 && $i1 < 0 && $L[$i] =~ /subs\(\[de;et\],s1,1\)==0\];/;
    $ist = $i if $L[$i] =~ /^st = lpisettings_sop\('light'\);/;
}
die "block not found\n" if $i0 < 0 || $i1 < 0;  die "st line not found\n" if $ist < 0;
my @new = ("ph = pde_var([s1;s2],[0,1;0,1],[2;0]);",
           "e2 = pde_var([s1;s2],[0,1;0,1],[0;1]);",
           "ps = pde_var([s1;s2],[0,1;0,1],[2;1]);",
           "de = pde_var(s1,[0,1],2);",
           "et = pde_var(s1,[0,1],2);",
           "PDE_t = [diff(ph,'t')==ps;",
           "         diff(ps,'t')==diff(ph,s1,2)+diff(e2,s2);",
           "         diff(e2,'t')==diff(ps,s2);",
           "         diff(de,'t')==et;",
           "         diff(et,'t')==-et-de-subs(ps,s2,0);",
           "         subs([ph;ps],s1,0)==0;     subs([ph;ps],s1,1)==0;",
           "         subs(ps,s2,1)==0;",
           "         subs(e2,s2,0)==-et;",
           "         subs([de;et],s1,0)==0;     subs([de;et],s1,1)==0];");
my @old = map { my $s = $_; $s =~ s/\r?\n$//; mark("%$s", ', was') . $eol } @L[$i0..$i1];
my @add = map { mark($_, '') . $eol } @new;
my @stl = map { mark($_, '') . $eol } ("st.lpi.ineq.dD = 0;");   # 12:24: deg_Q{2} [0 1 1] dropped (dD 0 alone certifies, DIAG2M)
splice(@L, $ist + 1, 0, @stl);                    # after the st line (below the block, so do it first)
splice(@L, $i0, $i1 - $i0 + 1, @old, @add);
open(my $w, '>', $f) or die; print $w @L; close $w;
print "edited $f: block lines ", $i0 + 1, "-", $i1 + 1, ", st line ", $ist + 1, "\n";
