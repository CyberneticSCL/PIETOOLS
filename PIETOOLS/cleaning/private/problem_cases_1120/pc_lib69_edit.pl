# pc_lib69_edit.pl FILE
# Worktree copy of the PDE 69 example (PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty.m) only: replaces the state
# declarations and PDE_t by the mixed position-strain form of DIAG69M ('69mix2': u of order (0,2), e1 = u_s1 of order
# (1,2) with e1 = 0 at s2 = 0, 1, the heat flux coupling as e1(0,s2) = -w_s1(0,s2)), the old lines kept commented as
# (was), and adds st.lpi.ineq.dD = 0 after st = lpisettings_sop('light'); so libcheck runs the file with the form and
# the setting certified in the harness. Keeps the line endings of the file.
use strict; use warnings;
my ($f) = @ARGV;
open(my $in, '<', $f) or die "no $f\n"; my @L = <$in>; close $in;
my $eol = $L[0] =~ /\r\n$/ ? "\r\n" : "\n";
sub mark { my ($code, $sfx) = @_; my $m = "% MMP, 11/20/2026 (worktree test$sfx)";
           return length($code) < 76 ? sprintf("%-76s%s", $code, $m) : "$code  $m"; }
my ($i0, $i1, $ist) = (-1, -1, -1);
for my $i (0..$#L) {
    $i0  = $i if $L[$i] =~ /^u = pde_var\(\[s1;s2\],\[0,1;0,1\],\[2;2\]\);/;
    $i1  = $i if $i0 >= 0 && $i1 < 0 && $L[$i] =~ /subs\(\[u;v;w\],s2,1\)==0\];/;
    $ist = $i if $L[$i] =~ /^st = lpisettings_sop\('light'\);/;
}
die "block not found\n" if $i0 < 0 || $i1 < 0;  die "st line not found\n" if $ist < 0;
my @new = ("u = pde_var([s1;s2],[0,1;0,1],[0;2]);",
           "e1 = pde_var([s1;s2],[0,1;0,1],[1;2]);",
           "v = pde_var([s1;s2],[0,1;0,1],[1;2]);",
           "w = pde_var([s1;s2],[0,1;0,1],[2;2]);",
           "PDE_t = [diff(u,'t')==v;",
           "         diff(v,'t')==diff(e1,s1)+diff(u,s2,2);",
           "         diff(e1,'t')==diff(v,s1);",
           "         diff(w,'t')==diff(w,s1,2)+diff(w,s2,2);",
           "         subs(u,s2,0)==0;           subs(u,s2,1)==0;",
           "         subs(e1,s1,0)==-subs(diff(w,s1),s1,0);",
           "         subs(e1,s2,0)==0;          subs(e1,s2,1)==0;",
           "         subs(v,s1,1)==0;           subs(v,s2,0)==0;     subs(v,s2,1)==0;",
           "         subs(w,s1,1)==0;           subs(w,s1,0)==subs(v,s1,0);",
           "         subs(w,s2,0)==0;           subs(w,s2,1)==0];");
my @old = map { my $s = $_; $s =~ s/\r?\n$//; mark("%$s", ', was') . $eol } @L[$i0..$i1];
my @add = map { mark($_, '') . $eol } @new;
my @stl = map { mark($_, '') . $eol } ("st.lpi.ineq.dD = 0;");
splice(@L, $ist + 1, 0, @stl);                    # after the st line (below the block, so do it first)
splice(@L, $i0, $i1 - $i0 + 1, @old, @add);
open(my $w, '>', $f) or die; print $w @L; close $w;
print "edited $f: block lines ", $i0 + 1, "-", $i1 + 1, ", st line ", $ist + 1, "\n";
