# pc_lib44x3_edit.pl VARIANT FILE   (VARIANT 44 or X3)
# Worktree copy of the PDE 44 example (PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay.m) or of X3
# (PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D.m) only: replaces the five states x1-x5 and PDE_t by the minimal form
# of DIAG44D/DIAG44E (u of order 2, v = u_t of order 1, the delayed velocity z; on s, or z on s2 for X3), the old
# lines kept commented as (was); runs the dual test (stability_dual = 1 with stability_form 'P', i.e. 'Pd', in place of
# stability = 1 and 'Q'); adds st.epneg after st = lpisettings_sop('light'); (44: 1e-1, X3: 1e-2, the settings
# certified in DIAG44F). Keeps the line endings of the file.
use strict; use warnings;
my ($var, $f) = @ARGV;  die "variant 44 or X3\n" unless $var =~ /^(44|X3)$/;
open(my $in, '<', $f) or die "no $f\n"; my @L = <$in>; close $in;
my $eol = $L[0] =~ /\r\n$/ ? "\r\n" : "\n";
sub mark { my ($code, $sfx) = @_; my $m = "% MMP, 11/20/2026 (worktree test$sfx)";
           return length($code) < 76 ? sprintf("%-76s%s", $code, $m) : "$code  $m"; }
my ($i0, $i1, $ist, $ifl, $ifm) = (-1, -1, -1, -1, -1);
for my $i (0..$#L) {
    $i0  = $i if $L[$i] =~ /^x1 = pde_var\(\);/;
    $i1  = $i if $i0 >= 0 && $i1 < 0 && $L[$i] =~ /subs\(x5,s2?,0\)==subs\(x3,s1?,(1|b)\)\];/;
    $ist = $i if $L[$i] =~ /^st = lpisettings_sop\('light'\);/;
    $ifl = $i if $L[$i] =~ /^evalin\('base','stability = 1;'\);/;
    $ifm = $i if $L[$i] =~ /^evalin\('base','stability_form = ''Q'';'\);/;
}
die "block not found\n" if $i0 < 0 || $i1 < 0;  die "st/flag/form line not found\n" if $ist < 0 || $ifl < 0 || $ifm < 0;
my @new = $var eq '44' ?
  ("u = pde_var(s,[0,1],2);", "v = pde_var(s,[0,1],1);", "z = pde_var(s,[0,1],1);",
   "PDE_t = [diff(u,'t')==v;",
   "         diff(v,'t')==diff(u,s,2);",
   "         diff(z,'t')==-(1/tau)*diff(z,s);",
   "         subs(u,s,0)==0;          subs(v,s,0)==0;",
   "         subs(diff(u,s),s,1)==-k*(1-mu)*subs(v,s,1)-k*mu*subs(z,s,1);",
   "         subs(z,s,0)==subs(v,s,1)];") :
  ("u = pde_var(s1,[a,b],2);", "v = pde_var(s1,[a,b],1);", "z = pde_var(s2,[0,1],1);",
   "PDE_t = [diff(u,'t')==v;",
   "         diff(v,'t')==diff(u,s1,2);",
   "         diff(z,'t')==-(1/tau)*diff(z,s2);",
   "         subs(u,s1,a)==0;         subs(v,s1,a)==0;",
   "         subs(diff(u,s1),s1,b)==-k*(1-mu)*subs(v,s1,b)-k*mu*subs(z,s2,1);",
   "         subs(z,s2,0)==subs(v,s1,b)];");
my $ep = $var eq '44' ? 'st.epneg = 1e-1;' : 'st.epneg = 1e-2;';
my $was = sub { my $s = shift; $s =~ s/\r?\n$//; return mark("%$s", ', was') . $eol; };
# bottom-up so the indices stay valid
splice(@L, $ist + 1, 0, mark($ep, '') . $eol);
splice(@L, $i0, $i1 - $i0 + 1, (map { $was->($_) } @L[$i0..$i1]), (map { mark($_, '') . $eol } @new));
splice(@L, $ifm, 1, $was->($L[$ifm]), mark("evalin('base','stability_form = ''P'';');", '') . $eol);
splice(@L, $ifl, 1, $was->($L[$ifl]), mark("evalin('base','stability_dual = 1;');", '') . $eol);
open(my $w, '>', $f) or die; print $w @L; close $w;
print "edited $f ($var): flag ", $ifl + 1, ", form ", $ifm + 1, ", block ", $i0 + 1, "-", $i1 + 1, ", st ", $ist + 1, "\n";
