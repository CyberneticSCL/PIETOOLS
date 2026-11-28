#!/usr/bin/perl
# MAKE_STUBS.PL ROOT - decision E1 of ../PLAN.md (maintainer, 11/28/2026: "Error, but with
# specific guidance on the correct calling format for the case which led to the specific
# error"), 11/28/2026. Writes ROOT/executives/retired_2025/<name>.m for each executive of
# PIETOOLS 2025 that became a form of an executive of this release: a call stops with the
# identifier <name>:retired and a message that gives the call of this release.
# The forms: the LPI in the 2025 headers, the map of lpiscript (2025 and now), and the
# numerical identification of 10/2026 and 11/27/2026 (../../rename_copvar_opvar_1125/phase2);
# in 2-D the default name is the coercive form and _non_coercive the T'*Q form.
use strict; use warnings; use File::Path qw(make_path);
my $root = shift or die "usage: ROOT\n";
my $dir = "$root/executives/retired_2025"; make_path($dir);
my $H2 = "The 2025 input OPTIONS.h2, a fixed bound, has no counterpart: PIETOOLS_H2_norm minimizes the bound.";
my $ND = "The executive takes a PIE in any number of spatial variables; the solver status is info.status.";
my $GN = "GAIN as in 2025: the LPI solved at gam = gain (0 or [] minimizes).";
# name, 2025 call, call now, form, notes
my @T = (
 ['PIETOOLS_PDEstability',            '[prog,P] = PIETOOLS_PDEstability(PIE,settings)',             q{[prog,P,info] = PIETOOLS_stability(PIE,settings,'P')},           'P',   ''],
 ['PIETOOLS_PDEstability_dual',       '[prog,P] = PIETOOLS_PDEstability_dual(PIE,settings)',        q{[prog,P,info] = PIETOOLS_stability(PIE,settings,'Pd')},          'Pd',  ''],
 ['PIETOOLS_PIE2PDEstability',        '[prog,P] = PIETOOLS_PIE2PDEstability(PIE,settings)',         q{[prog,P,info] = PIETOOLS_stability(PIE,settings,'Q')},           'Q',   ''],
 ['PIETOOLS_PIE2PDEstability_dual',   '[prog,P] = PIETOOLS_PIE2PDEstability_dual(PIE,settings)',    q{[prog,P,info] = PIETOOLS_stability(PIE,settings,'Qd')},          'Qd',  ''],
 ['PIETOOLS_Hinf_gain_coercive',      '[prog,P,gam] = PIETOOLS_Hinf_gain_coercive(PIE,settings)',   q{[prog,P,gam,info] = PIETOOLS_Hinf_gain(PIE,settings,'P')},       'P',   ''],
 ['PIETOOLS_Hinf_gain_dual',          '[prog,R,gam] = PIETOOLS_Hinf_gain_dual(PIE,settings)',       q{[prog,R,gam,info] = PIETOOLS_Hinf_gain(PIE,settings,'Qd')},      'Qd',  ''],
 ['PIETOOLS_Hinf_gain_dual_coercive', '[prog,P,gam] = PIETOOLS_Hinf_gain_dual_coercive(PIE,settings)', q{[prog,P,gam,info] = PIETOOLS_Hinf_gain(PIE,settings,'Pd')},   'Pd',  ''],
 ['PIETOOLS_H2_norm_o',               '[prog,W,gam,R,Q] = PIETOOLS_H2_norm_o(PIE,settings,options)', q{[prog,W,gam,R,Q,info] = PIETOOLS_H2_norm(PIE,settings,'o')},  'o',   $H2],
 ['PIETOOLS_H2_norm_c',               '[prog,W,gam,R,Q] = PIETOOLS_H2_norm_c(PIE,settings,options)', q{[prog,W,gam,R,Q,info] = PIETOOLS_H2_norm(PIE,settings,'c')},  'c',   $H2],
 ['PIETOOLS_H2_norm_o_coercive',      '[prog,Wo,gam] = PIETOOLS_H2_norm_o_coercive(PIE,settings,options)', q{[prog,Wo,gam,~,~,info] = PIETOOLS_H2_norm(PIE,settings,'oco')}, 'oco', $H2],
 ['PIETOOLS_H2_norm_c_coercive',      '[prog,Wc,gam] = PIETOOLS_H2_norm_c_coercive(PIE,settings,options)', q{[prog,Wc,gam,~,~,info] = PIETOOLS_H2_norm(PIE,settings,'cco')}, 'cco', $H2],
 ['PIETOOLS_stability_2D',            '[prog,Pop,Qop,solve_val] = PIETOOLS_stability_2D(PIE,settings)', q{[prog,P,info] = PIETOOLS_stability(PIE,settings,'P')},       'P',   $ND],
 ['PIETOOLS_stability_dual_2D',       '[prog,Pop,Qop,solve_val] = PIETOOLS_stability_dual_2D(PIE,settings)', q{[prog,P,info] = PIETOOLS_stability(PIE,settings,'Pd')}, 'Pd',  $ND],
 ['PIETOOLS_Hinf_gain_2D',            '[prog,Pop,gam,solve_val] = PIETOOLS_Hinf_gain_2D(PIE,settings,gain)', q{[prog,P,gam,info] = PIETOOLS_Hinf_gain(PIE,settings,'P',gain)}, 'P', "$ND $GN"],
 ['PIETOOLS_Hinf_gain_2D_non_coercive', '[prog,Rop,gam,solve_val] = PIETOOLS_Hinf_gain_2D_non_coercive(PIE,settings,gain)', q{[prog,R,gam,info] = PIETOOLS_Hinf_gain(PIE,settings,'Q',gain)}, 'Q', "$ND $GN"],
 ['PIETOOLS_Hinf_gain_dual_2D',       '[prog,Pop,gam,solve_val] = PIETOOLS_Hinf_gain_dual_2D(PIE,settings,gain)', q{[prog,P,gam,info] = PIETOOLS_Hinf_gain(PIE,settings,'Pd',gain)}, 'Pd', "$ND $GN"],
 ['PIETOOLS_H2_norm_2D_c',            '[prog,Wc,gam] = PIETOOLS_H2_norm_2D_c(PIE,settings,options)', q{[prog,Wc,gam,~,~,info] = PIETOOLS_H2_norm(PIE,settings,'cco')}, 'cco', "$ND $H2"],
 ['PIETOOLS_H2_norm_2D_c_non_coercive', '[prog,W,gam,R,Q] = PIETOOLS_H2_norm_2D_c_non_coercive(PIE,settings,options)', q{[prog,W,gam,R,Q,info] = PIETOOLS_H2_norm(PIE,settings,'c')}, 'c', "$ND $H2"],
 ['PIETOOLS_H2_norm_2D_o',            '[prog,Wo,gam] = PIETOOLS_H2_norm_2D_o(PIE,settings,options)', q{[prog,Wo,gam,~,~,info] = PIETOOLS_H2_norm(PIE,settings,'oco')}, 'oco', "$ND $H2"],
 ['PIETOOLS_H2_norm_2D_o_non_coercive', '[prog,W,gam,R,Q] = PIETOOLS_H2_norm_2D_o_non_coercive(PIE,settings,options)', q{[prog,W,gam,R,Q,info] = PIETOOLS_H2_norm(PIE,settings,'o')}, 'o', "$ND $H2"],
 ['PIETOOLS_Hinf_estimator_2D',       '[prog,Lop,gam,Pop,Zop] = PIETOOLS_Hinf_estimator_2D(PIE,settings,gain)', q{[prog,L,gam,P,Z,info] = PIETOOLS_Hinf_estimator(PIE,settings,gain)}, '', "$ND $GN"],
);
my $gpl = <<'EOT';
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
EOT
# wrap text into comment lines of at most 75 characters
sub wrapc { my ($text,$pre) = @_; my @w = split ' ', $text; my @out; my $l = $pre;
  for my $w (@w) { if (length($l) + length($w) + 1 > 75 && $l ne $pre) { push @out, $l; $l = $pre }
    $l .= ($l eq $pre ? '' : ' ') . $w } push @out, $l if $l ne $pre; return join("\n", @out) . "\n" }
# a MATLAB string of the message, cut into pieces joined by ...
sub mstr { my $s = shift; $s =~ s/'/''/g; my @w = split / /, $s; my @p; my $l = '';
  for my $w (@w) { if (length($l) + length($w) + 1 > 62 && $l ne '') { push @p, "$l "; $l = '' } $l .= ($l eq '' ? '' : ' ') . $w }
  push @p, $l; return join(" ...\n     ", map { "'$_'" } @p) }
my $n = 0;
for my $r (@T) { my ($name,$old,$new,$form,$note) = @$r;
  (my $exec = $new) =~ s/^.*= (\w+)\(.*$/$1/;
  my $why = $form ne '' ? "The form '$form' poses the LPI of $name." : "$exec poses the LPI of $name.";
  my $msg = "$name is an executive of PIETOOLS 2025. Call $new: "
          . ($form ne '' ? "the form '$form' poses the same LPI." : "it poses the same LPI.") . ($note ne '' ? " $note" : '');
  my $help = wrapc(uc($name) . " - an executive of PIETOOLS 2025, retired. A call stops with an error that gives the call of this release:", '% ')
           . "%\n%   PIETOOLS 2025:\n%     $old\n%   this release:\n%     $new\n%\n"
           . wrapc("$why " . ($note ne '' ? "$note " : '') . "See $exec for the forms and the outputs, and LPISETTINGS for SETTINGS.", '% ');
  my $s = "function varargout = $name(varargin) %#ok<STOUT> the call only errors\n"
        . ('%' x 75) . "\n" . $help . "%\n% PIETOOLS - $name\n%\n" . $gpl . "%\n" . ('%' x 75) . "\n"
        . "%\n% If you modify this code, document all changes carefully and include date\n"
        . "% authorship, and a brief description of modifications\n%\n"
        . "% Initial coding MMP, 11/28/2026\n\n"
        . "% The LPI of this executive of PIETOOLS 2025 is posed by an executive of\n"
        . "% this release: stop with the call that poses it.\n"
        . "error('$name:retired', ...\n    [" . mstr($msg) . "]);\nend\n";
  open(my $o, '>', "$dir/$name.m") or die "$dir/$name.m: $!"; print $o $s; close $o; $n++ }
print "$n stubs written to $dir\n";
