% Test_polyopvar_sosquadvar.m -- check polyopvar_sosquadvar structurally:
%
% PART 1 - 'free' option: block dimensions must match the sizes of the
%   input bases exactly, and all decision variables across the d^2 blocks
%   must be distinct (no accidental aliasing between unrelated blocks).
%
% PART 2 - 'pos' option, equal-size blocks: the declared Gram matrix must
%   be exactly symmetric, i.e. Pcell{i,i}-Pcell{i,i}.' must be identically
%   the zero dpvar.
%
% PART 3 - option handling:
%   3a: 'sym' with mismatched per-block dimensions (equal totals,
%       mismatched blocks) must error, as sosquadvar's own check only
%       compares the block TOTALS sum(m)==sum(n), not each block, since
%       Z1,Z2 are passed as a single cell {1};
%   3b: an option string in a different case (e.g. 'POS') must succeed,
%       as sosquadvar matches the option case-sensitively.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_polyopvar_sosquadvar
%
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ (CC), 10/15/2026: Initial coding;

clear; clear stateNameGenerator; close all; clc;

pvar s s_dum
dom = [0,1];
opdeg = 1;
Zmon = monomials([s,s_dum],0:opdeg);
Zop  = opvar();
Zop.var1 = s;   Zop.var2 = s_dum;   Zop.I = dom;
Zop.R.R0 = [0*Zmon;0*Zmon];
Zop.R.R1 = [Zmon;0*Zmon];
Zop.R.R2 = [0*Zmon;Zmon];

Z   = tensopmat(sopvar2tensopvar(opvar2sopvar(Zop)));
x   = polyopvar('x1',1,s,dom);
Zx  = Z*x;              % degree-1 basis
Zs2 = DMB(Zx,Zx);        % degree-2 (Def.4) basis, different size than Zx
n1 = size(Zx,1);
n2 = size(Zs2,1);

%% Part 1: 'free' option, two unequal-size blocks.
prog = piesos_program(x);
[prog,Pcell_free] = polyopvar_sosquadvar(prog,{Zx,Zs2},{Zx,Zs2},'free');

pass1 = isequal(size(Pcell_free{1,1}),[n1,n1]) && isequal(size(Pcell_free{1,2}),[n1,n2]) ...
     && isequal(size(Pcell_free{2,1}),[n2,n1]) && isequal(size(Pcell_free{2,2}),[n2,n2]);

dv = cell(4,1);
dv{1} = Pcell_free{1,1}.dvarname(:);    dv{2} = Pcell_free{1,2}.dvarname(:);
dv{3} = Pcell_free{2,1}.dvarname(:);    dv{4} = Pcell_free{2,2}.dvarname(:);
dv_all = unique(vertcat(dv{:}));
total_dvars = n1^2 + n1*n2 + n2*n1 + n2^2;
pass1 = pass1 && size(dv_all,1)==total_dvars;
fprintf(" --- Part 1 ('free' blocks sized correctly, no dvar aliasing): %s ---\n",pass_str(pass1));

%% Part 2: 'pos' option, equal-size block - must be exactly symmetric.
prog2 = piesos_program(x);
[prog2,Pcell_pos] = polyopvar_sosquadvar(prog2,{Zx},{Zx},'pos');
D = Pcell_pos{1,1} - Pcell_pos{1,1}.';
pass2 = ~any(any(full(D.C)));
fprintf(" --- Part 2 ('pos' block is exactly symmetric): %s ---\n",pass_str(pass2));

%% Part 3a: 'sym' with mismatched per-block dims must error.
prog3 = piesos_program(x);
try
    polyopvar_sosquadvar(prog3,{Zx,Zs2},{Zs2,Zx},'sym');
    pass3a = false;
catch
    pass3a = true;
end
fprintf(" --- Part 3a ('sym' + mismatched block dims now errors): %s ---\n",pass_str(pass3a));

%% Part 3b: option case-insensitivity must be honored end-to-end.
prog4 = piesos_program(x);
try
    polyopvar_sosquadvar(prog4,{Zx},{Zx},'POS');
    pass3b = true;
catch
    pass3b = false;
end
fprintf(" --- Part 3b ('POS', uppercase, now succeeds): %s ---\n",pass_str(pass3b));

all_pass = pass1 && pass2 && pass3a && pass3b;
fprintf("\n=== Test_polyopvar_sosquadvar: %s ===\n", pass_str(all_pass));


function str = pass_str(tf)
if tf
    str = 'PASS';
else
    str = 'FAIL';
end
end
