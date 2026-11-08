function test_tier1_fixes_sop()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_TIER1_FIXES_SOP() checks the substitution of values for decision
% variables in three kinds of LPI objects, against the action of the
% operators evaluated by quadrature from the 'sopvar' class definition
% (heatNd_apply), or against the value built by another route, and errors
% at the first failed check. With bb a random 'sopvar' block from
% L2^2[s] to L2^2[s], s in [0,1], of degree at most 1, bz the same block
% with its lower-integral parameter stored as the scalar 0, bf the block
% with that parameter stored as a zero matrix, gam a 'dpvar', and the
% action evaluated on f(s) = [1+s, s^2] at 7 points of [0.05,0.95]:
%   (1) subs_dvar_sop of gam*bz, bz*gam, -gam*bz and (gam*bz)' at
%       gam = 0.5 returns 'sopvar' objects acting as 0.5*bz, 0.5*bz,
%       -0.5*bz and 0.5*bf', and of (gam*bz)*bf one acting as
%       0.5*(bf*bf), each to 1e-12 relative;
%   (2) subs_dvar_sop of sopvar2sdopvar(bz,{'a'}) at a = 0 and a = 2.7
%       returns a 'sopvar' acting as bz to 1e-14 relative; for the
%       'sdopvar' gam*bb with the last column of its first B parameter
%       removed, it raises the error 'getsol_lpivar_sop:badB';
%   (3) for a 'cdopvar' Q on L2[s] from lpivar_cdopvar at random values,
%       lpigetsol_sop of the nested cell {{Q},Q} returns the 'copvar'
%       lpigetsol_sop(prog,Q) in both places.
%
% INPUT
% none.
%
% OUTPUT
% none; one printed line with the number of checks passed.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_tier1_fixes_sop
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
% Initial coding MMP, 09/29/2026
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and headers of the local functions; code
%                unchanged. The deleted header text was the history of the
%                three defects the checks target.
% MMP, 11/04/2026: Programs declared by LPIPROGRAM; LPIPROGRAM_SOP, whose
%                extra input forms (names, containers) only filled fields that
%                container routines do not read, is retired.
% MMP, 11/08/2026: The legacy half of (3) (09/29/2026) deleted: substitution
%                into an 'opvar' with 'dpvar' fields, and lpigetsol_sop and
%                LPIGETSOL on a random 'opvar' (RAND_OPVAR). The opvar
%                classes and LPIGETSOL are retired (cleaning/private), and
%                lpigetsol_sop rejects them. The local function MAXC, used
%                only there, deleted too. Help: entry (3).

% No warnings for non-canonical multipliers; check counter
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
nchk = 0;

% ------------------------------------------ a block with a scalar-0 cell
rng(1);
b0 = rand_copvar(struct('out',{{{'s'}}},'in',{{{'s'}}}),struct('out',2,'in',2),[0 1],1,0.9);
bb = b0.C{1,1};     prm = bb.params;    prm{2} = 0;
bz = sopvar(prm,bb.vars,bb.ZL,bb.ZR,bb.dom,bb.dims);
assert(isscalar(bz.params{2}) && nnz(bz.params{1}) && nnz(bz.params{3}),...
    'setup: the constructor must keep the scalar-0 shorthand');
% gam and its value; test functions, points, and the action of bz
gam = dpvar('gam');     g0 = 0.5;
f = @(X) [1+X(:,1), X(:,1).^2];     S = linspace(0.05,0.95,7)';
yb = heatNd_apply(bz,f,S);

% (1) scalar forms: subs(form, gam = g0) = c*bz by action
forms = {@() gam*bz, g0; @() bz*gam, g0; @() -gam*bz, -g0; @() (gam*bz)', NaN};
for k = 1:size(forms,1)
    % Substitute; the result must be fixed; compare the actions
    X = forms{k,1}();
    Y = subs_dvar_sop(X,{'gam'},g0);
    assert(isa(Y,'sopvar'),'(1) form %d: %s after substitution',k,class(Y));
    if isnan(forms{k,2})    % adjoint: compare with g0*bf' applied, bf = bz
        % with the zero cell stored at full size (@sopvar/ctranspose errors
        % on the scalar-0 shorthand itself, a pre-existing defect)
        prf = bb.params;    prf{2} = sparse(size(prf{2},1),size(prf{2},2));
        bf0 = sopvar(prf,bb.vars,bb.ZL,bb.ZR,bb.dom,bb.dims);
        yr = g0*heatNd_apply(bf0',f,S);
    else
        yr = forms{k,2}*yb;
    end
    e = norm(heatNd_apply(Y,f,S)-yr)/norm(yr);
    assert(e<1e-12,'(1) form %d: relative error %.2e',k,e);
    nchk = nchk+1;
end
% composition with the scaled block: subs((gam*bz)*bf) = g0*(bf*bf), bf
% the full-size copy of bz (the class compositions error on a shorthand
% OPERAND, pre-existing; here only the scaled factor is under test)
prf = bb.params;    prf{2} = sparse(size(prf{2},1),size(prf{2},2));
bf = sopvar(prf,bb.vars,bb.ZL,bb.ZR,bb.dom,bb.dims);
assert(norm(heatNd_apply(bf,f,S)-yb)<1e-14*norm(yb),'setup: bf is not bz');
Y = subs_dvar_sop((gam*bz)*bf,{'gam'},g0);
yr = g0*heatNd_apply(bf*bf,f,S);
e = norm(heatNd_apply(Y,f,S)-yr)/norm(yr);
assert(e<1e-12,'(1) (gam*bz)*bz: relative error %.2e',e);
nchk = nchk+1;

% (2) sopvar2sdopvar of the shorthand block, substituted, is the block
for a0 = [0, 2.7]
    Y = subs_dvar_sop(sopvar2sdopvar(bz,{'a'}),{'a'},a0);
    e = norm(heatNd_apply(Y,f,S)-yb)/norm(yb);
    assert(isa(Y,'sopvar') && e<1e-14,'(2) a = %g: relative error %.2e',a0,e);
    nchk = nchk+1;
end
% a nonzero mis-sized cell is still refused (not silently zeroed)
Xb = gam*bb;        prmb = Xb.params;   % B = b*vec(C)', nonzero
assert(nnz(prmb.B{1})>0,'setup: the control needs a nonzero B');
prmb.B{1} = prmb.B{1}(:,1:end-1);
Xb.params = prmb;
e = get_err(@() subs_dvar_sop(Xb,{'gam'},1));
assert(strcmp(e.identifier,'getsol_lpivar_sop:badB'),'(2) mis-sized nonzero B: %s',e.identifier);
nchk = nchk+1;

% --------------------------------------- (3) nested container cells        % MMP, 11/08/2026
% Deleted: the legacy half of (3) (MMP, 09/29/2026), substitution into an   % MMP, 11/08/2026
% 'opvar' with 'dpvar' fields and LPIGETSOL on a random 'opvar'; the opvar  % MMP, 11/08/2026
% classes and LPIGETSOL are retired.                                        % MMP, 11/08/2026
% nested cells of container decision operators
% prog = lpiprogram_sop({'s'},[0 1]);                                       % MMP, 11/04/2026 (was)
prog = lpiprogram(polynomial({'s'}),[0 1]);                                 % MMP, 11/04/2026
[prog,Q] = lpivar_cdopvar(prog,1,{{'s'}},[0 1],1);
prog.solinfo = struct('RRx',randn(numel(prog.decvartable),1),'x',[],'info',1);
Qs = lpigetsol_sop(prog,Q);
r = lpigetsol_sop(prog,{{Q},Q});
assert(isa(r{1}{1},'copvar') && isequal(r{1}{1},Qs) && isequal(r{2},Qs),'(3) nested cell {{Q},Q}');
nchk = nchk+1;

% Summary
fprintf('test_tier1_fixes_sop passed (%d checks).\n',nchk);
end


function e = get_err(f)
% E = GET_ERR(F) returns the error that F() raises, or a 'struct' with an
% empty identifier and message when it raises none.
% INPUT  f: 'function_handle', called with no argument.
% OUTPUT e: 'MException' object, or 'struct' with fields identifier and
%        message ('char').
e = struct('identifier','','message','');
try,        f();
catch ME,   e = ME;
end
end

% The local function MAXC (MMP, 09/29/2026), used only by the deleted       % MMP, 11/08/2026
% legacy half of (3), is deleted.                                           % MMP, 11/08/2026
