function Test_SOS_DP_v2()
% TEST_SOS_DP_V2 verifies SOS_DP_v2(prog,d,opdeg,x,dom,opts), which
% declares a distributed SOS polynomial
% DP = <Zop*Z_d(x),P*Zop*Z_d(x)>, P>=0, with Z_d(x)=[x;...;x^d], using
% piesos_sosvar.
%
% (1) Structure: DP is a scalar 'polyopvar' with monomials of degrees
%     2,...,2d, Pcell is d x d, and the program declares an SOS variable.
%     The options are passed on (errors for invalid psatz).
% (2) Semantics, by solving small feasibility programs DP == target, with
%     targets that are sums of squares, or are not positive, by
%     construction:
%     - target = ||T*x||^2 + ||R*x||^2 (T the map to u, R=du/ds: kernel
%       only, no delta term): feasible, and -target infeasible;
%     - target = <x,x> has a delta (multiplier) term: infeasible for the
%       default option exclude=[1;0;0] (no multiplier in the basis
%       operators, as in SOS_DP), feasible for opts.exclude=[0;0;0];
%     - target = int s(1-s)x(s)^2 ds: not SOS (the weight is negative
%       outside [0,1]), infeasible with psatz=0, feasible for psatz=2,
%       the multipliers (s-a)^k*(b-s)^(2-k) including s(1-s);
%     - degree 4 (d=2): target = <x,x> + 3||(T*x).*(R*x)||^2, a sum of
%       squares of the basis [x; (Tx)(Rx)], feasible for d=2, infeasible
%       for d=1 (no degree-4 terms in DP), and the same target with a
%       negative sign for the quartic term infeasible.
% The feasibility flags are those of mosek (pinf, dinf, numerr).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_SOS_DP_v2
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

clear stateNameGenerator
pvar s t
dom = [0,1];
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2); subs(x,s,dom(1))==0; subs(x,s,dom(2))==0];
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
Top = PIE.T;    xs = PIE.x;
Top_o = sopvar2opvar(Top.ops{1}.ops{1,1});
Rop = tensopmat(sopvar2tensopvar(opvar2sopvar(diff(Top_o,Top_o.var1,1,'pure'))));
u = Top*xs;     us = Rop*xs;

%% (1) Structure
for d=1:3
    prog = piesos_program(xs);
    nvar0 = numel(prog.var.type);
    [prog,DP,Pcell,Zop] = SOS_DP_v2(prog,d,1,xs,dom);
    assert(isa(DP,'polyopvar') && isequal(size(DP),[1,1]),'(1a) DP must be a scalar polyopvar (d=%d).',d);
    assert(isequal(unique(sum(DP.degmat,2)),(2:2*d)'),'(1b) DP must have the degrees 2,...,2d (d=%d).',d);
    assert(isequal(size(Pcell),[d,d]),'(1c) Pcell must be d x d (d=%d).',d);
    assert(isa(Zop,'tensopmat'),'(1d) Zop must be a tensopmat.');
    assert(numel(prog.var.type)>nvar0 && any(strcmp(prog.var.type(nvar0+1:end),'sos')),'(1e) SOS_DP must declare an SOS variable.');
end
% options
threw = false;
try,    SOS_DP_v2(piesos_program(xs),1,1,xs,dom,struct('psatz','a'));  catch,  threw = true;   end
assert(threw,'(1f) an error is expected for an invalid psatz.');
% psatz with order>0 declares more SOS blocks than psatz=0
progA = piesos_program(xs);     progA = SOS_DP_v2(progA,1,1,xs,dom);
progB = piesos_program(xs);     progB = SOS_DP_v2(progB,1,1,xs,dom,struct('psatz',2));
nsosA = sum(strcmp(progA.var.type,'sos'));      nsosB = sum(strcmp(progB.var.type,'sos'));
assert(nsosB>nsosA,'(1g) psatz=2 must declare more SOS variables (%d vs %d).',nsosB,nsosA);
% the options of the basis are passed on: degs and exclude
[~,~,~,Zop1] = SOS_DP_v2(piesos_program(xs),1,0,xs,dom);
[~,~,~,Zop2] = SOS_DP_v2(piesos_program(xs),1,2,xs,dom);
[~,~,~,Zop3] = SOS_DP_v2(piesos_program(xs),1,2,xs,dom,struct('exclude',[0;0;0]));
assert(size(Zop2,1)>size(Zop1,1),'(1h) a larger degree opdeg must give more basis operators.');
assert(size(Zop3,1)==size(Zop2,1)+1,'(1i) exclude=[0;0;0] must add the multiplier term (1 row): %d vs %d.',size(Zop3,1),size(Zop2,1));

%% (2) Semantics
target_k = innerprod(u,u) + innerprod(us,us);           % kernel terms only
target_d = innerprod(xs,xs);                            % delta term
target_w = innerprod(xs,s*(1-s),xs);                    % weighted delta term
w = times(u,us);
target_4 = innerprod(xs,xs) + 3*innerprod(w,w);         % degree 2 and 4
target_4n = innerprod(xs,xs) - 3*innerprod(w,w);

cases = {
 % name,                  d, opdeg, opts,                         target,     expected feasible
 'kernel terms',          1, 1, struct(),                         target_k,   true;
 '-kernel terms',         1, 1, struct(),                         -target_k,  false;
 'delta, default excl.',  1, 1, struct(),                         target_d,   false;
 'delta, no exclude',     1, 1, struct('exclude',[0;0;0]),        target_d,   true;
 'weighted, psatz=0',     1, 1, struct('exclude',[0;0;0]),        target_w,   false;
 'weighted, psatz=2',     1, 1, struct('exclude',[0;0;0],'psatz',2), target_w, true;
 'degree 4, d=2',         2, 2, struct('exclude',[0;0;0]),        target_4,   true;
 'degree 4, d=1',         1, 2, struct('exclude',[0;0;0]),        target_4,   false;
 'degree 4, wrong sign',  2, 2, struct('exclude',[0;0;0]),        target_4n,  false;
};
nfail = 0;
for k=1:size(cases,1)
    prog = piesos_program(xs);
    [prog,DP] = SOS_DP_v2(prog,cases{k,2},cases{k,3},xs,dom,cases{k,4});
    prog = piesos_eq(prog,DP-cases{k,5});
    info = solve_info(prog);
    feasible = ~info.pinf && ~info.dinf && ~info.numerr;
    status = 'ok';
    if feasible~=cases{k,6}
        status = 'FAIL';    nfail = nfail+1;
    end
    fprintf('%-5s %-22s expected feasible=%d: pinf=%d dinf=%d numerr=%d feasratio=%.3f\n', ...
        status,cases{k,1},cases{k,6},info.pinf,info.dinf,info.numerr,info.feasratio);
end
assert(nfail==0,'%d of the SOS_DP feasibility cases failed.',nfail);

disp('Test_SOS_DP_v2: all tests passed.')
end


function info = solve_info(prog)
% Solve the program with mosek, without output, and return the solver info.
opts.solver = 'mosek';
opts.simplify = true;
evalc('prog_sol = piesos_solve(prog,opts);');
info = prog_sol.solinfo.info;
end
