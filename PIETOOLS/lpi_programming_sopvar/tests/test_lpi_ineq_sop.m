function test_lpi_ineq_sop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_LPI_INEQ_SOP tests LPI_INEQ_SOP, solving with MOSEK when installed
% and SeDuMi otherwise, and errors at the first failed check, with rel_b =
% ||At'*x - b||/||b|| the relative equality residual of a solved program:
%   (1) the fixed 'copvar' Pt on L2[0,1] with kernel 2 - max(s,t), equal to
%       V'*V + L'*L for (V x) = int_0^1 x and (L x)(s) = int_0^s x, is
%       certified with the default degrees and terms: rel_b <= 1e-6;
%   (2) for the PIE T v_t = A0 v of u_t = u_s1s1 on [0,1], u = 0 at s1 = 0
%       and s1 = 1, and a free 'cdopvar' P of degree 1, the LPI
%       P'*T = T'*P, P'*T + T'*P - 2*ep^2*T'*T >= 0 (ep = 0.1) and
%       -(P'*A0 + A0'*P + kappa*(P'*T + T'*P)) >= 0, both inequalities
%       imposed by LPI_INEQ_SOP, is feasible at kappa = 9 (rel_b <= 1e-6
%       and pinf = 0) and rejected at kappa = 11 > pi^2 (pinf ~= 0 or
%       rel_b > 1e-6);
%   (3) for the PIE of u_t = u_s1s1 + u_s2s2 on [0,1]^2 (u = 0 at s1 = 0,
%       s1 = 1 and s2 = 0, u_s2 = 0 at s2 = 1), T'*T >= 0 is certified
%       with the default terms: rel_b <= 1e-6;
%   (4) for Pt, the options psatz 'none', psatz 'faces', a given deg,
%       prune false and dD = dw = 0 build (no solve), with fewer decision
%       variables than the default for 'none' and dD = dw = 0, more for
%       'faces', at least as many for prune false, and info.degrees.given
%       set for the given deg;
%   (5) for the 'dopvar' of POSLPIVAR on a 1-D program, LPI_INEQ_SOP adds
%       the decision variables and expressions LPI_INEQ adds;
%   (6) opts.like = Pt'*Pt gives the lift and weight that GET_LIFT_DEGS
%       reads from Pt'*Pt, with more decision variables than the default,
%       and a 'double' opts.like is refused;
%   (7) for X1 = -(T'*P1*A0 + A0'*P1*T) and X1 - T'*P2*T on the PIE of
%       (2), P1 and P2 of POSCOPVAR at degrees 2 and 3: info.ndec equals
%       the growth of prog.decvartable in the five calls of (7), and
%       sum(N.*(N+1)/2) the decision variables of Pop in (1), (3), the
%       given deg of (4), (6) and those calls; X1 is sized from P1, with
%       the degrees and decision variables of opts.like = P1, and from X1
%       itself, with more decision variables, without prog.sop_pos;
%       X1 - T'*P2*T is sized from P1 and P2, with the degrees and decision
%       variables of opts.like = {P1,P2}.
%
% INPUT
% none.
%
% OUTPUT
% none; one printed line per check.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_lpi_ineq_sop
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
% Initial coding MMP, 10/08/2026
% MMP, 10/21/2026: Part (7), for the count before the build and the support
%                read from the operators of 'poscopvar' (lpi_ineq_sop).
% MMP, 10/28/2026: Header to the header protocol (what is checked, inputs,
%                outputs), the licence block added, step comments in the
%                body; code unchanged. The deleted header text was the
%                paper and generator of the heat LPI and the degrees the
%                reader chose in (1) and (2); part (6) is now listed. The
%                help text of the 10/21/2026 entry is in (7), without
%                markers.

% Silence the canonical-form warnings; MOSEK, else SeDuMi; the pass count
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
sopts = struct('solver','mosek','simplify',false);
if isempty(which('mosekopt')),  sopts.solver = 'sedumi';   end
npass = 0;

%% (1) V'*V + L'*L
% Pt on L2[0,1]: lower kernel 2 - s, upper kernel 2 - t
vars = {'s1'};  dom = [0 1];
ZL = {[0;1]};   ZR = {[0;1]};
params = reshape({sparse(2,2),sparse([2 0;-1 0]),sparse([2 -1;0 0])},[3 1 1]);
Pt = copvar({sopvar(params,struct('out',{vars},'in',{vars}),ZL,ZR,struct('out',dom,'in',dom),[1 1])});
% Pt >= 0, solved; the residual
prog = lpiprogram_sop(Pt);
[prog,Pop,inf1] = lpi_ineq_sop(prog,Pt);
prog = lpisolve(prog,sopts);
rel_b = cx_resid(prog);
if rel_b>1e-6
    error('test_lpi_ineq_sop: (1) V''V + L''L is not certified (rel_b %.2e).',rel_b);
end
fprintf('  passed: (1) V''V + L''L certified, rel_b %.1e, lift %s weight %s, terms %s offsets %s, ndec %d\n',...
        rel_b,mat2str(inf1.degrees.D),mat2str(inf1.degrees.w),mat2str(inf1.terms.codes),mat2str(inf1.terms.offsets),numel(Pop.Zd));
npass = npass+1;

%% (2) the 1-D heat LPI through two inequalities
% The heat PIE; kappa below and above pi^2
pie = heatNd_pie(1,0);  vars = pie.vars;    dom = pie.dom;    ep = 0.1;
res = zeros(1,2);   kap = [9, 11];
for i = 1:2
    % P free, P'*T self-adjoint, the two inequalities; solve
    prog = lpiprogram_sop(vars,dom);
    [prog,P] = lpivar_cdopvar(prog,1,{vars},dom,1);
    PT = P'*pie.T;  PA = P'*pie.A0;  TT = pie.T'*pie.T;
    E1 = PT + PT' - 2*ep^2*TT;                  % 2 (P*T - ep^2 T*T), self-adjoint form
    prog = lpi_eq_sop(prog,PT-PT','symmetric');
    [prog,~,iR] = lpi_ineq_sop(prog,E1);
    X = -(PA + PA' + kap(i)*(PT + PT'));        % -(P*A + A*P + 2 kappa P*T) >= 0
    [prog,~,iQ] = lpi_ineq_sop(prog,X);
    prog = lpisolve(prog,sopts);
    res(i) = cx_resid(prog);
    pinf(i) = prog.solinfo.info.pinf; %#ok<AGROW>
    % The degrees and terms of both certificates, printed once
    if i==1
        fprintf('  (2) degrees: R lift %s weight %s; Q lift %s weight %s; terms %s offsets %s\n',...
                mat2str(iR.degrees.D),mat2str(iR.degrees.w),mat2str(iQ.degrees.D),mat2str(iQ.degrees.w),...
                mat2str(iQ.terms.codes),mat2str(iQ.terms.offsets));
    end
end
% kappa = 9 feasible, kappa = 11 rejected
if ~(res(1)<=1e-6 && pinf(1)==0)
    error('test_lpi_ineq_sop: (2) kappa = 9 should be feasible (rel_b %.2e, pinf %d).',res(1),pinf(1));
end
if pinf(2)==0 && res(2)<=1e-6
    error('test_lpi_ineq_sop: (2) kappa = 11 > pi^2 was accepted (rel_b %.2e).',res(2));
end
fprintf('  passed: (2) heat LPI: kappa 9 feasible (rel_b %.1e), kappa 11 rejected (pinf %d, rel_b %.1e)\n',res(1),pinf(2),res(2));
npass = npass+1;

%% (3) 2-D: T'*T >= 0
% T'*T of the 2-D heat PIE, solved; the residual
pie2 = heatNd_pie(2,0);
TT2 = pie2.T'*pie2.T;
prog = lpiprogram_sop(pie2.vars,pie2.dom);
[prog,Pop2,inf3] = lpi_ineq_sop(prog,TT2);
prog = lpisolve(prog,sopts);
rel_b = cx_resid(prog);
if rel_b>1e-6
    error('test_lpi_ineq_sop: (3) 2-D T''T is not certified (rel_b %.2e).',rel_b);
end
fprintf('  passed: (3) 2-D T''T certified, rel_b %.1e, lift %s weight %s, terms %s, ndec %d\n',...
        rel_b,mat2str(inf3.degrees.D),mat2str(inf3.degrees.w),mat2str(inf3.terms.codes),numel(Pop2.Zd));
npass = npass+1;

%% (4) options
% Each option variant on Pt, built and not solved
prog = lpiprogram_sop(Pt);
[~,Pa] = lpi_ineq_sop(prog,Pt,struct('psatz','none'));
[~,Pb] = lpi_ineq_sop(prog,Pt,struct('psatz','faces'));
[~,Pc,ic] = lpi_ineq_sop(prog,Pt,struct('deg',struct('int',1,'mult',1)));
[~,Pd] = lpi_ineq_sop(prog,Pt,struct('prune',false));
[~,Pe] = lpi_ineq_sop(prog,Pt,struct('dD',0,'dw',0));
% The decision counts ordered; a given deg bypasses the reader
if ~(numel(Pa.Zd)<numel(Pop.Zd) && numel(Pb.Zd)>numel(Pop.Zd) && numel(Pe.Zd)<numel(Pop.Zd) && numel(Pd.Zd)>=numel(Pop.Zd))
    error('test_lpi_ineq_sop: (4) the decision counts of the option variants are not ordered as expected.');
end
if ~isfield(ic.degrees,'given')
    error('test_lpi_ineq_sop: (4) a given ''deg'' should bypass the reader.');
end
fprintf('  passed: (4) options: none %d < default %d < faces %d; given deg %d; prune off %d; dD = dw = 0 %d\n',...
        numel(Pa.Zd),numel(Pop.Zd),numel(Pb.Zd),numel(Pc.Zd),numel(Pd.Zd),numel(Pe.Zd));
npass = npass+1;

%% (5) legacy dispatch
% A positive 'dopvar' through lpi_ineq and through lpi_ineq_sop
pvar s th;
prog1 = lpiprogram(s,th,[0 1]);
[prog1,Pl] = poslpivar(prog1,[1;1],1);
p1 = lpi_ineq(prog1,Pl);
p2 = lpi_ineq_sop(prog1,Pl);
if ~isequal(numel(p1.decvartable),numel(p2.decvartable)) || p1.expr.num~=p2.expr.num
    error('test_lpi_ineq_sop: (5) the legacy dispatch differs from lpi_ineq.');
end
fprintf('  passed: (5) legacy dopvar: lpi_ineq and lpi_ineq_sop agree (%d decision variables, %d expressions)\n',...
        numel(p2.decvartable),p2.expr.num);
npass = npass+1;

%% (6) opts.like: the support of another operator, laid out over P's spaces
% Pt sized from the support of Pt'*Pt
Pq = Pt'*Pt;                                % a larger support than Pt's
[~,~,iq] = get_lift_degs(Pq);
[~,Pl,il] = lpi_ineq_sop(prog,Pt,struct('like',Pq));
if ~(isequal(il.degrees.D,iq.D) && isequal(il.degrees.w,iq.w) && numel(Pl.Zd)>numel(Pop.Zd))
    error('test_lpi_ineq_sop: (6) ''like'' should size P from the support of the given operator.');
end
% A 'double' opts.like is refused
try
    lpi_ineq_sop(prog,Pt,struct('like',1));     e6 = '';
catch ME,   e6 = ME.message;
end
if isempty(e6)
    error('test_lpi_ineq_sop: (6) a non-container ''like'' should be refused.');
end
fprintf('  passed: (6) like: lift %s weight %s from the given operator (own %s %s), ndec %d > %d; non-container refused\n',...
        mat2str(il.degrees.D),mat2str(il.degrees.w),mat2str(inf1.degrees.D),mat2str(inf1.degrees.w),numel(Pl.Zd),numel(Pop.Zd));
npass = npass+1;

%% (7) the count before the build; the support of poscopvar's operators     % MMP, 10/21/2026
% (a) info.ndec, counted before anything is declared, is the growth of
%     prog.decvartable (N^2 per Gram block) in the calls below, and Pop
%     has N(N+1)/2 distinct names per block in (1), (3), (4), (6) and below;
% (b) a target that depends on an operator of 'poscopvar' is sized from it
%     (prog.sop_pos), exactly as with opts.like = that operator; without
%     the record it is sized from the target, with more variables;
% (c) a target that depends on two such operators is sized from both,
%     exactly as with opts.like = {P1,P2}.
% Two operators of poscopvar; targets on P1 and on P1 and P2
pie1 = heatNd_pie(1,0);
prog = lpiprogram_sop(pie1.vars,pie1.dom);
[prog,P1] = poscopvar(prog,1,{pie1.vars},pie1.dom,2);
[prog,P2] = poscopvar(prog,1,{pie1.vars},pie1.dom,3);
X1 = -(pie1.T'*P1*pie1.A0 + pie1.A0'*P1*pie1.T);
X12 = X1 - pie1.T'*P2*pie1.T;
% The five calls: from the record, from opts.like, without the record
n0 = numel(prog.decvartable);
[p7a,N7a,i7a] = lpi_ineq_sop(prog,X1);
[p7b,N7b,i7b] = lpi_ineq_sop(prog,X1,struct('like',P1));
[p7c,N7c,i7c] = lpi_ineq_sop(rmfield(prog,'sop_pos'),X1);
[p7d,N7d,i7d] = lpi_ineq_sop(prog,X12);
[p7e,N7e,i7e] = lpi_ineq_sop(prog,X12,struct('like',{{P1,P2}}));
% (a) the growth and the distinct names against the counts
grow = cellfun(@(p) numel(p.decvartable),{p7a,p7b,p7c,p7d,p7e}) - n0;
nd = [i7a.ndec, i7b.ndec, i7c.ndec, i7d.ndec, i7e.ndec];
pairs = {inf1,Pop; inf3,Pop2; ic,Pc; il,Pl; i7a,N7a; i7b,N7b; i7c,N7c; i7d,N7d; i7e,N7e};
cnt = cellfun(@(i,P) sum(i.N.*(i.N+1)/2)==numel(P.Zd),pairs(:,1),pairs(:,2));
if ~isequal(grow,nd) || ~all(cnt)
    error('test_lpi_ineq_sop: (7a) the count before the build is not what was declared: grew %s, counted %s; names wrong in cases %s.',...
          mat2str(grow),mat2str(nd),mat2str(find(~cnt)'));
end
% (b) and (c): the sources and the degrees
if ~(strcmp(i7a.sized_from,'1 operator(s) of poscopvar') && isequal(i7a.deg,i7b.deg) ...
        && numel(N7a.Zd)==numel(N7b.Zd) && strcmp(i7c.sized_from,'the operator') && numel(N7c.Zd)>numel(N7a.Zd))
    error('test_lpi_ineq_sop: (7b) a target on an operator of poscopvar should be sized from it.');
end
if ~(strcmp(i7d.sized_from,'2 operator(s) of poscopvar') && isequal(i7d.deg,i7e.deg) && numel(N7d.Zd)==numel(N7e.Zd))
    error('test_lpi_ineq_sop: (7c) a target on two operators of poscopvar should be sized from both.');
end
fprintf('  passed: (7) counted = declared %s; from P1 %d (= like %d) against %d from the target; from P1, P2 %d (= like %d)\n',...
        mat2str(nd),i7a.ndec,i7b.ndec,i7c.ndec,i7d.ndec,i7e.ndec);
npass = npass+1;

% The summary line
fprintf('\ntest_lpi_ineq_sop passed (%d checks).\n',npass);
end
