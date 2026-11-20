function [prog,K,gam,P,Z,W,info] = PIETOOLS_H2_control_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,K,GAM,P,Z,W,INFO] = PIETOOLS_H2_CONTROL_SOP(PIE,SETTINGS) computes a
% state feedback u = K x for the 1-D PIE
%
%       T x_t = A x + B1 w + B2 u,      z = C1 x + D12 u
%
% and a bound GAM on the H2 norm from w to z of the closed loop, by the LPI,
% posed on 'cdopvar' operators and solved as an SDP that minimizes gam,
%
%   P >= eppos I,  Z and W free,  gam >= trace(W),
%   [-gam I, B1'; B1, (T*P)*A' + (A*P)*T' + (B2*Z)*T' + (T*Z')*B2'] <= 0,
%   [W, C1*P + D12*Z; (C1*P + D12*Z)', P] >= 0,
%
% with P a positive operator, Z and W free operators, and K the gain of
% bounded degree that minimizes ||Z - K*P||^2, an approximation of
% Z*P^{-1}. While the solver returns UNKNOWN or INFEASIBLE, the program is
% rebuilt with the degree of the slacks raised by one, up to
% settings.sop.max_raise times.
%
% INPUT
% - PIE:      'pie_struct' object of a PIE in one spatial variable,
%             initialized when its T is not yet a 'copvar' object; T, A,
%             B1, B2, C1 and D12 'copvar' objects; Tw, Tu and D11 zero or
%             absent; K empty (a constraint K x = 0 on the state, e.g.
%             periodic BCs, is an error); its other operators are not
%             used; a PIE in two or more variables is an error;
% - settings: (optional) 'struct' of LPI settings:
%               lpi.pos   the storage row of P;
%               lpi.var   field deg, a 1x3 'cell': deg{1} the degrees of
%                         Z and W;
%               eppos     the margin of P on every space (eppos2, when
%                         set, the margin on L2);
%               sos_opts  solver options;
%               sop       the options of the executives, among them dw_P,
%                         the degree raise of the slacks, and max_raise.
%             Default: the 'light' tier.
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - K:        'copvar' object, the feedback gain from x to u ([] when the
%             solver returned no solution);
% - gam:      double, the certified closed-loop bound; NaN when the solver
%             returned no value;
% - P, Z, W:  'copvar' objects, the solved P, Z and W;
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded', 'unknown' or 'residual' (cuADMM,
%                          which decides no feasibility: never certified);
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               form ('ctrl'), gam, gam_lpi (equal to gam), nv;
%               P, Z, W    as the outputs;
%               G          'copvar' object, the gain K; inv  'struct', the
%                          diagnostics of the computation of K;
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          the last raise and the loop history, the slack
%                          mode, its degrees and terms, the size of the SDP,
%                          and the dual kernel of the negativity constraint
%                          with its 1-norm.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_H2_control_sop
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

% Default settings; check the PIE
if nargin<2 || isempty(settings),  settings = lpisettings_sop('light'); end
if ~isa(PIE,'pie_struct')
    error('PIETOOLS_H2_control_sop:pie','The PIE should be a ''pie_struct'' object.');
end
PIE = exec_pie_sop(PIE);                    % every operator a 'copvar' or 0
label = 'H2_control_sop';

% Complete the options of the executives and read the tools of the
% dimension of the PIE (storage of P, margins, degrees of Z); a PIE in two
% or more variables, boundary inputs, D11 and a constraint K are errors
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension
if ~(PIE.Tw==0),    error('PIETOOLS_H2_control_sop:Tw','The H2 LPI cannot be solved with disturbances at the boundary.');  end
if nv>=2,           error('PIETOOLS_H2_control_sop:dim2','H2 control of PIEs in two or more variables is not supported.');    end
if ~(PIE.Tu==0),    error('PIETOOLS_H2_control_sop:Tu','Inputs at the boundary are not supported.');   end
if ~(PIE.D11==0),   error('PIETOOLS_H2_control_sop:D11','Feedthrough D11 is not supported.');          end
% A constraint K x = 0 on the state (e.g. periodic BCs) cannot be imposed
% in this dual LPI, and a dual LPI without it is not sufficient
if isa(PIE.K,'copvar') && size(PIE.K,1)>0
    error('PIETOOLS_H2_control_sop:K','PIEs with an algebraic constraint K are not supported.');
end
eppos = tl.eppos;

% The PIE operators, 'copvar' objects, and the identity on w
Tm = PIE.T;     Am = PIE.A;
Bw = PIE.B1;    Cz = PIE.C1;
B2m = PIE.B2;   D12m = PIE.D12;
Iw = tl.ident(PIE.B1,2);

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
gam = cl.obj;

% The solved P, Z, W, as 'copvar' objects, and the gain K, an
% approximation of Z*P^{-1}
K = [];     P = [];     Z = [];     W = [];     Gc = [];    iinv = [];
try
    Psol = getsol_lpivar_sop(prog,aux.Pm);  Zsol = getsol_lpivar_sop(prog,aux.Zm);  Wsol = getsol_lpivar_sop(prog,aux.Wm);
    % K of bounded degree minimizing ||Z - K*P||^2
    P = Psol;      Z = Zsol;      W = Wsol;
    [Gc,~,iinv] = getController_sop(Psol,Zsol);     K = Gc;
catch ME
    if st.sop.verbose,  fprintf('  [%s] output extraction: %s\n',label,ME.message);  end
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info_sop(prog,aux,st,label);
info.form = 'ctrl';     info.gam = gam;     info.gam_lpi = gam;     info.nv = nv;
info.P = P;     info.Z = Z;     info.W = W;     info.G = Gc;    info.inv = iinv;
if st.sop.verbose
    fprintf('  [%s] %s: value %.8g\n',label,cl.status,gam);
end

    function [prog,aux] = assemble(dw)
        % [PROG,AUX] = ASSEMBLE(DW) builds the H2 synthesis LPI with the slack
        % degrees raised by settings.sop.dw_P + DW.
        % INPUT  dw:    nonnegative integer, the raise of this pass.
        % OUTPUT prog:  'struct', the LPI program;
        %        aux:   'struct' with fields obj (the 'dpvar' gam), Pm, Zm, Wm
        %               ('cdopvar' P, Z, W), li (the info of the slack of
        %               the negativity constraint), form.
        % New program; record the equality rows when the dual is kept
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);
        if st.sop.keepdual,     prog.sopeq = {};    end
        % gam >= 0 as the objective
        [prog,gamv] = lpidecvar(prog,'gam');
        prog = lpi_ineq_sop(prog,gamv);
        prog = lpisetobj(prog,gamv);
        % P >= eppos I; Z free from the state to u; W free on the space of z
        [spx,dmx] = copvar_space_list(Am,'out');
        [prog,Pm] = tl.storage(prog,Am,'out');
        Pm = Pm + tl.eye(1,eppos);
        [spu,dmu] = copvar_space_list(B2m,'in');
        [prog,Zm] = lpivar_cdopvar(prog,struct('out',dmu,'in',dmx),struct('out',{spu},'in',{spx}),tl.dom,tl.zdeg);
        [spz,dmz] = copvar_space_list(Cz,'out');
        [prog,Wm] = lpivar_cdopvar(prog,dmz,spz,tl.dom,tl.zdeg);
        % The negativity (Dneg) and positivity (Dpos) constraints
        Dneg = [-(gamv*Iw),  Bw';
                Bw,          (Tm*Pm)*Am' + (Am*Pm)*Tm' + (B2m*Zm)*Tm' + (Tm*Zm')*B2m'];
        Dp12 = Cz*Pm + on_registry_sop(D12m,Zm)*Zm;
        Dpos = [Wm,     Dp12;
                Dp12',  Pm];
        % Dneg <= 0 and Dpos >= 0: the 'stock' slacks and equalities, or the
        % slacks sized from P
        if strcmp(st.sop.slack,'stock')
            [prog,Nneg] = tl.slack(prog,Dneg);
            [prog,tag] = lpi_eq_sop(prog,Nneg + Dneg,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
            [prog,Npos] = tl.slack(prog,Dpos);
            prog = lpi_eq_sop(prog,Npos - Dpos,'symmetric');
        else
            o = struct('like',Pm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin,'cap',tl.capfor(Dneg));
            [prog,~,li] = lpi_ineq_sop(prog,-Dneg,o);
            op = o;     op.cap = tl.capfor(Dpos);       % its own variables
            [prog,~,~] = lpi_ineq_sop(prog,Dpos,op);
        end
        % gam >= trace(W)
        prog = lpi_ineq_sop(prog,gamv - trace_rn_sop(Wm));
        aux = struct('obj',gamv,'Rm',[],'Qm',[],'Pm',Pm,'Zm',Zm,'Wm',Wm,'li',li,'form','ctrl');
    end
end
