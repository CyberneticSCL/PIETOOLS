function [prog,L,gam,P,Z,W,info] = PIETOOLS_H2_estimator_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,L,GAM,P,Z,W,INFO] = PIETOOLS_H2_ESTIMATOR_SOP(PIE,SETTINGS) computes
% an observer gain L for the 1-D PIE
%
%       T x_t = A x + B1 w,      z = C1 x,      y = C2 x + D21 w
%
% and a bound GAM on the H2 norm from w to the estimation error of z, by
% the LPI, posed on 'cdopvar' operators and solved as an SDP that minimizes
% gam,
%
%   P >= eppos I,  Z and W free,  gam >= trace(W),
%   [-gam I, C1; 
%   C1', (T'*P)*A + (A'*P)*T + (T'*Z)*C2 + (C2'*Z')*T] <= 0,
%   [W, -(B1'*P + D21'*Z');
%   -(B1'*P + D21'*Z')', P] >= 0,
%
% with P a positive operator, Z and W free operators, and L the gain of
% bounded degree that minimizes ||P*L - Z||^2, an approximation of
% P^{-1}*Z. While the solver returns UNKNOWN or INFEASIBLE, the program is
% rebuilt with the degree of the slacks raised by one, up to
% settings.sop.max_raise times. PIEs with an algebraic constraint (K, Kw or
% Ku with rows, e.g. from periodic boundary conditions) are not supported.
%
% INPUT
% - PIE:      'pie_struct' object of a PIE in one spatial variable,
%             initialized when its T is not yet a 'copvar' object; T, A,
%             B1, C1, C2 and D21 'copvar' objects; Tw and D11 zero or
%             absent; K, Kw, Ku empty or without rows; its other
%             operators are not used; a PIE in
%             two or more variables is an error;
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
% - L:        'copvar' object, the observer gain from y to the state ([] when
%             the solver returned no solution);
% - gam:      double, the certified bound; NaN when the solver returned no
%             value;
% - P, Z, W:  'copvar' objects, the solved P, Z and W;
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded', 'unknown' or 'residual' (cuADMM,
%                          which decides no feasibility: never certified);
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               form ('est'), gam, gam_lpi (equal to gam), nv;
%               P, Z, W    as the outputs;
%               G          'copvar' object, the gain L; inv  'struct', the
%                          diagnostics of the computation of L (output INFO
%                          of GETOBSERVER_SOP, method 'qp': the field qp;
%                          no inverse of P is formed);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          the last raise and the loop history, the slack
%                          mode, its degrees and terms, the size of the SDP,
%                          and the dual kernel of the negativity constraint
%                          with its 1-norm.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_H2_estimator_sop
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
% DJ, 11/21/2026: Error for a PIE with an algebraic constraint (K, Kw, Ku
%                with rows);

% Default settings; check the PIE
if nargin<2 || isempty(settings),  settings = lpisettings_sop('light'); end
if ~isa(PIE,'pie_struct')
    error('PIETOOLS_H2_estimator_sop:pie','The PIE should be a ''pie_struct'' object.');
end
PIE = exec_pie_sop(PIE);                    % every operator a 'copvar' or 0
label = 'H2_estimator_sop';

% Complete the options of the executives and read the tools of the
% dimension of the PIE (storage of P, margins, degrees of Z); boundary
% disturbances, a PIE in two or more variables and D11 are errors
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension
if ~(PIE.Tw==0),    error('PIETOOLS_H2_estimator_sop:Tw','The H2 LPI cannot be solved with disturbances at the boundary.');  end
if nv>=2,           error('PIETOOLS_H2_estimator_sop:dim2','H2 estimation of PIEs in two or more variables is not supported.'); end
if ~(PIE.D11==0),   error('PIETOOLS_H2_estimator_sop:D11','Feedthrough D11 is not supported.');          end
eppos = tl.eppos;

% The PIE operators, 'copvar' objects, and the identity on z
Tm = PIE.T;     Am = PIE.A;
% The constraint K x = 0 on the state, if any (empty for most PIEs)
Kcon = PIE.K;   hasK = isa(Kcon,'copvar') && size(Kcon,1)>0;
% PIEs with an algebraic constraint (K, Kw, Ku with rows) are not supported
if any(cellfun(@(n) isa(PIE.(n),'copvar') && size(PIE.(n),1)>0,{'K','Kw','Ku'})) % DJ, 11/21/2026
    error('PIETOOLS_H2_estimator_sop:algebraicConstraint','The H2 estimator synthesis LPI is not supported for PIEs with algebraic constraints.');
end
% A constraint involving the input (Kw nonzero) is not supported; u is zero
if hasK && isa(PIE.Kw,'copvar') && any(cellfun(@(B) ~isempty(B) && any(cellfun(@nnz,B.params(:))),PIE.Kw.C(:)))
    error('PIETOOLS_H2_estimator_sop:constraintInputs','Algebraic constraints involving inputs are currently not supported.');
end
Bw = PIE.B1;    Cz = PIE.C1;
C2m = PIE.C2;   D21m = PIE.D21;
Iz = tl.ident(PIE.C1,1);

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
gam = cl.obj;

% The solved P, Z, W, as 'copvar' objects, and the gain L, an
% approximation of P^{-1}*Z
L = [];     P = [];     Z = [];     W = [];     Gc = [];    iinv = [];
try
    Psol = getsol_lpivar_sop(prog,aux.Pm);  Zsol = getsol_lpivar_sop(prog,aux.Zm);  Wsol = getsol_lpivar_sop(prog,aux.Wm);
    % L of bounded degree minimizing ||P*L - Z||^2
    P = Psol;      Z = Zsol;      W = Wsol;
    [Gc,~,iinv] = getObserver_sop(Psol,Zsol);       L = Gc;
catch ME
    if st.sop.verbose,  fprintf('  [%s] output extraction: %s\n',label,ME.message);  end
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info_sop(prog,aux,st,label);
info.form = 'est';      info.gam = gam;     info.gam_lpi = gam;     info.nv = nv;
info.P = P;     info.Z = Z;     info.W = W;     info.G = Gc;    info.inv = iinv;
if st.sop.verbose
    fprintf('  [%s] %s: value %.8g\n',label,cl.status,gam);
end

    function [prog,aux] = assemble(dw)
        % [PROG,AUX] = ASSEMBLE(DW) builds the H2 estimator LPI with the
        % slack degrees raised by settings.sop.dw_P + DW.
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
        % P >= eppos I; Z free from y to the state; W free on the space of w
        [spx,dmx] = copvar_space_list(Tm,'out');
        [prog,Pm] = tl.storage(prog,Tm,'out');
        Pm = Pm + tl.eye(1,eppos);
        [spy,dmy] = copvar_space_list(C2m,'out');
        [prog,Zm] = lpivar_cdopvar(prog,struct('out',dmx,'in',dmy),struct('out',{spx},'in',{spy}),tl.dom,tl.zdeg);
        [spw,dmw] = copvar_space_list(Bw,'in');
        [prog,Wm] = lpivar_cdopvar(prog,dmw,spw,tl.dom,tl.zdeg);
        % The negativity (Dneg) and positivity (Dpos) constraints
        % The state block, with K'*Z_K + Z_K'*K if K is given
        X22 = (Tm'*Pm)*Am + (Am'*Pm)*Tm + (Tm'*Zm)*C2m + (C2m'*Zm')*Tm;
        if hasK,    [prog,X22] = exec_kmult_sop(prog,X22,Kcon,Tm,tl.dom,Pm);    end
        Dneg = [-(gamv*Iz),  Cz;
                Cz',         X22];
        D12 = Bw'*Pm + on_registry_sop(D21m',Zm')*Zm';
        Dpos = [Wm,       -D12;
                -D12',    Pm];
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
        aux = struct('obj',gamv,'Rm',[],'Qm',[],'Pm',Pm,'Zm',Zm,'Wm',Wm,'li',li,'form','est');
    end
end
