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
%   [-gam I, C1; C1', (T'*P)*A + (A'*P)*T + (T'*Z)*C2 + (C2'*Z')*T] <= 0,
%   [W, -(B1'*P + D21'*Z'); -(B1'*P + D21'*Z')', P] >= 0,
%
% with L = P^{-1}*Z. While the solver returns UNKNOWN or INFEASIBLE, the
% program is rebuilt with the degree of the slacks raised by one, up to
% settings.sop.max_raise times.
%
% INPUT
% - PIE:      'pie_struct' object of a 1-D PIE, with T, A, B1, C1, C2 and D21
%             'opvar' objects, or 'copvar' objects (a 1-D container PIE,
%             EXEC_PIE_SOP); Tw and D11 zero or absent; a PIE in two or
%             more variables is an error;
% - settings: (optional) 'struct' of LPI settings, as returned by lpisettings
%             or lpisettings_sop:
%               dd1, dd12, options1, options12, override1   degrees of P;
%               ddZ     degrees of Z and W;
%               eppos, eppos2   margins of P;
%               sos_opts  solver options;
%               sop     container options (EXEC_SOP_SETTINGS): dw_P, the
%                       degree raise of the slacks, and max_raise.
%             Default lpisettings('light').
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - L:        'opvar' object, the observer gain from y to the state ([] when
%             the solver returned no solution);
% - gam:      double, the certified bound; NaN when the solver returned no
%             value;
% - P, Z, W:  'opvar' objects, the solved P, Z and W;
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded' or 'unknown';
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               form ('est'), gam, gam_lpi (equal to gam), nv;
%               P, Z, W    as the outputs;
%               G          'copvar' object, the gain L; inv  'struct' of the
%                          inversion of P (GETOBSERVER_SOP);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          as returned by EXEC_INFO_SOP.
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
% MMP, 10/28/2026: The body (a call of private/h2_build_sop and a local
%                safe_op whose branch, a call of the nonexistent
%                cop2opvar_sop, was never reached) is replaced by the
%                builder's 'est' code (MMP, 10/08/2026, 10/09/2026), the
%                program unchanged. Header to the header protocol.
% MMP, 10/29/2026: The PIE goes through EXEC_PIE_SOP (INITIALIZE for an
%                opvar PIE, as before; a 1-D container PIE accepted), and
%                every PIE in two or more variables is rejected, not only a
%                2-D one: the gain is extracted in 1-D only.
% MMP, 10/29/2026: Both slacks take opts.cap of EXEC_TOOLS_SOP capfor, each
%                laid over its own target's variables (the per-direction
%                bounds of settings.sop.cap_delta); empty by default, so
%                every program is unchanged unless cap_delta is set.
% MMP, 10/29/2026: st.sop takes tl.sop of EXEC_TOOLS_SOP: each per-dimension
%                vector (dw_feas, dw_P, dw_Q, max_raise, max_raise_feas,
%                cap_delta) reduced to its entry for the dimension of the
%                PIE, so the settings tiers may hold [1-D 2-D 3-D] values;
%                scalar settings are unchanged.

% % % BEGIN body replaced by MMP, 10/28/2026 - everything from here to the
% % % END marker at the foot of the file is new; see the header entries
% % % above for what was deleted.
% Default settings; check the PIE
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
if ~isa(PIE,'pie_struct')
    error('PIETOOLS_H2_estimator_sop:pie','The PIE should be a ''pie_struct'' object.');
end
% PIE = initialize(PIE);                                                    % MMP, 10/29/2026 (was)
PIE = exec_pie_sop(PIE);                    % opvar or container PIE        % MMP, 10/29/2026
label = 'H2_estimator_sop';

% Complete the container options and read the tools of the PIE's dimension
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension (EXEC_TOOLS_SOP)        % MMP, 10/29/2026
if ~(PIE.Tw==0),    error('PIETOOLS_H2_estimator_sop:Tw','The H2 LPI cannot be solved with disturbances at the boundary.');  end
% if nv==2,           error('PIETOOLS_H2_estimator_sop:dim2','H2 estimation of 2-D PIEs is not supported.'); end % MMP, 10/29/2026 (was)
if nv>=2,           error('PIETOOLS_H2_estimator_sop:dim2','H2 estimation of PIEs in two or more variables is not supported.'); end % MMP, 10/29/2026
if ~(PIE.D11==0),   error('PIETOOLS_H2_estimator_sop:D11','Feedthrough D11 is not supported.');          end
eppos = tl.eppos;

% The PIE operators as containers, and the identity on z
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);
Bw = op2copvar_sop(PIE.B1); Cz = op2copvar_sop(PIE.C1);
C2m = op2copvar_sop(PIE.C2);    D21m = op2copvar_sop(PIE.D21);
Iz = tl.ident(PIE.C1,1);

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
gam = cl.obj;

% The solved P, Z, W, and the gain L = P^{-1} Z on the container path
L = [];     P = [];     Z = [];     W = [];     Gc = [];    iinv = [];
try
    Psol = getsol_lpivar_sop(prog,aux.Pm);  Zsol = getsol_lpivar_sop(prog,aux.Zm);  Wsol = getsol_lpivar_sop(prog,aux.Wm);
    P = tl.toop(Psol);      Z = tl.toop(Zsol);      W = tl.toop(Wsol);
    [Gc,L,iinv] = getObserver_sop(Psol,Zsol);
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
        %               ('cdopvar' P, Z, W), li (slack info of LPI_INEQ_SOP
        %               for the negativity constraint), form.
        % New program; record the equality rows when the dual is kept
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        if st.sop.keepdual,     prog.sopeq = {};    end
        % gam >= 0 as the objective
        [prog,gamv] = lpidecvar(prog,'gam');
        prog = lpi_ineq(prog,gamv);
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
        Dneg = [-(gamv*Iz),  Cz;
                Cz',         (Tm'*Pm)*Am + (Am'*Pm)*Tm + (Tm'*Zm)*C2m + (C2m'*Zm')*Tm];
        D12 = Bw'*Pm + on_registry_sop(D21m',Zm')*Zm';
        Dpos = [Wm,       -D12;
                -D12',    Pm];
        % Dneg <= 0 and Dpos >= 0: the stock slacks and equalities, or the
        % slacks sized from P
        if strcmp(st.sop.slack,'stock')
            [prog,Nneg] = tl.slack(prog,Dneg);
            [prog,tag] = lpi_eq_sop(prog,Nneg + Dneg,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
            [prog,Npos] = tl.slack(prog,Dpos);
            prog = lpi_eq_sop(prog,Npos - Dpos,'symmetric');
        else
%           o = struct('like',Pm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin); % MMP, 10/29/2026 (was)
            o = struct('like',Pm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin,'cap',tl.capfor(Dneg)); % MMP, 10/29/2026
            [prog,~,li] = lpi_ineq_sop(prog,-Dneg,o);
%           [prog,~,~] = lpi_ineq_sop(prog,Dpos,o);                         % MMP, 10/29/2026 (was)
            op = o;     op.cap = tl.capfor(Dpos);       % its own variables % MMP, 10/29/2026
            [prog,~,~] = lpi_ineq_sop(prog,Dpos,op);                        % MMP, 10/29/2026
        end
        % gam >= trace(W)
        prog = lpi_ineq(prog,gamv - trace_rn_sop(Wm));
        aux = struct('obj',gamv,'Rm',[],'Qm',[],'Pm',Pm,'Zm',Zm,'Wm',Wm,'li',li,'form','est');
    end
% % % END body replaced by MMP, 10/28/2026
end
