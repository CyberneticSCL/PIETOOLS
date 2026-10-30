function [prog,K,gam,P,Z,info] = PIETOOLS_Hinf_control_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,K,GAM,P,Z,INFO] = PIETOOLS_HINF_CONTROL_SOP(PIE,SETTINGS) computes
% a state feedback u = K x for the 1-D PIE
%
%       T x_t = A x + B1 w + B2 u,      z = C1 x + D11 w + D12 u
%
% and a bound GAM on the L2 gain from w to z of the closed loop, by the
% LPI, posed on 'cdopvar' operators and solved as an SDP that minimizes gam,
%
%   P >= eppos I,  Z free,  CZ = C1*P + D12*Z,  AZ = A*P + B2*Z,
%   [-gam I, D11, CZ*T'; D11', -gam I, B1'; T*CZ', B1, AZ*T' + T*AZ'] <= 0,
%
% with K = Z*P^{-1}. While the solver returns UNKNOWN or INFEASIBLE, the
% program is rebuilt with the degree of the negativity slack raised by
% one, up to settings.sop.max_raise times.
%
% INPUT
% - PIE:      'pie_struct' object of a 1-D PIE, with T, A, B1, B2, C1, D11
%             and D12 'opvar' objects, or 'copvar' objects (a 1-D container
%             PIE, EXEC_PIE_SOP); Tw and Tu zero or absent (no disturbance
%             or input at the boundary); a PIE in two or more variables is
%             an error;
% - settings: (optional) 'struct' of LPI settings, as returned by lpisettings
%             or lpisettings_sop:
%               dd1, dd12, options1, options12, override1   degrees of P;
%               ddZ     degrees of Z;
%               eppos, eppos2   margins of P;
%               sos_opts  solver options;
%               sop     container options (EXEC_SOP_SETTINGS): dw_P, the
%                       degree raise of the slack, max_raise, and
%                       gam_fixed (a positive double: solve the
%                       feasibility LPI at gam = gam_fixed).
%             Default lpisettings('light').
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - K:        'opvar' object, the feedback gain from x to u ([] when the
%             solver returned no solution);
% - gam:      double, the certified closed-loop bound; with gam_fixed,
%             gam_fixed when certified and Inf otherwise; NaN when the
%             solver returned no value;
% - P, Z:     'opvar' objects, the solved P and Z;
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded' or 'unknown';
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               form ('ctrl'), gam, gam_fixed ([]), nv;
%               P, Z       'opvar' objects; Psol, Zsol the 'copvar' objects;
%               K          'cdopvar' object, the LPI operator;
%               G          'copvar' object, the gain K; inv  'struct' of the
%                          inversion of P (GETCONTROLLER_SOP);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          as returned by EXEC_INFO_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_Hinf_control_sop
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
% MMP, 10/28/2026: The body (a call of private/synth_build_sop) is replaced
%                by the builder's 'ctrl' code (MMP, 10/08/2026, 10/09/2026),
%                the program unchanged. Header to the header protocol.
% MMP, 10/29/2026: The PIE goes through EXEC_PIE_SOP (INITIALIZE for an
%                opvar PIE, as before; a 1-D container PIE accepted), and
%                every PIE in two or more variables is rejected, not only a
%                2-D one: the controller gain is extracted in 1-D only
%                (GETCONTROLLER_SOP).
% MMP, 10/29/2026: The slack takes opts.cap = capfor(K) of EXEC_TOOLS_SOP,
%                the per-direction bounds of settings.sop.cap_delta; empty
%                by default, so every program is unchanged unless cap_delta
%                is set.
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
    error('PIETOOLS_Hinf_control_sop:pie','The PIE should be a ''pie_struct'' object.');
end
% PIE = initialize(PIE);                                                    % MMP, 10/29/2026 (was)
PIE = exec_pie_sop(PIE);                    % opvar or container PIE        % MMP, 10/29/2026
label = 'Hinf_control_sop';

% Complete the container options and read the tools of the PIE's dimension
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension (EXEC_TOOLS_SOP)        % MMP, 10/29/2026
% if nv==2,   error('PIETOOLS_Hinf_control_sop:dim2','Optimal control of 2D PIEs is currently not supported.');   end % MMP, 10/29/2026 (was)
if nv>=2,   error('PIETOOLS_Hinf_control_sop:dim2','Optimal control of PIEs in two or more variables is currently not supported.');   end % MMP, 10/29/2026
if ~(PIE.Tw==0) || ~(PIE.Tu==0)
    error('PIETOOLS_Hinf_control_sop:Tw','The LPI cannot be solved for systems with disturbances or inputs at the boundary.');
end
eppos = tl.eppos;
gfix = [];  if isfield(st.sop,'gam_fixed') && ~isempty(st.sop.gam_fixed) && st.sop.gam_fixed~=0,    gfix = st.sop.gam_fixed;    end

% The PIE operators as containers, and the identities on w and z
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);
Bw = op2copvar_sop(PIE.B1); Cz = op2copvar_sop(PIE.C1);  Dzw = op2copvar_sop(PIE.D11);
Bu = op2copvar_sop(PIE.B2); Dzu = op2copvar_sop(PIE.D12);
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
certified = any(strcmp(cl.status,{'optimal','inaccurate'}));
if isempty(gfix),   gam = cl.obj;
elseif certified,   gam = gfix;
else,               gam = Inf;
end

% The solved P and Z, and the gain K = Z P^{-1} on the container path
K = [];     P = [];     Z = [];     Psol = [];  Zsol = [];  Gc = [];    iinv = [];
try
    Psol = getsol_lpivar_sop(prog,aux.Pm);  Zsol = getsol_lpivar_sop(prog,aux.Zm);
    P = tl.toop(Psol);      Z = tl.toop(Zsol);
    [Gc,K,iinv] = getController_sop(Psol,Zsol);
catch ME
    if st.sop.verbose,  fprintf('  [%s] gain extraction: %s\n',label,ME.message);  end
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info_sop(prog,aux,st,label);
info.form = 'ctrl';     info.gam = gam;     info.gam_fixed = gfix;  info.nv = nv;
info.P = P;     info.Z = Z;     info.Psol = Psol;   info.Zsol = Zsol;   info.K = aux.Km;
info.G = Gc;    info.inv = iinv;
if st.sop.verbose
    fprintf('  [%s] %s: gamma %.8g (certified: %d); gain class %s\n',label,cl.status,gam,info.certified,class(K));
end

    function [prog,aux] = assemble(dw)
        % [PROG,AUX] = ASSEMBLE(DW) builds the synthesis LPI with the slack
        % degree raised by settings.sop.dw_P + DW.
        % INPUT  dw:    nonnegative integer, the raise of this pass.
        % OUTPUT prog:  'struct', the LPI program;
        %        aux:   'struct' with fields obj (the 'dpvar' gam), Pm, Zm
        %               ('cdopvar' P and Z), Km (the LPI operator), Nm (the
        %               slack), li (slack info of LPI_INEQ_SOP), form.
        % New program; record the equality rows when the dual is kept
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        if st.sop.keepdual,     prog.sopeq = {};    end
        % gam >= 0 as the objective, or the fixed gain
        if isempty(gfix)
            [prog,gamv] = lpidecvar(prog,'gam');
            prog = lpi_ineq(prog,gamv);
            prog = lpisetobj(prog,gamv);
        else
            gamv = gfix;
        end
        % P >= eppos I, and Z free from the state to u
        [prog,Pm] = tl.storage(prog,Tm,'out');
        if ~all(eppos==0),  Pm = Pm + tl.eye(1,eppos);  end
        [spu,dmu] = copvar_space_list(Bu,'in');
        [spx,dmx] = copvar_space_list(Bu,'out');
        [prog,Zm] = lpivar_cdopvar(prog,struct('out',dmu,'in',dmx),struct('out',{spu},'in',{spx}),tl.dom,tl.zdeg);
        % The LPI operator
        Dzu_r = on_registry_sop(Dzu,Zm);
        CZ = Cz*Pm + Dzu_r*Zm;      AZ = Am*Pm + Bu*Zm;
        Km = [-(gamv*Iz),   Dzw,         CZ*Tm';
               Dzw',        -(gamv*Iw),  Bw';
               Tm*CZ',      Bw,          AZ*Tm' + Tm*AZ'];
        % Km <= 0: the stock slack and an equality, or the slack sized from P
        if strcmp(st.sop.slack,'stock')
            [prog,Nm] = tl.slack(prog,Km);
            [prog,tag] = lpi_eq_sop(prog,Nm+Km,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
%           o = struct('like',Pm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin); % MMP, 10/29/2026 (was)
            o = struct('like',Pm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin,'cap',tl.capfor(Km)); % MMP, 10/29/2026
            [prog,Nm,li] = lpi_ineq_sop(prog,-Km,o);
        end
        obj = gamv;     if ~isempty(gfix),   obj = [];  end
        aux = struct('obj',obj,'Pm',Pm,'Zm',Zm,'Km',Km,'Nm',Nm,'li',li,'form','ctrl');
    end
% % % END body replaced by MMP, 10/28/2026
end
