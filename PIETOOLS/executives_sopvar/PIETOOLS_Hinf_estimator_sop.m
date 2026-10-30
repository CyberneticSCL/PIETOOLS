function [prog,L,gam,P,Z,info] = PIETOOLS_Hinf_estimator_sop(PIE,settings,gain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,L,GAM,P,Z,INFO] = PIETOOLS_HINF_ESTIMATOR_SOP(PIE,SETTINGS,GAIN)
% computes an observer gain L for the 1-D or 2-D PIE
%
%       T x_t + Tw w_t = A x + B1 w,   z = C1 x + D11 w,   y = C2 x + D21 w
%
% and a bound GAM on the L2 gain from w to the estimation error of z, by
% the LPI, posed on 'cdopvar' operators and solved as an SDP that minimizes
% gam,
%
%   P >= eppos I,  Z free,  PB = P*B1 + Z*D21,  PA = P*A + Z*C2,
%   [K11, -D11', K13; -D11, -gam I, C1; K13', C1', PA'*T + T'*PA] <= 0,
%   K11 = -gam I + Tw'*PB + PB'*Tw,  K13 = -PB'*T - Tw'*PA,
%
% with L = P^{-1}*Z. While the solver returns UNKNOWN or INFEASIBLE, the
% program is rebuilt with the degree of the negativity slack raised by
% one, up to settings.sop.max_raise times.
%
% INPUT
% - PIE:      'pie_struct' object of a 1-D or 2-D PIE, with T, Tw, A, B1, C1,
%             D11, C2 and D21 'opvar' (1-D) or 'opvar2d' (2-D) objects, or of
%             a container PIE in any number of variables, with these
%             operators 'copvar' objects (EXEC_PIE_SOP; Tw may be absent);
%             in two or more variables, Tw and B1 not both nonzero;
% - settings: (optional) 'struct' of LPI settings, as returned by lpisettings
%             or lpisettings_sop:
%               dd1, dd12, options1, options12, override1   1-D degrees of P;
%               ddZ     1-D degrees of Z;
%               settings_2d   2-D degrees of P (LF_deg, ...) and of Z
%                       (the largest entry of Zop_deg);
%               eppos, eppos2 (1-D), settings_2d.eppos (2-D)  margins of P;
%               settings_nd    degrees and margin of P and Z in 3 or more
%                              variables (EXEC_TOOLS_SOP);
%               sos_opts  solver options;
%               sop     container options (EXEC_SOP_SETTINGS): dw_P, the
%                       degree raise of the slack, max_raise, gam_fixed.
%             Default lpisettings('light');
% - gain:     (optional) positive double: solve the feasibility LPI at
%             gam = gain instead of minimizing gam (0 or [] minimizes).
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - L:        'opvar' (1-D) or 'opvar2d' (2-D) object, the observer gain from
%             y to the state ([] when the solver returned no solution, and
%             in 3 or more variables, where no gain extraction exists: the
%             LPI is solved and P, Z are returned);
% - gam:      double, the certified bound; with GAIN, GAIN when certified and
%             Inf otherwise; NaN when the solver returned no value;
% - P, Z:     'opvar' (1-D) or 'opvar2d' (2-D) objects, the solved P and Z;
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded' or 'unknown';
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               form ('est'), gam, gam_fixed, nv;
%               P, Z       as the outputs; Psol, Zsol the 'copvar' objects;
%               K          'cdopvar' object, the LPI operator;
%               G          'copvar' object, the gain L (1-D); inv  'struct' of
%                          the inversion of P (GETOBSERVER_SOP; 1-D);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          as returned by EXEC_INFO_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_Hinf_estimator_sop
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
% MMP, 10/28/2026: One executive for the 1-D and 2-D estimators of
%                10/08/2026: the body of this file (a call of
%                private/synth_build_sop) is replaced by the builder's 'est'
%                code (MMP, 10/08/2026, 10/09/2026), with the third argument
%                GAIN of PIETOOLS_Hinf_estimator_2D_sop; the programs
%                unchanged. Header to the header protocol.
% MMP, 10/29/2026: PIEs in 3 or more variables and container PIEs
%                (maintainer: the executives must support n variables, at
%                least 3-D). The PIE goes through EXEC_PIE_SOP (INITIALIZE
%                for an opvar PIE, as before); the 2-D Tw/B1 check holds for
%                every nv >= 2. In 3 or more variables the LPI is solved and
%                no gain is extracted (GETOBSERVER_SOP is 1-D, getObserver_2D
%                2-D): L = [] with a message.
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
% Default settings and fixed gain; check the PIE
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
if nargin>=3 && ~isempty(gain) && gain~=0,  settings.sop.gam_fixed = gain;  end
if ~isa(PIE,'pie_struct')
    error('PIETOOLS_Hinf_estimator_sop:pie','The PIE should be a ''pie_struct'' object.');
end
% PIE = initialize(PIE);                                                    % MMP, 10/29/2026 (was)
PIE = exec_pie_sop(PIE);                    % opvar or container PIE        % MMP, 10/29/2026
label = 'Hinf_estimator_sop';

% Complete the container options and read the tools of the PIE's dimension
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension (EXEC_TOOLS_SOP)        % MMP, 10/29/2026
Tw0 = (PIE.Tw==0);
% if nv==2 && ~Tw0 && ~(PIE.B1==0)                                          % MMP, 10/29/2026 (was)
if nv>=2 && ~Tw0 && ~(PIE.B1==0)                                            % MMP, 10/29/2026
    error('PIETOOLS_Hinf_estimator_sop:TwB','The PIE takes both the input w and its derivative dw/dt; LPI based estimator synthesis is currently not supported.');
end
gfix = [];  if isfield(st.sop,'gam_fixed') && ~isempty(st.sop.gam_fixed) && st.sop.gam_fixed~=0,    gfix = st.sop.gam_fixed;    end
eppos = tl.eppos;

% The PIE operators as containers, and the identities on w and z
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);
Bw = op2copvar_sop(PIE.B1); Cz = op2copvar_sop(PIE.C1);  Dzw = op2copvar_sop(PIE.D11);
Cy = op2copvar_sop(PIE.C2); Dyw = op2copvar_sop(PIE.D21);
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);
Twm = [];   if ~Tw0,    Twm = op2copvar_sop(PIE.Tw);    end

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
certified = any(strcmp(cl.status,{'optimal','inaccurate'}));
if isempty(gfix),   gam = cl.obj;
elseif certified,   gam = gfix;
else,               gam = Inf;
end

% The solved P and Z, and the gain L = P^{-1} Z (container path in 1-D)
L = [];     P = [];     Z = [];     Psol = [];  Zsol = [];  Gc = [];    iinv = [];
try
    Psol = getsol_lpivar_sop(prog,aux.Pm);  Zsol = getsol_lpivar_sop(prog,aux.Zm);
    P = tl.toop(Psol);      Z = tl.toop(Zsol);
    if nv==1,   [Gc,L,iinv] = getObserver_sop(Psol,Zsol);
%   else,       L = getObserver_2D(P,Z);                                    % MMP, 10/29/2026 (was)
    elseif nv==2,   L = getObserver_2D(P,Z);                                % MMP, 10/29/2026
    elseif st.sop.verbose                   % no N-D gain extraction        % MMP, 10/29/2026
        fprintf('  [%s] gain extraction: none in %d variables; P and Z returned.\n',label,nv); % MMP, 10/29/2026
    end
catch ME
    if st.sop.verbose,  fprintf('  [%s] gain extraction: %s\n',label,ME.message);  end
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info_sop(prog,aux,st,label);
info.form = 'est';      info.gam = gam;     info.gam_fixed = gfix;  info.nv = nv;
info.P = P;     info.Z = Z;     info.Psol = Psol;   info.Zsol = Zsol;   info.K = aux.Km;
info.G = Gc;    info.inv = iinv;
if st.sop.verbose
    fprintf('  [%s] %s: gamma %.8g (certified: %d); gain class %s\n',label,cl.status,gam,info.certified,class(L));
end

    function [prog,aux] = assemble(dw)
        % [PROG,AUX] = ASSEMBLE(DW) builds the estimator LPI with the slack
        % degree raised by settings.sop.dw_P + DW.
        % INPUT  dw:    nonnegative integer, the raise of this pass.
        % OUTPUT prog:  'struct', the LPI program;
        %        aux:   'struct' with fields obj (the 'dpvar' gam, [] at a
        %               fixed gain), Pm, Zm ('cdopvar' P and Z), Km (the LPI
        %               operator), Nm (the slack), li (slack info of
        %               LPI_INEQ_SOP), form.
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
        % P >= eppos I, and Z free from y to the state
        [prog,Pm] = tl.storage(prog,Tm,'out');
        if ~all(eppos==0),  Pm = Pm + tl.eye(1,eppos);  end
        [spx,dmx] = copvar_space_list(Cy,'in');
        [spy,dmy] = copvar_space_list(Cy,'out');
        [prog,Zm] = lpivar_cdopvar(prog,struct('out',dmx,'in',dmy),struct('out',{spx},'in',{spy}),tl.dom,tl.zdeg);
        % The LPI operator, with the terms of a boundary disturbance
        Dyw_r = on_registry_sop(Dyw,Zm);
        PB = Pm*Bw + Zm*Dyw_r;      PA = Pm*Am + Zm*Cy;
        K11 = -(gamv*Iw);
        K13 = -PB'*Tm;              K31 = -Tm'*PB;
        if ~Tw0
            X   = Twm'*PB;
            K11 = on_registry_sop(K11,X) + X + PB'*Twm;
            K13 = K13 - Twm'*PA;    K31 = K31 - PA'*Twm;
        end
        Km = [K11,      -Dzw',       K13;
              -Dzw,     -(gamv*Iz),  Cz;
              K31,      Cz',         PA'*Tm + Tm'*PA];
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
        aux = struct('obj',obj,'Pm',Pm,'Zm',Zm,'Km',Km,'Nm',Nm,'li',li,'form','est');
    end
% % % END body replaced by MMP, 10/28/2026
end
