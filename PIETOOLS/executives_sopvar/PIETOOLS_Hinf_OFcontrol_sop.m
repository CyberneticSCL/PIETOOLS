function [prog,Kc,gam,info] = PIETOOLS_Hinf_OFcontrol_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,KC,GAM,INFO] = PIETOOLS_HINF_OFCONTROL_SOP(PIE,SETTINGS)
% computes a dynamic output feedback controller for the PIE
%
%       T x_t = A x + B1 w + B2 u,
%       z = C1 x + D11 w + D12 u,       y = C2 x + D21 w      (D22 = 0),
%
% in the form of the PIE
%
%       TK xK_t = AK xK + BK y,         u = CK xK + DK y,
%
% and a bound GAM on the L2 gain from w to z of the closed loop, by an LPI
% from Jungbae Chun: the primal KYP lemma (Theorem 11(a) of Shivakumar
% et al.) applied to the closed loop with state [x; xK], after a congruence
% transformation and the change of variables An, Bn, Cn, Dn of the free
% operators Q, U. The LPI minimizes gam over the positive operators R, V and
% the free operators Q, U, An, Bn, Cn, Dn such that
%       TK = T'(QU - T), with
%     T'*Q = R,   T'T*U = V,   [R T'T; T'T V] >= 0,   Phi <= 0,
%     Phi = [ He(Q'A+Bn C2)        An+(A+B2 Dn C2)'T    Q'B1+Bn D21         C1'+C2'Dn'D12'
%             .                    He(T'(A U+B2 Cn))    T'(B1+B2 Dn D21)    U'C1'+Cn'D12'
%             .                    .                    -gam I              D11'+D21'Dn'D12'
%             .                    .                    .                   -gam I ],
%     and the controller is recovered by:
%     DK = Dn,  CK = Dn C2 U - Cn,  BK = Bn - Q'B2 Dn,
%     AK = BK C2 U + Q'(A U + B2 Cn) - An.
% No PI operator is inverted, and TK is not fixed to T: it is recovered from
% the solution. Here, He(X) = X + X' and the lower triangle of Phi is its
% adjoint. The operators are 'copvar' containers, so the PIE may be in any
% number of variables (1-D, 2-D, ...).
% While the solver returns UNKNOWN or INFEASIBLE, the program is rebuilt with
% the degree of the negativity slack raised by one, up to settings.sop.max_raise
% times.
%
% INPUT
% - PIE:      'pie_struct' object, with T, A, B1, B2, C1, D11, D12, C2 and
%             D21 'opvar' (or 'opvar2d') objects, or 'copvar' objects (a
%             container PIE, EXEC_PIE_SOP), and D22 zero; Tw and Tu zero or
%             absent (no disturbance or input at the boundary), and no
%             algebraic constraint K;
% - settings: (optional) 'struct' of LPI settings, as returned by lpisettings
%             or lpisettings_sop:
%               dd1, dd12, options1, options12, override1   degrees of R, V;
%               eppos, eppos2   margins of R and V; the free operators Q and U
%                       rescale R and V, and with the default margins of
%                       lpisettings (1e-4) the solver stalls (status 'unknown')
%                       on the examples tried, so eppos = eppos2 = 1e-7 is
%                       advised;
%               sos_opts  solver options;
%               sop     container options (EXEC_SOP_SETTINGS): dw_Q, the
%                       degree raise of the slack, max_raise,
%                       gam_fixed (a positive double: solve the
%                       feasibility LPI at gam = gam_fixed), and
%                       coupling_margin (default 0: the coupling LPI is
%                       required to hold with R and V reduced by this
%                       multiple of the identity).
%             Default lpisettings('light').
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - Kc:       'struct' with fields TK, AK, BK, CK, DK, the 'copvar' objects
%             of the controller PIE above, and G, a 'struct' with the same
%             five fields; [] when the solver returned no solution or the
%             controller could not be recovered;
% - gam:      double, the certified closed-loop bound; with gam_fixed,
%             gam_fixed when certified and Inf otherwise; NaN when the
%             solver returned no value;
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded', 'unknown' or 'residual' (cuADMM,
%                          which decides no feasibility: never certified);
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%                gam, gam_fixed ([]), nv;
%               sol        'struct' of the solved 'copvar' operators Q, U, An,
%                          Bn, Cn, Dn, R and V;
%               K          'cdopvar' object, the LPI operator Phi;
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          as returned by EXEC_INFO_SOP.
%
% NOTES
% - The closed loop of the PIE and the controller is the PIE returned by
%   CLOSEDLOOPPIE(PIE,KC,'dynamic'), to which the analysis executives apply; its L2
%   gain should be tested with PIETOOLS_Hinf_gain_sop in form 'Q', as the controller carries factors of T.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_Hinf_OFcontrol_sop
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
% Initial coding: DB 10/31/2026, based on formulation and code of Jungbae Chun
% DB 11/01/2026 - Removed the form 'inv' (TK = T, inverts PI operators); only
%                 the inversion-free form remains, with TK recovered
% MMP 11/05/2026 - Kc holds the 'copvar' operators: the conversion by
%                 tl.toop needs the opvar classes in 1-D and 2-D, and without
%                 them the recovery failed inside the try/catch, so Kc was [].
%                 gam >= 0 by LPI_INEQ_SOP instead of the stock LPI_INEQ. No
%                 OP2COPVAR_SOP on the PIE operators, which EXEC_PIE_SOP has
%                 already converted (it returned its input).
% MMP 11/07/2026 - Help text: status 'residual', which LPI_CLASSIFY_SOP
%                 returns for a cuADMM exit since 11/04/2026. Code unchanged.
% MMP 11/08/2026 - Default settings lpisettings_sop('light') (was the stock
%                 lpisettings('light'), whose folder leaves the path): the
%                 same storage degrees and margins, with the container
%                 options of the light tier (maintainer, 11/08/2026).

% Default settings; check the PIE
% if nargin<2 || isempty(settings),   settings = lpisettings('light');    end % MMP, 11/08/2026 (was)
if nargin<2 || isempty(settings),  settings = lpisettings_sop('light'); end % MMP, 11/08/2026

if ~isa(PIE,'pie_struct')
    error('PIETOOLS_Hinf_OFcontrol_sop:pie','The PIE should be a ''pie_struct'' object.');
end
PIE = exec_pie_sop(PIE);                    % opvar or container PIE
label = 'Hinf_OFcontrol_sop';

% Complete the container options and read the tools of the PIE's dimension
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension (EXEC_TOOLS_SOP)
if ~(PIE.Tw==0) || ~(PIE.Tu==0)
    error('PIETOOLS_Hinf_OFcontrol_sop:Tw','The LPI cannot be solved for systems with disturbances or inputs at the boundary.');
end
if ~(PIE.D22==0)
    error('PIETOOLS_Hinf_OFcontrol_sop:D22','The LPI is posed for D22 = 0 (the controller depends on y only through the state).');
end
if isa(PIE.K,'copvar') && size(PIE.K,1)>0
    error('PIETOOLS_Hinf_OFcontrol_sop:K','PIEs with an algebraic constraint K are not supported.');
end
eppos = tl.eppos;
cmar = 0;   if isfield(st.sop,'coupling_margin') && ~isempty(st.sop.coupling_margin),  cmar = st.sop.coupling_margin;   end
gfix = [];  if isfield(st.sop,'gam_fixed') && ~isempty(st.sop.gam_fixed) && st.sop.gam_fixed~=0,    gfix = st.sop.gam_fixed;    end

% The PIE operators as containers, and the identities on w and z
% EXEC_PIE_SOP has converted them: no OP2COPVAR_SOP                         % MMP, 11/05/2026
% Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);                    % MMP, 11/05/2026 (was)
% Bw = op2copvar_sop(PIE.B1); Bu = op2copvar_sop(PIE.B2);                   % MMP, 11/05/2026 (was)
% Cz = op2copvar_sop(PIE.C1); Dzw = op2copvar_sop(PIE.D11);   Dzu = op2copvar_sop(PIE.D12); % MMP, 11/05/2026 (was)
% Cy = op2copvar_sop(PIE.C2); Dyw = op2copvar_sop(PIE.D21);                 % MMP, 11/05/2026 (was)
Tm = PIE.T;     Am = PIE.A;                                                 % MMP, 11/05/2026
Bw = PIE.B1;    Bu = PIE.B2;                                                % MMP, 11/05/2026
Cz = PIE.C1;    Dzw = PIE.D11;      Dzu = PIE.D12;                          % MMP, 11/05/2026
Cy = PIE.C2;    Dyw = PIE.D21;                                              % MMP, 11/05/2026
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);

% The spaces of the free operators, as lists of the variables of each space
% (Q in -> out, U, An in -> in, Bn y -> in, Cn in -> u, Dn y -> u)
[sp_in,dm_in]   = copvar_space_list(Tm,'in');
[sp_out,dm_out] = copvar_space_list(Tm,'out');
[sp_u,dm_u]     = copvar_space_list(Bu,'in');       % the actuator input u
[sp_y,dm_y]     = copvar_space_list(Cy,'out');      % the observed output y

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
certified = any(strcmp(cl.status,{'optimal','inaccurate'}));
if isempty(gfix),   gam = cl.obj;
elseif certified,   gam = gfix;
else,               gam = Inf;
end

% The solved operators, and the controller recovered from them
Kc = [];    sol = struct();
try
    for f = {'Q','U','An','Bn','Cn','Dn','R','V'}
        sol.(f{1}) = getsol_lpivar_sop(prog,aux.(f{1}));
    end
    % The controller from the change of variables, with no inverse
    DK = sol.Dn;
    CK = DK*Cy*sol.U - sol.Cn;
    BK = sol.Bn - sol.Q'*Bu*DK;
    AK = BK*Cy*sol.U + sol.Q'*(Am*sol.U + Bu*sol.Cn) - sol.An;
    TK = Tm'*(sol.Q*sol.U - Tm);
%   Kc = struct('TK',tl.toop(TK),'AK',tl.toop(AK),'BK',tl.toop(BK),'CK',tl.toop(CK),'DK',tl.toop(DK)); % MMP, 11/05/2026 (was)
    % The containers: tl.toop needs the opvar classes in 1-D and 2-D        % MMP, 11/05/2026
    Kc = struct('TK',TK,'AK',AK,'BK',BK,'CK',CK,'DK',DK);                   % MMP, 11/05/2026
    Kc.G = struct('TK',TK,'AK',AK,'BK',BK,'CK',CK,'DK',DK);
catch ME
    if st.sop.verbose,  fprintf('  [%s] controller recovery: %s\n',label,ME.message);  end
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info_sop(prog,aux,st,label);
info.form = 'OF';   info.gam = gam;     info.gam_fixed = gfix;  info.nv = nv;
info.sol = sol;
info.K = aux.Km;
if st.sop.verbose
    fprintf('  [%s] %s: gamma %.8g (certified: %d); controller %s\n',label,cl.status,gam,info.certified,class(Kc));
end

    function [prog,aux] = assemble(dw)
        % [PROG,AUX] = ASSEMBLE(DW) builds the synthesis LPI with the slack
        % degree raised by settings.sop.dw_Q + DW.
        % INPUT  dw:    nonnegative integer, the raise of this pass.
        % OUTPUT prog:  'struct', the LPI program;
        %        aux:   'struct' with fields obj (the 'dpvar' gam), Q, U, An,
        %               Bn, Cn, Dn, R, V ('cdopvar' objects), Km (the LPI
        %               operator), Nm (the slack), li (slack info of
        %               LPI_INEQ_SOP).
        % New program; record the equality rows when the dual is kept
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        if st.sop.keepdual,     prog.sopeq = {};    end
        % gam >= 0 as the objective, or the fixed gain
        if isempty(gfix)
            [prog,gamv] = lpidecvar(prog,'gam');
%           prog = lpi_ineq(prog,gamv);                                     % MMP, 11/05/2026 (was)
            prog = lpi_ineq_sop(prog,gamv);                                 % MMP, 11/05/2026
            prog = lpisetobj(prog,gamv);
        else
            gamv = gfix;
        end
        % The positive operators R and V, both on the input space of T
        % (V = T'T U)
        [prog,Rm] = tl.storage(prog,Tm,'in');
        [prog,Vm] = tl.storage(prog,Tm,'in');
        if ~all(eppos==0)
            Rm = Rm + tl.eye(2,eppos);      Vm = Vm + tl.eye(2,eppos);
        end
        % The free operators
        dQ = get_lpivar_degs_sop(Rm);       dU = get_lpivar_degs_sop(Vm);
        % Q: in -> out, U and An: in -> in, Bn: y -> in, Cn: in -> u, Dn: y -> u
        sUo = {sp_in,dm_in,sp_in,dm_in};    sCn = {sp_u,dm_u,sp_in,dm_in};
        [prog,Qm]  = lpivar_cdopvar(prog,struct('out',dm_out,'in',dm_in),struct('out',{sp_out},'in',{sp_in}),tl.dom,dQ);
        [prog,Um]  = lpivar_cdopvar(prog,struct('out',sUo{2},'in',sUo{4}),struct('out',sUo(1),'in',sUo(3)),tl.dom,dU);
        [prog,Anm] = lpivar_cdopvar(prog,struct('out',sUo{2},'in',sUo{4}),struct('out',sUo(1),'in',sUo(3)),tl.dom,dU);
        [prog,Bnm] = lpivar_cdopvar(prog,struct('out',dm_in,'in',dm_y),struct('out',{sp_in},'in',{sp_y}),tl.dom,dQ);
        [prog,Cnm] = lpivar_cdopvar(prog,struct('out',sCn{2},'in',sCn{4}),struct('out',sCn(1),'in',sCn(3)),tl.dom,dU);
        [prog,Dnm] = lpivar_cdopvar(prog,struct('out',dm_u,'in',dm_y),struct('out',{sp_u},'in',{sp_y}),tl.dom,dQ);
        % The equalities T'Q = R and T'T U = V
        Tc = Tm'*Tm;
        prog = lpi_eq_sop(prog,Tm'*Qm-Rm);
        prog = lpi_eq_sop(prog,Tc*Um-Vm);
        % The LPI operator Phi, and the coupling [R T'T; T'T V] >= 0
        Dzu_r = on_registry_sop(Dzu,Cnm);
        Phi11 = Qm'*Am + Bnm*Cy;                Phi11 = Phi11 + Phi11';
        Phi13 = Qm'*Bw + Bnm*Dyw;
        Phi14 = Cz' + Cy'*Dnm'*Dzu_r';
        Phi24 = Um'*Cz' + Cnm'*Dzu_r';
        Phi34 = Dzw' + Dyw'*Dnm'*Dzu_r';
        Phi12 = Anm + (Am + Bu*Dnm*Cy)'*Tm;
        Phi22 = Tm'*(Am*Um + Bu*Cnm);       Phi22 = Phi22 + Phi22';
        Phi23 = Tm'*(Bw + Bu*Dnm*Dyw);
        Km = [Phi11,    Phi12,      Phi13,          Phi14;
              Phi12',   Phi22,      Phi23,          Phi24;
              Phi13',   Phi23',     -(gamv*Iw),     Phi34;
              Phi14',   Phi24',     Phi34',         -(gamv*Iz)];
        Cp = [Rm-tl.eye(2,cmar), Tc; Tc', Vm-tl.eye(2,cmar)];   % margin: [R T'T; T'T V] >= cmar*I
        % Km <= 0 and Cp >= 0: the stock slack and an equality, or the slack
        % sized from R
        if strcmp(st.sop.slack,'stock')
            [prog,Nm] = tl.slack(prog,Km);
            [prog,tag] = lpi_eq_sop(prog,Nm+Km,'symmetric');
            [prog,Nc] = tl.slack(prog,-Cp);
            prog = lpi_eq_sop(prog,Nc-Cp,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
            o = struct('like',Rm,'dw',st.sop.dw_Q+dw,'margin',st.sop.margin,'cap',tl.capfor(Km));
            [prog,Nm,li] = lpi_ineq_sop(prog,-Km,o);
            o = struct('like',Rm,'dw',st.sop.dw_Q+dw,'margin',st.sop.margin,'cap',tl.capfor(Cp));
            prog = lpi_ineq_sop(prog,Cp,o);
        end
        obj = gamv;     if ~isempty(gfix),   obj = [];  end
        aux = struct('obj',obj,'Q',Qm,'U',Um,'An',Anm,'Bn',Bnm,'Cn',Cnm,'Dn',Dnm,'R',Rm,'V',Vm, ...
                     'Km',Km,'Nm',Nm,'li',li,'form','OF');
    end
end
