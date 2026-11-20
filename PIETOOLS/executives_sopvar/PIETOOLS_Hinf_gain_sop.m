function [prog,R,gam,info] = PIETOOLS_Hinf_gain_sop(PIE,settings,form,gain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,R,GAM,INFO] = PIETOOLS_HINF_GAIN_SOP(PIE,SETTINGS,FORM,GAIN)
% computes an upper bound GAM on the L2 gain from w to z of the PIE
%
%       T x_t + Tw w_t = A x + B1 w,      z = C1 x + D11 w,
%
% by a KYP LPI on 'cdopvar' operators, solved as an SDP that minimizes gam.
% FORM selects the LPI; R is a positive operator (in 'Q' the generalized
% storage T'*Q), Q a free operator:
%
%   'Q'   R >= 0, T'*Q = R,
%         [-gam I, D11', B1'*Q; D11, -gam I, C1; Q'*B1, C1', A'*Q + Q'*A] <= 0;
%   'Qd'  R >= 0, T*Q = R,
%         [-gam I, D11, C1*Q; D11', -gam I, B1'; Q'*C1', B1, Q'*A' + A*Q] <= 0;
%   'P'   R >= eppos I,
%         [K11, D11', K13; D11, -gam I, C1; K13', C1', (R*A)'*T + T'*(R*A)] <= 0,
%         K11 = -gam I + Tw'*R*B1 + B1'*R*Tw,  K13 = (R*B1)'*T + Tw'*R*A;
%   'Pd'  R >= eppos I,
%         [-gam I, D11, C1*R*T'; D11', -gam I, B1'; T*R*C1', B1, T*R*A' + A*R*T'] <= 0.
%
% 'Q' is posed as 'P' when Tw is nonzero; 'Qd' and 'Pd' require Tw = 0; in
% two or more variables, 'Qd' is posed as 'Pd'. In 'Pd', settings without
% eppos give R the margin 1e-4 on R^n and 1e-6 on the space in every
% variable, 0 on the other spaces. While the solver returns UNKNOWN or
% INFEASIBLE, the program is rebuilt with the degree of the negativity
% slack raised by one, up to settings.sop.max_raise times.
%
% If PIE.K is nonempty, the PIE state satisfies K x = 0 (e.g. periodic
% boundary conditions): in 'Q' and 'P' the block of the KYP operator on the
% state carries + K'*Z + Z'*K, with Z a free operator, which vanishes on
% ker K; 'Qd' and 'Pd' are posed as 'Q' and 'P'.
%
% INPUT
% - PIE:      'pie_struct' object of a PIE in any number of spatial
%             variables, initialized when its T is not yet a 'copvar'
%             object; T, Tw, A, B1, C1 and D11 'copvar' objects (Tw may be
%             absent), and K a 'copvar' object (the constraint K x = 0 on
%             the state) or empty; its other operators are not used; in two
%             or more variables, Tw and B1 not both nonzero;
% - settings: (optional) 'struct' of LPI settings:
%               lpi.pos   the storage row of R, in every dimension;
%                         settings_2d (optional) the 16-basis 2-D degrees,
%                         in place of lpi.pos in two variables;
%               eppos     the margin of R in the 'P' forms, in every
%                         dimension;
%               sos_opts  solver options;
%               sop       the options of the executives, among them dw_Q
%                         (the degree raise of the slack in the 'Q' forms),
%                         dw_P (in the 'P' forms) and max_raise.
%             Default: the 'light' tier;
% - form:     (optional) char, 'Q' (default), 'Qd', 'P' or 'Pd';
% - gain:     (optional) positive double: solve the feasibility LPI at
%             gam = gain instead of minimizing gam (0 or [] minimizes).
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - R:        'copvar' object, the solved positive operator R ([] when the
%             solver returned no solution);
% - gam:      double, the certified bound; with GAIN, GAIN when certified
%             and Inf otherwise; NaN when the solver returned no value;
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded', 'unknown' or 'residual' (cuADMM,
%                          which decides no feasibility: never certified);
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               form, gam, gam_fixed, nv   the LPI posed, the bound, GAIN
%                          ([] if none), the dimension;
%               storage    'copvar' object, the solved R;
%               K, N       'cdopvar' objects, the KYP operator and the slack;
%               gain_lb    double, in 1-D with settings.sop.witness the
%                          largest gain of the Legendre discretization of
%                          the PIE over the sampled frequencies (a lower
%                          bound), else NaN; omega its frequency; gap
%                          gam/gain_lb - 1;
%               dual_alignment  double in [0,1], alignment of the dual
%                          kernel with the worst input of that
%                          discretization (NaN if not computed);
%               witness    'struct' of the discretization ([] if none);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          the last raise and the loop history, the slack
%                          mode, its degrees and terms, the size of the SDP,
%                          and the dual kernel of the negativity constraint
%                          with its 1-norm.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_Hinf_gain_sop
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

% Default settings, form and fixed gain
if nargin<2 || isempty(settings),  settings = lpisettings_sop('light'); end
if nargin<3 || isempty(form),       form = 'Q';     end
if nargin>=4 && ~isempty(gain) && gain~=0,  settings.sop.gam_fixed = gain;  end
if ~isa(PIE,'pie_struct')
    error('PIETOOLS_Hinf_gain_sop:pie','The PIE should be a ''pie_struct'' object.');
end
PIE = exec_pie_sop(PIE);                    % every operator a 'copvar' or 0
label = ['Hinf_gain_sop ' form];

% Complete the options of the executives and read the tools of the
% dimension of the PIE (storage of R, margins)
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension
Tw0 = (PIE.Tw==0);

% Check the form against the PIE: boundary disturbances and the dimension
if nv>=2 && ~Tw0 && ~(PIE.B1==0)
    error('PIETOOLS_Hinf_gain_sop:TwB','The PIE takes both the input w and its derivative dw/dt; not supported in two or more variables.');
end
% The constraint K x = 0 on the state, if any (empty for most PIEs). The
% dual forms have no term for it: posed as the primal ones
Kcon = PIE.K;   hasK = isa(Kcon,'copvar') && size(Kcon,1)>0;
if hasK && strcmp(form,'Qd'),   form = 'Q';     end
if hasK && strcmp(form,'Pd'),   form = 'P';     end
switch form
    case 'Q'
        if ~Tw0,    form = 'P';     label = [label ' -> P (Tw ~= 0)'];   end
    case {'Qd','Pd'}
        if ~Tw0
            error('PIETOOLS_Hinf_gain_sop:Tw','The dual LPIs cannot be solved for systems with disturbances at the boundary.');
        end
        if strcmp(form,'Qd') && nv>=2,  form = 'Pd';    end     % no N-D Qd form
    case 'P'
    otherwise
        error('PIETOOLS_Hinf_gain_sop:form','FORM should be ''Q'', ''Qd'', ''P'' or ''Pd''.');
end
coercive = any(strcmp(form,{'P','Pd'}));

% Degree raise of the slack and the fixed gain, if any
dwbase = st.sop.dw_Q;   if coercive,    dwbase = st.sop.dw_P;   end
gfix = [];  if isfield(st.sop,'gam_fixed') && ~isempty(st.sop.gam_fixed) && st.sop.gam_fixed~=0,    gfix = st.sop.gam_fixed;    end

% The PIE operators, 'copvar' objects, and the identities on w and z
Tm = PIE.T;     Am = PIE.A;
Bw = PIE.B1;    Cz = PIE.C1;    Dzw = PIE.D11;
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);
Twm = [];   if ~Tw0,    Twm = PIE.Tw;    end

% Margin of R: eppos I in the 'P' forms in every dimension (in 'Pd' without
% settings eppos: 1e-4 on R^n, 1e-6 on the space in every variable, 0
% between, over the 2^nv variable subsets)
e = [];
if coercive
    e = tl.eppos;
    if strcmp(form,'Pd') && tl.eppos_default
        e = zeros(2^nv,1);  e(1) = 1e-4;    e(end) = 1e-6;  % [1e-4;0;0;1e-6] in 2-D
    end
end

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
certified = any(strcmp(cl.status,{'optimal','inaccurate'}));
if isempty(gfix),   gam = cl.obj;
elseif certified,   gam = gfix;
else,               gam = Inf;
end

% The solved R, as a 'copvar'
Rsol = [];  R = [];
try
    Rsol = getsol_lpivar_sop(prog,aux.Rm);
    R = Rsol;
catch
    R = Rsol;
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info_sop(prog,aux,st,label);
info.form = form;   info.gam = gam;  info.gam_fixed = gfix;   info.nv = nv;
info.storage = Rsol;    info.K = aux.Km;    info.N = aux.Nm;

% Numerical witness (1-D): the gain of the discretized PIE, a lower bound
info.witness = [];  info.gain_lb = NaN;     info.omega = NaN;   info.gap = NaN;     info.dual_alignment = NaN;
if st.sop.witness && nv==1
    % In 'Q' and 'P', with the dual kernel and its block of w: the
    % alignment of the dual with the worst input
    if any(strcmp(form,{'Q','P'})),     Wt = pie_witness_sop(PIE,'gain',st.sop,info.dual,aux.wblk);
    else,                               Wt = pie_witness_sop(PIE,'gain',st.sop);
    end
    info.witness = Wt;
    if Wt.ok
        info.gain_lb = Wt.gain;     info.omega = Wt.omega;  info.dual_alignment = Wt.dual_alignment;
        if ~isnan(gam) && Wt.gain>0,    info.gap = gam/Wt.gain - 1;     end
    end
end
if st.sop.verbose
    fprintf('  [%s] %s: gamma %.8g (certified upper bound: %d); numerical lower bound %.8g at omega %.4g; gap %.2e; dual alignment %.3f\n', ...
            label,cl.status,gam,info.certified,info.gain_lb,info.omega,info.gap,info.dual_alignment);
end

    function [prog,aux] = assemble(dw)
        % [PROG,AUX] = ASSEMBLE(DW) builds the KYP LPI of FORM with the
        % slack degree raised by DWBASE + DW.
        % INPUT  dw:    nonnegative integer, the raise of this pass.
        % OUTPUT prog:  'struct', the LPI program;
        %        aux:   'struct' with fields obj (the 'dpvar' gam, [] at a
        %               fixed gain), Rm ('cdopvar' R), Km (the KYP
        %               operator), Nm (the slack), li (the info of the
        %               slack), wblk (the block of K that w enters), form.
        % New program; record the equality rows when the dual is kept
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);
        if st.sop.keepdual,     prog.sopeq = {};    end
        % gam >= 0 as the objective, or the fixed gain
        if isempty(gfix)
            [prog,gamv] = lpidecvar(prog,'gam');
            prog = lpi_ineq_sop(prog,gamv);
            prog = lpisetobj(prog,gamv);
        else
            gamv = gfix;
        end
        % The positive operator R, with its margin
        [prog,Rm] = tl.storage(prog,Tm,'out');
        if ~isempty(e) && ~all(e==0),   Rm = Rm + tl.eye(1,e);     end
        % The KYP operator of FORM
        switch form
            case 'Q'
                [sp,dm] = copvar_space_list(Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,dm,sp,tl.dom,get_lpivar_degs_sop(Rm));
                prog = lpi_eq_sop(prog,Tm'*Qm-Rm);
                % The state block, with K'*Z + Z'*K if K is given
                X33 = Am'*Qm + Qm'*Am;
                if hasK,    [prog,X33] = exec_kmult_sop(prog,X33,Kcon,Tm,tl.dom,Rm);    end
                Km = [-(gamv*Iw),   Dzw',       Bw'*Qm;
                       Dzw,         -(gamv*Iz), Cz;
                       Qm'*Bw,      Cz',        X33];
                wblk = 1;
            case 'Qd'
                % R = T*Q; the dual KYP operator
                [sp,dm] = copvar_space_list(Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,dm,sp,tl.dom,get_lpivar_degs_sop(Rm));
                prog = lpi_eq_sop(prog,Tm*Qm-Rm);
                Km = [-(gamv*Iz),   Dzw,        Cz*Qm;
                       Dzw',        -(gamv*Iw), Bw';
                       Qm'*Cz',     Bw,         Qm'*Am' + Am*Qm];
                wblk = 2;
            case 'P'
                % R in place of Q, with the terms of Tw
                PB = Rm*Bw;     PA = Rm*Am;
                K11 = -(gamv*Iw);   K13 = PB'*Tm;   K31 = Tm'*PB;
                if ~Tw0         % the terms of the boundary disturbance
                    X = Twm'*PB;
                    K11 = on_registry_sop(K11,X) + X + PB'*Twm;
                    K13 = K13 + Twm'*PA;    K31 = K31 + PA'*Twm;
                end
                % The state block, with K'*Z + Z'*K if K is given
                X33 = PA'*Tm + Tm'*PA;
                if hasK,    [prog,X33] = exec_kmult_sop(prog,X33,Kcon,Tm,tl.dom,Rm);    end
                Km = [K11,      Dzw',       K13;
                      Dzw,      -(gamv*Iz), Cz;
                      K31,      Cz',        X33];
                wblk = 1;
            case 'Pd'
                CP = Cz*Rm;     AP = Am*Rm;
                Km = [-(gamv*Iz),   Dzw,        CP*Tm';
                       Dzw',        -(gamv*Iw), Bw';
                       Tm*CP',      Bw,         Tm*AP' + AP*Tm'];
                wblk = 2;
        end
        % K <= 0: the 'stock' slack and an equality, or the slack sized from R
        if strcmp(st.sop.slack,'stock')
            [prog,Nm] = tl.slack(prog,Km);
            [prog,tag] = lpi_eq_sop(prog,Nm+Km,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
            o = struct('like',Rm,'dw',dwbase+dw,'margin',st.sop.margin,'cap',tl.capfor(Km));
            [prog,Nm,li] = lpi_ineq_sop(prog,-Km,o);
        end
        % The objective gam, none at a fixed gain
        obj = gamv;     if ~isempty(gfix),   obj = [];  end
        aux = struct('obj',obj,'Rm',Rm,'Km',Km,'Nm',Nm,'li',li,'wblk',wblk,'form',form);
    end
end
