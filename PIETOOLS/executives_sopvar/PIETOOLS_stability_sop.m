function [prog,P,info] = PIETOOLS_stability_sop(PIE,settings,form)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,INFO] = PIETOOLS_STABILITY_SOP(PIE,SETTINGS,FORM) tests stability
% of the PIE
%
%       T x_t = A x
%
% by searching for a Lyapunov certificate, posed as an LPI on 'cdopvar'
% operators and solved as an SDP. FORM selects the LPI:
%
%   'P'   P >= eppos I,          T'*P*A + A'*P*T + 2*epneg*T'*P*T <= 0;
%   'Pd'  P >= eppos I,          T*P*A' + A*P*T' + 2*epneg*T*P*T' <= 0;
%   'Q'   P >= eppos2*T'*T,  T'*Q = P,   A'*Q + Q'*A + 2*epneg*P <= 0;
%         (posed directly: T'*Q = Q'*T, P = T'*Q, no storage)               % MMP, 11/18/2026 (was); MMP, 11/17/2026
%         (posed directly: T'*Q = Q'*T, P = T'*Q the storage, Q at          % MMP, 11/18/2026
%         settings.sop.deg_Q; 'P' is the case Q = P*T)                      % MMP, 11/18/2026
%   'Qd'  P >= eppos2*T*T',  T*Q = P,    A*Q + Q'*A' + 2*epneg*P <= 0;
%   'FP'  P >= 0,                T'*P*A + A'*P*T + 2*epneg*T'*T <= 0;
%   'FQ'  P >= 0,            T'*Q = P,   A'*Q + Q'*A + 2*epneg*T'*T <= 0;
%
% with P a positive operator, Q a free operator, eppos I the identity
% scaled by the margin of each space (eppos2 the margin on L2), and epneg
% the decay rate, in every dimension: a solution of 'P' gives
% ||T x(t)|| <= M exp(-epneg t) ||T x(0)|| for some M > 0. The forms 'FP', 'FQ' test finite
% energy, int_0^inf ||T x(t)||^2 dt <= C ||T x(0)||^2 for some C > 0, and
% need epneg > 0. 'Q' is posed in every dimension, eppos2 the least L2      % MMP, 11/17/2026
% margin; 'Qd' and 'FQ' exist in 1-D only, and for a PIE in two or more     % MMP, 11/17/2026
% variables are posed as 'Pd' and 'FP'. While the solver                    % MMP, 11/17/2026
% returns UNKNOWN or INFEASIBLE, the program is rebuilt with the degree of
% the negativity slack raised by one, up to settings.sop.max_raise_feas (its
% entry for the dimension of the PIE, the last entry for a higher one; or
% max_raise) times.
%
% If PIE.K is nonempty, the PIE state satisfies K x = 0 (e.g. periodic
% boundary conditions), and the negativity condition of 'P', 'Q', 'FP' and
% 'FQ' carries + K'*Z + Z'*K, with Z a free operator, which vanishes on
% ker K: the condition is imposed on ker K only. 'Pd' and 'Qd' are then
% posed as 'P' and 'Q': without K they are not sufficient for stability.
%
% INPUT
% - PIE:      'pie_struct' object of a 1-D or 2-D PIE, with T and A 'opvar'
%             (1-D) or 'opvar2d' (2-D) objects, or of a container PIE in any
%             number of variables, with T and A 'copvar' objects
%             (EXEC_PIE_SOP), and K a 'copvar' object (the constraint
%             K x = 0 on the state) or empty; its other operators are not
%             used;
% - settings: (optional) 'struct' of LPI settings, as returned by
%             lpisettings_sop:
%               dd1, dd12, options1, options12, override1   1-D degrees of P;
%               lpi.pos        the storage row of P, in every dimension     % MMP, 11/17/2026
%                              (POSCOPVAR); settings_2d (optional, no       % MMP, 11/17/2026
%                              tier field) the 16-basis 2-D degrees,        % MMP, 11/17/2026
%                              from SETTINGS_2D_STOCK;                      % MMP, 11/17/2026
%               eppos   the margin of P on every space, in every dimension
%                       (EXEC_SOP_SETTINGS);
%               epneg   the decay rate, >= 0;
%               sos_opts  solver options;
%               sop     container options (EXEC_SOP_SETTINGS), among them
%                       dw_feas, max_raise_feas and eppos_stab.
%                       For 'Q': deg_Q (degrees of Q; [] = lpi.var.deg),  % MMP, 11/18/2026
%                       max_raise_Q ([] = max_raise_feas) and raise_Q,     % MMP, 11/18/2026
%                       'slack', 'deg' or 'both', what a raise adds.       % MMP, 11/18/2026
%             Default: in 1-D lpisettings_sop('heavy'), in 2-D and higher
%             lpisettings_sop('light') with sos_opts.simplify 1, each with
%             the margin of its tier, epneg 0, and epneg 1e-5 for the forms
%             'FP', 'FQ';
% - form:     (optional) char, 'P' (default), 'Pd', 'Q', 'Qd', 'FP' or 'FQ'.
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - P:        'copvar' object, the solved positive operator ([] when the
%             solver returned no solution);
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded', 'unknown' or 'residual' (cuADMM,
%                          which decides no feasibility: never certified);
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               form, epneg, eppos, nv   the LPI posed and the dimension;
%               storage    'copvar' object, the solved P;
%               D, N       'cdopvar' objects, the derivative operator
%                          (with K'*Z + Z'*K when K is used) and the slack;
%               maxre      double, in 1-D with settings.sop.witness the
%                          largest real part of the eigenvalues of the
%                          Legendre discretization of (A,T), else NaN;
%               witness    'struct' of that discretization ([] if none);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          as returned by EXEC_INFO_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_stability_sop
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
% Initial coding MMP, 10/28/2026: one executive for the six stability
%                executives of 10/08/2026 (PIETOOLS_PDEstability_sop, _dual,
%                PIETOOLS_PIE2PDEstability_sop, _dual,
%                PIETOOLS_stability_2D_sop, _dual_2D), which called
%                private/stability_build_sop; their LPIs and default settings
%                are the forms here, the program of each form unchanged.
% MMP, 10/28/2026: Forms 'FP' and 'FQ', the finite-energy stability LPIs of
%                PIETOOLS_FEPDEstability and PIETOOLS_FEPIE2PDEstability (and
%                their 2-D versions), which had no container counterpart. As
%                the stock: P >= 0 with no margin, the decay term epneg*T'*T
%                (2*epneg in 2-D), epneg > 0 required, epneg 1e-5 in the
%                default settings; 'FQ' posed as 'FP' in 2-D.
% MMP, 10/29/2026: settings.sop.max_raise_feas read per dimension, entry
%                min(nv,end), as dw_feas is (the 2-D raises cost far more than
%                the 1-D ones and certified nothing in the 10/28 example
%                library baseline); a scalar serves every dimension as before.
% MMP, 10/29/2026: PIEs in 3 or more variables and container PIEs: the PIE
%                goes through EXEC_PIE_SOP, and the 2-D rules hold for every
%                nv >= 2 (the Q forms posed as the P forms, 2*epneg, the 2-D
%                default settings); EXEC_TOOLS_SOP supplies the N-D storage.
%                The slack takes opts.cap = capfor(D) of EXEC_TOOLS_SOP (the
%                per-direction bounds of settings.sop.cap_delta; empty by
%                default, every program unchanged), and st.sop takes tl.sop,
%                each per-dimension vector of the settings tiers reduced to
%                its entry for the dimension of the PIE.
% MMP, 11/04/2026: The PIE operators are read as they are (EXEC_PIE_SOP
%                returns 'copvar' objects, which OP2COPVAR_SOP returned
%                unchanged), and P is returned as a 'copvar' object, not an
%                opvar / opvar2d (maintainer: the executives return
%                containers). Programs unchanged.
% MMP, 11/07/2026: The constraint K x = 0 of the PIE (PIE.K, set by
%                convert_PIETOOLS_PDE for rank-deficient boundary conditions,
%                e.g. periodic, with auxiliary states) is used: in 'P', 'Q',
%                'FP', 'FQ' the negativity condition, which acts on the state
%                x, carries + K'*Z + Z'*K, Z free, so it holds on ker K only.
%                Without it the LPI was posed off ker K too, where the
%                auxiliary dynamics miss a boundary condition: x_t = x_ss - x
%                with periodic or Neumann BCs was infeasible, now certified
%                up to the exact decay rate. 'Pd' and 'Qd' ignore K and
%                certified unstable PIEs: with K they are posed as 'P' and
%                'Q'. Programs without K unchanged. Help: status 'residual'
%                (LPI_CLASSIFY_SOP, a cuADMM exit since 11/04/2026).
% MMP, 11/08/2026: Default settings lpisettings_sop('light') in 2-D and
%                higher, lpisettings_sop('heavy') in 1-D (was the stock
%                lpisettings of the same levels, retired with its folder):
%                the same storage degrees, with the container options of
%                those tiers (maintainer, 11/08/2026).
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): log condensed, help
%                defaults named after lpisettings_sop, a stale body note
%                (conversion back to opvar) removed; code unchanged.
% MMP, 11/09/2026: One decay rate and one margin in every dimension
%                (maintainer: the executives apply eppos and epneg the same
%                way in every dimension). The decay term is 2*epneg in every
%                form and dimension, as in two or more variables and in
%                PIETOOLS_WELL_POSEDNESS_SOP: in 1-D it was epneg, so a
%                given epneg certified a decay rate of epneg/2 for ||T x||
%                in 1-D and epneg in 2-D. The default settings take the
%                margin of their tier in every dimension (were eppos 1e-4,
%                eppos2 1e-6 in 1-D and settings_2d.eppos 1e-2 in 2-D and
%                higher), and the finite-energy default epneg 1e-5 is set
%                once. The entries of 10/28/2026 and 10/29/2026 on the 2-D
%                factor 2 and the 2-D default settings no longer apply.
% MMP, 11/12/2026: The program takes the settings of the run (the last
%                input of LPIPROGRAM), so LPI_INEQ_SOP reads its degree
%                defaults from them instead of the light tier (maintainer,
%                11/12/2026: "Settings on the program"). Programs unchanged:
%                every tier sets lpi.ineq.dD = 1, and with dw and like
%                given wmin does not apply.
% MMP, 11/17/2026: Form 'Q' in the direct PIE-to-PDE form, in every
%                dimension (settings revision, work package 4; settings
%                report 11/16/2026, Sec. 8.3): Q free with the shape of T
%                (LPIVAR_CDOPVAR at the entry of lpi.var.deg for the number
%                of variables of the PIE on every block, as the direct test
%                of the report; with DEG omitted the 1-D blocks of a mixed
%                PIE take the 1-D entry, and the slack rule carries that
%                degree into every space), the
%                symmetry T'*Q = Q'*T as an equality, P = T'*Q >= eppos2*T'*T
%                and -(A'*Q + Q'*A + 2*epneg*P) >= 0 through LPI_INEQ_SOP,
%                each slack sized from its own support (no 'like': sized
%                like P = T'*Q the lift degree misses the support of
%                A'*Q + Q'*A; S17 at 0.5 r2* UNKNOWN at 1.31M variables
%                where the own-support slack of the report certifies at
%                435k) with the dw, margin and cap of the 'P' form; no
%                storage; the base raise of the direct form is 0 (the
%                own-sized slack reaches the support; dw_feas, dw_P, dw_Q
%                are allowances for slacks sized like the storage) and the
%                raise loop adds. The storage-plus-equality form
%                (storage P, T'*Q = P, Q degrees from P), kept commented out,
%                fails S3 in two variables and S17 under every storage; the
%                direct form certifies S3, S17, S18 (report, Sec. 8.3). In
%                two or more variables 'Q' is no longer posed as 'P'; eppos2
%                is the least L2 margin (entry 2 of tl.eppos in 1-D, as
%                before). Output, witness and dual as in the other forms;
%                info.storage is the solved T'*Q. 'Qd' and 'FQ' unchanged.
%                The help lines on the 'Q' form are updated; they read
%                before: "The 'Q' forms exist in 1-D only; for a PIE in two
%                or more variables, 'Q', 'Qd' and 'FQ' are posed as 'P',
%                'Pd' and 'FP'."
% MMP, 11/17/2026: Help: settings_2d and settings_nd are no tier fields
%                since work package 5; the storage is the row lpi.pos, and
%                the stock 2-D struct is SETTINGS_2D_STOCK. Code unchanged.
%                The lines read before: "settings_2d  2-D degrees of P;
%                settings_nd  degrees and margin of P in 3 or more
%                variables (EXEC_TOOLS_SOP);".
% MMP, 11/18/2026: The direct form 'Q' separated from the storage form 'P'
%                (maintainer: "make the separation"; "Q is a generalized
%                storage where P form is a special case when Q=PT"): Q at
%                settings.sop.deg_Q (a cell per dimension; [] = lpi.var.deg,
%                which stays the degree of the synthesis free operators),
%                its raises settings.sop.max_raise_Q ([] = max_raise_feas)
%                and, per raise, what settings.sop.raise_Q names: the slack
%                weight ('slack', as before), the first entry of deg_Q
%                ('deg') or both. The tiers carry the values of 11/17/2026,
%                so every program of that date is unchanged (ab_sep.m).
%                Measured reason: at light, Q at [2 2 2] cannot represent
%                P*T of the light storage (kernel degree 2) and a T of
%                degree 1, and no slack raise recovers the five DDE
%                thresholds that [4 2 2] restores (DECISIONS D92, D94).
%                Help: the 'Q' line no longer says "no storage".

% Default settings and form, as the stock executive of the PIE's dimension
if nargin<3 || isempty(form),   form = 'P';     end
if ~isa(PIE,'pie_struct')
    error('PIETOOLS_stability_sop:pie','The PIE should be a ''pie_struct'' object.');
end
% PIE = initialize(PIE);                                                    % MMP, 10/29/2026 (was)
PIE = exec_pie_sop(PIE);                    % opvar or container PIE        % MMP, 10/29/2026
fe = any(strcmp(form,{'FP','FQ'}));         % finite-energy forms           % MMP, 10/28/2026
if nargin<2 || isempty(settings)
%   if PIE.dim==2                                                           % MMP, 10/29/2026 (was)
    if PIE.dim>=2                           % the 2-D defaults in N-D too   % MMP, 10/29/2026
%       settings = lpisettings('light');                                    % MMP, 11/08/2026 (was)
        settings = lpisettings_sop('light');                                % MMP, 11/08/2026
        settings.sos_opts.simplify = 1;
%       settings.settings_2d.eppos = 1e-2*ones(4,1);    settings.settings_2d.epneg = 0; % MMP, 11/09/2026 (was)
%       if fe,  settings.settings_2d.epneg = 1e-5;  end % the stock FE default % MMP, 11/09/2026 (was)
    else
%       settings = lpisettings('heavy');                                    % MMP, 11/08/2026 (was)
        settings = lpisettings_sop('heavy');                                % MMP, 11/08/2026
%       settings.eppos = 1e-4;  settings.eppos2 = 1e-6;     settings.epneg = 0; % MMP, 11/09/2026 (was)
%       if fe,  settings.epneg = 1e-5;  end             % the stock FE default % MMP, 11/09/2026 (was)
    end
    % The tier's margin and decay rate in every dimension; finite energy    % MMP, 11/09/2026
    % needs a positive decay rate, 1e-5 as in the stock executives          % MMP, 11/09/2026
    if fe,  settings.epneg = 1e-5;  end                                     % MMP, 11/09/2026
end
label = ['stability_sop ' form];

% Complete the container options and read the tools of the PIE's dimension
% (storage of P, stock slack, margins)
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension (EXEC_TOOLS_SOP)        % MMP, 10/29/2026
% if ~any(strcmp(form,{'P','Pd','Q','Qd'}))                                 % MMP, 10/28/2026 (was)
%     error('PIETOOLS_stability_sop:form','FORM should be ''P'', ''Pd'', ''Q'' or ''Qd''.'); % MMP, 10/28/2026 (was)
if ~any(strcmp(form,{'P','Pd','Q','Qd','FP','FQ'}))                         % MMP, 10/28/2026
    error('PIETOOLS_stability_sop:form','FORM should be ''P'', ''Pd'', ''Q'', ''Qd'', ''FP'' or ''FQ''.'); % MMP, 10/28/2026
end
% if nv==2 && strcmp(form,'Q'),   form = 'P';     end     % no 2-D Q forms  % MMP, 10/29/2026 (was)
% if nv==2 && strcmp(form,'Qd'),  form = 'Pd';    end                       % MMP, 10/29/2026 (was)
% if nv==2 && strcmp(form,'FQ'),  form = 'FP';    end                       % MMP, 10/29/2026 (was)
% if nv>=2 && strcmp(form,'Q'),   form = 'P';     end   % no Q forms in N-D % MMP, 11/17/2026 (was); MMP, 10/29/2026
% 'Q' is posed in every dimension (direct form, header 11/17/2026)          % MMP, 11/17/2026
if nv>=2 && strcmp(form,'Qd'),  form = 'Pd';    end                         % MMP, 10/29/2026
if nv>=2 && strcmp(form,'FQ'),  form = 'FP';    end                         % MMP, 10/29/2026
% The constraint K x = 0 on the state, if any (empty for most PIEs)         % MMP, 11/07/2026
Km = PIE.K;     hasK = isa(Km,'copvar') && size(Km,1)>0;                    % MMP, 11/07/2026
% Dual forms test the PIE without K, which is not sufficient: 'Pd', 'Qd'    % MMP, 11/07/2026
% certified x_t = x_ss + x, periodic BCs (unstable). Posed as 'P', 'Q'.     % MMP, 11/07/2026
if hasK && strcmp(form,'Pd'),   form = 'P';     end                         % MMP, 11/07/2026
if hasK && strcmp(form,'Qd'),   form = 'Q';     end                         % MMP, 11/07/2026
% coercive = any(strcmp(form,{'P','Pd'}));                                  % MMP, 10/28/2026 (was)
coercive = any(strcmp(form,{'P','Pd','FP'}));   % 'FP' sizes its slack as 'P' % MMP, 10/28/2026

% Degree raise of the slack, margins of P and the decay term
dwbase = st.sop.dw_Q;   if coercive,    dwbase = st.sop.dw_P;   end
if ~isempty(st.sop.dw_feas),    dwbase = st.sop.dw_feas(min(nv,end));  end
% 'Q': the slacks are sized from their own support, which the rule reaches % MMP, 11/17/2026
% by construction; dw_feas, dw_P and dw_Q are allowances for slacks sized  % MMP, 11/17/2026
% like the storage, so the base is 0 and the raise loop adds (base 1 in    % MMP, 11/17/2026
% two variables gave 2-5x the variables: S3 1.19M against 383k, S18 1.18M % MMP, 11/17/2026
% against 228k, S17 refused at 5.2M against 435k certified)                % MMP, 11/17/2026
if strcmp(form,'Q'),    dwbase = 0;     end                                 % MMP, 11/17/2026
eppos = tl.eppos;   epneg = tl.epneg;
if ~isempty(st.sop.eppos_stab),     eppos(:) = st.sop.eppos_stab;   end
if ~isempty(st.sop.max_raise_feas)
%   st.sop.max_raise = st.sop.max_raise_feas;                               % MMP, 10/29/2026 (was)
    st.sop.max_raise = st.sop.max_raise_feas(min(nv,end));  % per dimension % MMP, 10/29/2026
end
% The direct form: its own raise count and the degrees of Q (header         % MMP, 11/18/2026
% 11/18/2026); the storage form keeps max_raise_feas, and lpi.var.deg       % MMP, 11/18/2026
% stays with the synthesis LPIs                                             % MMP, 11/18/2026
if strcmp(form,'Q') && ~isempty(st.sop.max_raise_Q)                         % MMP, 11/18/2026
    st.sop.max_raise = st.sop.max_raise_Q;      % reduced by EXEC_TOOLS_SOP % MMP, 11/18/2026
end                                                                         % MMP, 11/18/2026
qdeg = st.sop.deg_Q;    if isempty(qdeg),   qdeg = st.lpi.var.deg;  end     % MMP, 11/18/2026
if ~iscell(qdeg),   qdeg = {qdeg};  end                                     % MMP, 11/18/2026
qdeg = qdeg{min(nv,end)};                                                   % MMP, 11/18/2026
% epfac = 1;  if nv==2,   epfac = 2;  end           % the 2-D LPI carries 2 epneg % MMP, 10/29/2026 (was)
% epfac = 1;  if nv>=2,   epfac = 2;  end     % 2 epneg, as the 2-D LPI     % MMP, 11/09/2026 (was)
epfac = 2;          % 2 epneg in every form and dimension: decay rate epneg % MMP, 11/09/2026
% Finite energy needs a positive decay term, as in the stock executives     % MMP, 10/28/2026
if fe && epneg==0                                                           % MMP, 10/28/2026
    error('PIETOOLS_stability_sop:epneg','Finite energy stability check requires settings.epneg >0. Please change your settings.'); % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026

% The PIE operators as containers
% EXEC_PIE_SOP returns every operator as a 'copvar': read as they are       % MMP, 11/04/2026
% Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);                    % MMP, 11/04/2026 (was)
Tm = PIE.T;     Am = PIE.A;                                                 % MMP, 11/04/2026

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;

% The solved storage, as a 'copvar'                                         % MMP, 11/04/2026
Psol = [];  P = [];
try
    Psol = getsol_lpivar_sop(prog,aux.Pm);
%   P = tl.toop(Psol);                                                      % MMP, 11/04/2026 (was)
    P = Psol;                                                               % MMP, 11/04/2026
catch
    P = Psol;
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info_sop(prog,aux,st,label);
info.form = form;   info.epneg = epneg;     info.eppos = eppos;     info.nv = nv;
info.storage = Psol;    info.D = aux.Dm;    info.N = aux.Nm;

% Numerical witness (1-D): the spectrum of the discretized pencil (A,T)
info.witness = [];  info.maxre = NaN;
if st.sop.witness && nv==1
    Wt = pie_witness_sop(PIE,'rate',st.sop);
    info.witness = Wt;
    if Wt.ok,   info.maxre = Wt.maxre;  end
end
if st.sop.verbose
    fprintf('  [%s] %s: certificate at epneg %g: %d; numerical spectrum max Re %.6g\n', ...
            label,cl.status,epneg,info.certified,info.maxre);
end

    function [prog,aux] = assemble(dw)
        % [PROG,AUX] = ASSEMBLE(DW) builds the LPI of FORM with the slack
        % degree raised by DWBASE + DW.
        % INPUT  dw:    nonnegative integer, the raise of this pass.
        % OUTPUT prog:  'struct', the LPI program;
        %        aux:   'struct' with fields obj ([] here: a feasibility
        %               LPI), Pm ('cdopvar' P), Dm (the derivative
        %               operator), Nm (the slack), li (slack info of
        %               LPI_INEQ_SOP), form.
        % New program; record the equality rows when the dual is kept
%       prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);             % MMP, 11/12/2026 (was)
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);          % MMP, 11/12/2026
        if st.sop.keepdual,     prog.sopeq = {};    end
        % Raise of this pass: in the direct form what raise_Q names, the    % MMP, 11/18/2026
        % slack weight ('slack', 'both') and the first entry of deg_Q       % MMP, 11/18/2026
        % ('deg', 'both'); the other forms raise the slack weight           % MMP, 11/18/2026
        dws = dw;   dwq = 0;                                                % MMP, 11/18/2026
        if strcmp(form,'Q')                                                 % MMP, 11/18/2026
            if strcmp(st.sop.raise_Q,'deg'),    dws = 0;    end             % MMP, 11/18/2026
            if ~strcmp(st.sop.raise_Q,'slack'), dwq = dw;   end             % MMP, 11/18/2026
        end                                                                 % MMP, 11/18/2026
        % The positive operator P and the derivative operator Dm of FORM
        switch form
            case 'P'
                [prog,Pm] = tl.storage(prog,Tm,'out');
                if ~all(eppos==0),  Pm = Pm + tl.eye(1,eppos);  end
                Dm = Tm'*Pm*Am + Am'*Pm*Tm;
                if epneg~=0,    Dm = Dm + epfac*epneg*(Tm'*Pm*Tm);  end
            case 'Pd'
                [prog,Pm] = tl.storage(prog,Tm,'in');
                if ~all(eppos==0),  Pm = Pm + tl.eye(2,eppos);  end
                Dm = Tm*Pm*Am' + Am*Pm*Tm';
                if epneg~=0,    Dm = Dm + epfac*epneg*(Tm*Pm*Tm');  end
            case 'Q'
                % BEGIN MMP, 11/17/2026: the direct form (header 11/17/2026;
                % settings report 11/16/2026, Sec. 8.3). Replaced: the      % MMP, 11/17/2026
                % storage P on the input spaces of T, the equality T'*Q = P % MMP, 11/17/2026
                % and Q at the degrees of P; it failed S3 in 2-D and S17.   % MMP, 11/17/2026
%               [prog,Pm] = tl.storage(prog,Tm,'in');                       % MMP, 11/17/2026 (was)
%               Pm = Pm + eppos(2)*(Tm'*Tm);                                % MMP, 11/17/2026 (was)
%               % Q free, from the output space of T to its input space     % MMP, 11/17/2026 (was)
%               [spO,dmO] = copvar_space_list(Tm,'out');    [spI,dmI] = copvar_space_list(Tm,'in'); % MMP, 11/17/2026 (was)
%               [prog,Qm] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs_sop(Pm)); % MMP, 11/17/2026 (was)
%               prog = lpi_eq_sop(prog,Tm'*Qm - Pm);                        % MMP, 11/17/2026 (was)
%               Dm = Am'*Qm + Qm'*Am;                                       % MMP, 11/17/2026 (was)
%               if epneg~=0,    Dm = Dm + epneg*Pm;     end                 % MMP, 11/09/2026 (was)
%               if epneg~=0,    Dm = Dm + epfac*epneg*Pm;   end             % MMP, 11/17/2026 (was); MMP, 11/09/2026
                % Q free with the shape of T, at the entry of lpi.var.deg  % MMP, 11/17/2026
                % for the number of variables of the PIE on every block,    % MMP, 11/17/2026
                % the convention of the direct test of the report (Sec.     % MMP, 11/17/2026
                % 9.2). DEG omitted gives the 1-D blocks of a mixed PIE the % MMP, 11/17/2026
                % 1-D entry, which the registry-wide slack rule then carries% MMP, 11/17/2026
                % into every space: S3 654k variables against 383k, S17     % MMP, 11/17/2026
                % UNKNOWN at 1.31M against 435k certified                   % MMP, 11/17/2026
%               [prog,Qm] = lpivar_cdopvar(prog,Tm);                        % MMP, 11/17/2026 (was)
%               [prog,Qm] = lpivar_cdopvar(prog,Tm,st.lpi.var.deg{min(nv,end)}); % MMP, 11/18/2026 (was); MMP, 11/17/2026
                q = qdeg;   if isnumeric(q),   q(1) = q(1) + dwq;   end     % MMP, 11/18/2026
                [prog,Qm] = lpivar_cdopvar(prog,Tm,q);   % deg_Q + ladder   % MMP, 11/18/2026
                % T'*Q self-adjoint; P its symmetric part, = T'*Q on the    % MMP, 11/17/2026
                % equality, written symmetric for LPI_INEQ_SOP              % MMP, 11/17/2026
                prog = lpi_eq_sop(prog,Tm'*Qm - Qm'*Tm);                    % MMP, 11/17/2026
                Pm = 0.5*(Tm'*Qm + Qm'*Tm);                                 % MMP, 11/17/2026
                % P >= eppos2*T'*T, eppos2 the least L2 margin (eppos(2)    % MMP, 11/17/2026
                % in 1-D, as before); slack options as for Dm below         % MMP, 11/17/2026
%               o = struct('like',Pm,'dw',dwbase+dw,'margin',st.sop.margin,'cap',tl.capfor(Pm)); % MMP, 11/17/2026 (was)
                % No storage: each slack is sized from its own support,     % MMP, 11/17/2026
                % as the direct test of the report (Sec. 9.2); sized like   % MMP, 11/17/2026
                % P = T'*Q the lift degree misses the support of A'*Q+Q'*A  % MMP, 11/17/2026
                % (S17: UNKNOWN at 1.31M variables where 435k certify)      % MMP, 11/17/2026
%               o = struct('dw',dwbase+dw,'margin',st.sop.margin,'cap',tl.capfor(Pm)); % MMP, 11/18/2026 (was); MMP, 11/17/2026
                o = struct('dw',dwbase+dws,'margin',st.sop.margin,'cap',tl.capfor(Pm)); % MMP, 11/18/2026
                prog = lpi_ineq_sop(prog,Pm - min(eppos(min(2,end):end))*(Tm'*Tm),o); % MMP, 11/17/2026
                Dm = Am'*Qm + Qm'*Am;                                       % MMP, 11/17/2026
                if epneg~=0,    Dm = Dm + epfac*epneg*Pm;   end             % MMP, 11/17/2026
                % END MMP, 11/17/2026
            case 'Qd'
                [prog,Pm] = tl.storage(prog,Tm,'out');
                Pm = Pm + eppos(2)*(Tm*Tm');
                [spO,dmO] = copvar_space_list(Tm,'in');     [spI,dmI] = copvar_space_list(Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs_sop(Pm));
                prog = lpi_eq_sop(prog,Tm*Qm - Pm);
                Dm = Am*Qm + Qm'*Am';
%               if epneg~=0,    Dm = Dm + epneg*Pm;     end                 % MMP, 11/09/2026 (was)
                if epneg~=0,    Dm = Dm + epfac*epneg*Pm;   end             % MMP, 11/09/2026
            case 'FP'                                                       % MMP, 10/28/2026
                % P >= 0 with no margin; the decay term on T'T              % MMP, 10/28/2026
                [prog,Pm] = tl.storage(prog,Tm,'out');                      % MMP, 10/28/2026
                Dm = Tm'*Pm*Am + Am'*Pm*Tm + epfac*epneg*(Tm'*Tm);          % MMP, 10/28/2026
            case 'FQ'                                                       % MMP, 10/28/2026
                % P >= 0 with no margin, Q free with T'*Q = P, as 'Q'; the  % MMP, 10/28/2026
                % decay term on T'T                                         % MMP, 10/28/2026
                [prog,Pm] = tl.storage(prog,Tm,'in');                       % MMP, 10/28/2026
                [spO,dmO] = copvar_space_list(Tm,'out');    [spI,dmI] = copvar_space_list(Tm,'in'); % MMP, 10/28/2026
                [prog,Qm] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs_sop(Pm)); % MMP, 10/28/2026
                prog = lpi_eq_sop(prog,Tm'*Qm - Pm);                        % MMP, 10/28/2026
%               Dm = Am'*Qm + Qm'*Am + epneg*(Tm'*Tm);                      % MMP, 11/09/2026 (was)
                Dm = Am'*Qm + Qm'*Am + epfac*epneg*(Tm'*Tm);                % MMP, 11/09/2026
        end
        % BEGIN MMP, 11/07/2026: Dm <= 0 on ker K only
        % In the primal forms Dm acts on the state x; with Z free from the
        % state space to the output space of K, x'(K'*Z + Z'*K)x = 0 when
        % K x = 0, so Dm + K'*Z + Z'*K <= 0 gives Dm <= 0 on ker K.
        if hasK && any(strcmp(form,{'P','Q','FP','FQ'}))
            [spK,dmK] = copvar_space_list(Km,'out');    [spI,dmI] = copvar_space_list(Tm,'in');
            [prog,Zm] = lpivar_cdopvar(prog,struct('out',dmK,'in',dmI),struct('out',{spK},'in',{spI}),tl.dom,get_lpivar_degs_sop(Pm));
            Dm = Dm + Km'*Zm + Zm'*Km;
        end
        % END MMP, 11/07/2026
        % Dm <= 0: the stock slack and an equality, or the slack sized from P
        if strcmp(st.sop.slack,'stock')
            [prog,Nm] = tl.slack(prog,Dm);
            [prog,tag] = lpi_eq_sop(prog,Dm+Nm,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
%           o = struct('like',Pm,'dw',dwbase+dw,'margin',st.sop.margin);    % MMP, 10/29/2026 (was)
%           o = struct('like',Pm,'dw',dwbase+dw,'margin',st.sop.margin,'cap',tl.capfor(Dm)); % MMP, 11/18/2026 (was); MMP, 10/29/2026
            o = struct('like',Pm,'dw',dwbase+dws,'margin',st.sop.margin,'cap',tl.capfor(Dm)); % MMP, 10/29/2026
            % 'Q': no storage, sized from the support of Dm itself (the    % MMP, 11/17/2026
            % 'Q' branch above)                                             % MMP, 11/17/2026
            if strcmp(form,'Q'),    o = rmfield(o,'like');  end             % MMP, 11/17/2026
            [prog,Nm,li] = lpi_ineq_sop(prog,-Dm,o);
        end
        aux = struct('obj',[],'Pm',Pm,'Dm',Dm,'Nm',Nm,'li',li,'form',form);
    end
end
