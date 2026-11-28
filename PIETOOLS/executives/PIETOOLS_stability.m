function [prog,P,info] = PIETOOLS_stability(PIE,settings,form)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,INFO] = PIETOOLS_STABILITY(PIE,SETTINGS,FORM) tests stability
% of the PIE
%
%       T x_t = A x
%
% by searching for a Lyapunov certificate, posed as an LPI on 'dopvar'
% operators and solved as an SDP. FORM selects the LPI:
%
%   'P'   P >= eppos I,          T'*P*A + A'*P*T + 2*epneg*T'*P*T <= 0;
%   'Pd'  P >= eppos I,          T*P*A' + A*P*T' + 2*epneg*T*P*T' <= 0;
%   'Q'   P >= eppos2*T'*T,  T'*Q = P,   A'*Q + Q'*A + 2*epneg*P <= 0;
%         (posed directly: T'*Q = Q'*T, P = T'*Q the storage, Q at
%         settings.sop.deg_Q; 'P' is the case Q = P*T)
%   'Qd'  P >= eppos2*T*T',  T*Q = P,    A*Q + Q'*A' + 2*epneg*P <= 0;
%   'FP'  P >= 0,                T'*P*A + A'*P*T + 2*epneg*T'*T <= 0;
%   'FQ'  P >= 0,            T'*Q = P,   A'*Q + Q'*A + 2*epneg*T'*T <= 0;
%
% with P a positive operator, Q a free operator, eppos I the identity
% scaled by the margin of each space (eppos2 the margin on L2), and epneg
% the decay rate, in every dimension: a solution of 'P' gives
% ||T x(t)|| <= M exp(-epneg t) ||T x(0)|| for some M > 0. The forms 'FP', 'FQ' test finite
% energy, int_0^inf ||T x(t)||^2 dt <= C ||T x(0)||^2 for some C > 0, and
% need epneg > 0. 'Q' is posed in every dimension, eppos2 the least L2
% margin; 'Qd' and 'FQ' with the declared storage P in 1-D, and in two or   % MMP, 11/21/2026
% more variables as 'Q' is, with Q free of the shape of T' ('Qd': the LPI   % MMP, 11/21/2026
% of 'Q' on the pair (T',A')) or of T ('FQ'). While the solver              % MMP, 11/21/2026
% returns UNKNOWN or INFEASIBLE, the program is rebuilt with the degree of
% the negativity slack raised by one, up to settings.sop.max_raise_feas (its
% entry for the dimension of the PIE, the last entry for a higher one; or
% max_raise) times.
%
% If PIE.K is nonempty, the PIE state satisfies K x = 0 (e.g. periodic
% boundary conditions), and the negativity condition of 'P', 'Q', 'FP' and
% 'FQ' carries + K'*Z + Z'*K, with Z a free operator, which vanishes on
% ker K: the condition is imposed on ker K only. 'Pd' and 'Qd' are then
% errors: they test the PIE without K, which is not sufficient for stability. % MMP, 11/21/2026
% No form is ever replaced by another.                                      % MMP, 11/21/2026
%
% INPUT
% - PIE:      'pie_struct' object of a PIE in any number of spatial
%             variables, initialized when its T is not yet an 'opvar'
%             object; T and An 'opvar' objects, and K an 'opvar' object
%             (the constraint K x = 0 on the state) or empty; its other
%             operators are not used;
% - settings: (optional) 'struct' of LPI settings:
%               lpi.pos        the storage row of P, in every dimension;
%                              settings_2d (optional) the 16-basis 2-D
%                              degrees, in place of lpi.pos in two
%                              variables;
%               eppos   the margin of P on every space, in every dimension;
%               epneg   the decay rate, >= 0;
%               sos_opts  solver options;
%               sop     the options of the executives, among them
%                       dw_feas, max_raise_feas and eppos_stab.
%                       For 'Q': deg_Q (degrees of Q; [] = lpi.var.deg),
%                       max_raise_Q ([] = max_raise_feas) and raise_Q,
%                       'slack', 'deg' or 'both', what a raise adds.
%             Default: in 1-D the 'heavy' tier, in 2-D and higher the
%             'light' tier with sos_opts.simplify 1, each with the margin
%             of its tier, epneg 0, and epneg 1e-5 for the forms
%             'FP', 'FQ';
% - form:     (optional) char, 'P' (default), 'Pd', 'Q', 'Qd', 'FP' or 'FQ'.
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - P:        'opvar' object, the solved positive operator ([] when the
%             solver returned no solution);
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded', 'unknown' or 'residual' (cuADMM,
%                          which decides no feasibility: never certified);
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               form, epneg, eppos, nv   the LPI posed and the dimension;
%               storage    'opvar' object, the solved P;
%               D, N       'dopvar' objects, the derivative operator
%                          (with K'*Z + Z'*K when K is used) and the slack;
%               maxre      double, in 1-D with settings.sop.witness the
%                          largest real part of the eigenvalues of the
%                          Legendre discretization of (A,T), else NaN;
%               witness    'struct' of that discretization ([] if none);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          the last raise and the loop history, the slack
%                          mode, its degrees and terms, the size of the SDP,
%                          and the dual kernel of the negativity constraint
%                          with its 1-norm;
%               aux        'struct', the operators of the last pass: Pm (P), % MMP, 11/21/2026
%                          Dm, Nm, li (the slack degrees and terms), form;  % MMP, 11/21/2026
%                          Npm, lip: the slack of the margin of P and its   % MMP, 11/21/2026
%                          info in the direct forms, [] in the others.      % MMP, 11/21/2026
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_stability
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
% Initial coding MMP, 10/28/2026
% MMP, 11/21/2026: 'Qd' and 'FQ' are posed in two or more variables, in the
%                direct form of 'Q' (Q free; for 'Qd' the LPI of 'Q' on the
%                pair (T',A'); for 'FQ' P >= 0 and the decay term
%                2*epneg*T'*T). They were replaced by 'Pd' and 'FP' without
%                a message, so a request for either ran another LPI
%                (maintainer, 11/21: a bug). 1-D unchanged. With K, 'Pd' and
%                'Qd' are now errors in place of their replacement by 'P'
%                and 'Q': a form certifies its own notion of stability, and
%                the executive never poses a form other than FORM
%                (maintainer, 11/21: "there should NEVER be a swap").
%                Help text before: "'Qd' and
%                'FQ' exist in 1-D only, and for a PIE in two or more
%                variables are posed as 'Pd' and 'FP'." and "'Pd' and 'Qd'
%                are then posed as 'P' and 'Q': without K they are not
%                sufficient for stability."
% MMP, 11/21/2026: info.aux returns the operators of the last pass, which an
%                a posteriori check of the solve needs (the storage fixed,
%                the slack sought; maintainer, 11/21: in the toolbox).
% MMP, 11/21/2026: info.aux also returns the slack of the margin of P in the
%                direct forms (Npm, lip): sized from the numeric P instead,
%                the check built a different slack and failed a solve whose
%                -D passed (PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input).

% Default settings and form, by the dimension of the PIE
if nargin<3 || isempty(form),   form = 'P';     end
if ~isa(PIE,'pie_struct')
    error('PIETOOLS_stability:pie','The PIE should be a ''pie_struct'' object.');
end
PIE = exec_pie(PIE);                        % every operator an 'opvar' or 0
fe = any(strcmp(form,{'FP','FQ'}));         % finite-energy forms
if nargin<2 || isempty(settings)
    if PIE.dim>=2                           % the 2-D defaults in N-D too
        settings = lpisettings('light');
        settings.sos_opts.simplify = 1;
    else
        settings = lpisettings('heavy');
    end
    % The margin and decay rate of the tier in every dimension; finite
    % energy needs a positive decay rate, 1e-5
    if fe,  settings.epneg = 1e-5;  end
end
label = ['stability ' form];

% Complete the options of the executives and read the tools of the
% dimension of the PIE (storage of P, 'stock' slack, margins)
st = exec_sop_settings(settings);
tl = exec_tools(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension
if ~any(strcmp(form,{'P','Pd','Q','Qd','FP','FQ'}))
    error('PIETOOLS_stability:form','FORM should be ''P'', ''Pd'', ''Q'', ''Qd'', ''FP'' or ''FQ''.');
end
% % 'Q' is posed in every dimension; 'Qd' and 'FQ' in 1-D only              % MMP, 11/21/2026 (was)
% if nv>=2 && strcmp(form,'Qd'),  form = 'Pd';    end                       % MMP, 11/21/2026 (was)
% if nv>=2 && strcmp(form,'FQ'),  form = 'FP';    end                       % MMP, 11/21/2026 (was)
% The constraint K x = 0 on the state, if any (empty for most PIEs)
Km = PIE.K;     hasK = isa(Km,'opvar') && size(Km,1)>0;
% The dual forms test the PIE without K, which is not sufficient for
% % stability: with K, 'Pd' and 'Qd' are posed as 'P' and 'Q'               % MMP, 11/21/2026 (was)
% stability: with K, 'Pd' and 'Qd' are errors (another form would certify   % MMP, 11/21/2026
% another notion of stability)                                              % MMP, 11/21/2026
% if hasK && strcmp(form,'Pd'),   form = 'P';     end                       % MMP, 11/21/2026 (was)
% if hasK && strcmp(form,'Qd'),   form = 'Q';     end                       % MMP, 11/21/2026 (was)
if hasK && any(strcmp(form,{'Pd','Qd'}))                                    % MMP, 11/21/2026
    error('PIETOOLS_stability:formK',['PIE.K is nonempty (the state satisfies K x = 0), ', ...
        'which the dual form ''%s'' does not test; use ''P'' or ''Q''.'],form); % MMP, 11/21/2026
end                                                                         % MMP, 11/21/2026
% 'Qd' and 'FQ' in two or more variables: the direct form of 'Q', Q free    % MMP, 11/21/2026
% with the storage T'*Q (T*Q for 'Qd', the LPI of 'Q' on (T',A')); in 1-D   % MMP, 11/21/2026
% the declared storage P with T'*Q = P (T*Q = P)                            % MMP, 11/21/2026
dq = nv>=2 && any(strcmp(form,{'Qd','FQ'}));                                % MMP, 11/21/2026
coercive = any(strcmp(form,{'P','Pd','FP'}));   % 'FP' sizes its slack as 'P'

% Degree raise of the slack, margins of P and the decay term
dwbase = st.sop.dw_Q;   if coercive,    dwbase = st.sop.dw_P;   end
if ~isempty(st.sop.dw_feas),    dwbase = st.sop.dw_feas(min(nv,end));  end
% 'Q': the slacks are sized from their own support, which the rule reaches
% by construction; dw_feas, dw_P and dw_Q are allowances for slacks sized
% like the storage, so the base is 0 and the raise loop adds (base 1 in
% two variables gave 2 to 5 times the variables)
% if strcmp(form,'Q'),    dwbase = 0;     end                               % MMP, 11/21/2026 (was)
if strcmp(form,'Q') || dq,  dwbase = 0;     end     % every direct form     % MMP, 11/21/2026
eppos = tl.eppos;   epneg = tl.epneg;
if ~isempty(st.sop.eppos_stab),     eppos(:) = st.sop.eppos_stab;   end
if ~isempty(st.sop.max_raise_feas)
    st.sop.max_raise = st.sop.max_raise_feas(min(nv,end));  % per dimension
end
% The direct form: its own raise count and the degrees of Q; the storage
% form keeps max_raise_feas, and lpi.var.deg stays with the synthesis LPIs
% if strcmp(form,'Q') && ~isempty(st.sop.max_raise_Q)                       % MMP, 11/21/2026 (was)
if (strcmp(form,'Q') || dq) && ~isempty(st.sop.max_raise_Q)                 % MMP, 11/21/2026
    st.sop.max_raise = st.sop.max_raise_Q;      % already the entry of this dimension
end
qdeg = st.sop.deg_Q;    if isempty(qdeg),   qdeg = st.lpi.var.deg;  end
if ~iscell(qdeg),   qdeg = {qdeg};  end
qdeg = qdeg{min(nv,end)};
epfac = 2;          % 2 epneg in every form and dimension: decay rate epneg
% Finite energy needs a positive decay term
if fe && epneg==0
    error('PIETOOLS_stability:epneg','Finite energy stability check requires settings.epneg >0. Please change your settings.');
end

% The PIE operators, 'opvar' objects
Tm = PIE.T;     Am = PIE.A;
% 'Qd' in two or more variables: 'Q' on the dual pair, since (T')'*Q = T*Q  % MMP, 11/21/2026
% and (A')'*Q + Q'*A' = A*Q + Q'*A'                                         % MMP, 11/21/2026
if dq && strcmp(form,'Qd'),  Tm = Tm';  Am = Am';  end                      % MMP, 11/21/2026

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop(@assemble,st,label);
cl = aux.classify;

% The solved storage, as an 'opvar'
Psol = [];  P = [];
try
    Psol = getsol_lpivar(prog,aux.Pm);
    P = Psol;
catch
    P = Psol;
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info(prog,aux,st,label);
info.form = form;   info.epneg = epneg;     info.eppos = eppos;     info.nv = nv;
info.storage = Psol;    info.D = aux.Dm;    info.N = aux.Nm;
info.aux = aux;     % the internals of the last pass, for a posteriori checks % MMP, 11/21/2026

% Numerical witness (1-D): the spectrum of the discretized pencil (A,T)
info.witness = [];  info.maxre = NaN;
if st.sop.witness && nv==1
    Wt = pie_witness(PIE,'rate',st.sop);
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
        %               LPI), Pm ('dopvar' P), Dm (the derivative
        %               operator), Nm (the slack), li (the info of the
        %               slack), form.
        % New program; record the equality rows when the dual is kept
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);
        if st.sop.keepdual,     prog.sopeq = {};    end
        % Raise of this pass: in the direct form what raise_Q names, the
        % slack weight ('slack', 'both') and the first entry of deg_Q
        % ('deg', 'both'); the other forms raise the slack weight
        dws = dw;   dwq = 0;
%       if strcmp(form,'Q')                                                 % MMP, 11/21/2026 (was)
        if strcmp(form,'Q') || dq                                           % MMP, 11/21/2026
            if strcmp(st.sop.raise_Q,'deg'),    dws = 0;    end
            if ~strcmp(st.sop.raise_Q,'slack'), dwq = dw;   end
        end
        % The positive operator P and the derivative operator Dm of FORM
        % (the direct 'Qd' and 'FQ' run the code of 'Q')                    % MMP, 11/21/2026
        fq = form;  if dq,  fq = 'Q';  end                                  % MMP, 11/21/2026
        Npm = [];   lip = [];   % slack of the margin of P (direct forms)   % MMP, 11/21/2026
%       switch form                                                         % MMP, 11/21/2026 (was)
        switch fq                                                           % MMP, 11/21/2026
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
                % The direct form: Q free with the shape of T, its degrees
                % deg_Q (the entry for the number of variables of the PIE)
                % on every block, the first entry raised by the ladder
                % Also 'Qd' (T, A their adjoints, set above) and 'FQ' in    % MMP, 11/21/2026
                % two or more variables                                     % MMP, 11/21/2026
                q = qdeg;   if isnumeric(q),   q(1) = q(1) + dwq;   end
                [prog,Qm] = lpivar(prog,Tm,q);           % deg_Q + ladder
                % T'*Q self-adjoint; P its symmetric part, = T'*Q on the
                % equality, written symmetric for the inequality
                prog = lpi_eq(prog,Tm'*Qm - Qm'*Tm);
                Pm = 0.5*(Tm'*Qm + Qm'*Tm);
                % P >= eppos2*T'*T, eppos2 the least L2 margin (eppos(2)
                % in 1-D). Each slack is sized from its own support: sized
                % like P = T'*Q, the lift degree misses the support of
                % A'*Q + Q'*A
                o = struct('dw',dwbase+dws,'margin',st.sop.margin,'cap',tl.capfor(Pm));
%               prog = lpi_ineq_sop(prog,Pm - min(eppos(min(2,end):end))*(Tm'*Tm),o); % MMP, 11/21/2026 (was)
                if fe       % 'FQ': P >= 0, no margin                       % MMP, 11/21/2026
%                   prog = lpi_ineq_sop(prog,Pm,o);                         % MMP, 11/21/2026 (was)
                    [prog,Npm,lip] = lpi_ineq(prog,Pm,o);                   % MMP, 11/21/2026
                else                                                        % MMP, 11/21/2026
%                   prog = lpi_ineq_sop(prog,Pm - min(eppos(min(2,end):end))*(Tm'*Tm),o); % MMP, 11/21/2026 (was)
                    [prog,Npm,lip] = lpi_ineq(prog,Pm - min(eppos(min(2,end):end))*(Tm'*Tm),o); % MMP, 11/21/2026
                end                                                         % MMP, 11/21/2026
                Dm = Am'*Qm + Qm'*Am;
%               if epneg~=0,    Dm = Dm + epfac*epneg*Pm;   end             % MMP, 11/21/2026 (was)
                if fe,  Dm = Dm + epfac*epneg*(Tm'*Tm);     % 'FQ': decay on T'*T % MMP, 11/21/2026
                elseif epneg~=0,    Dm = Dm + epfac*epneg*Pm;   end         % MMP, 11/21/2026
            case 'Qd'
                [prog,Pm] = tl.storage(prog,Tm,'out');
                Pm = Pm + eppos(2)*(Tm*Tm');
                [spO,dmO] = opvar_space_list(Tm,'in');     [spI,dmI] = opvar_space_list(Tm,'out');
                [prog,Qm] = lpivar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs(Pm));
                prog = lpi_eq(prog,Tm*Qm - Pm);
                Dm = Am*Qm + Qm'*Am';
                if epneg~=0,    Dm = Dm + epfac*epneg*Pm;   end
            case 'FP'
                % P >= 0 with no margin; the decay term on T'T
                [prog,Pm] = tl.storage(prog,Tm,'out');
                Dm = Tm'*Pm*Am + Am'*Pm*Tm + epfac*epneg*(Tm'*Tm);
            case 'FQ'
                % P >= 0 with no margin, Q free with T'*Q = P, as 'Q'; the
                % decay term on T'T
                [prog,Pm] = tl.storage(prog,Tm,'in');
                [spO,dmO] = opvar_space_list(Tm,'out');    [spI,dmI] = opvar_space_list(Tm,'in');
                [prog,Qm] = lpivar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs(Pm));
                prog = lpi_eq(prog,Tm'*Qm - Pm);
                Dm = Am'*Qm + Qm'*Am + epfac*epneg*(Tm'*Tm);
        end
        % Dm <= 0 on ker K only
        % In the primal forms Dm acts on the state x; with Z free from the
        % state space to the output space of K, x'(K'*Z + Z'*K)x = 0 when
        % K x = 0, so Dm + K'*Z + Z'*K <= 0 gives Dm <= 0 on ker K.
        if hasK && any(strcmp(form,{'P','Q','FP','FQ'}))
            [spK,dmK] = opvar_space_list(Km,'out');    [spI,dmI] = opvar_space_list(Tm,'in');
            [prog,Zm] = lpivar(prog,struct('out',dmK,'in',dmI),struct('out',{spK},'in',{spI}),tl.dom,get_lpivar_degs(Pm));
            Dm = Dm + Km'*Zm + Zm'*Km;
        end
        % Dm <= 0: the 'stock' slack and an equality, or the slack sized from P
        if strcmp(st.sop.slack,'stock')
            [prog,Nm] = tl.slack(prog,Dm);
            [prog,tag] = lpi_eq(prog,Dm+Nm,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
            o = struct('like',Pm,'dw',dwbase+dws,'margin',st.sop.margin,'cap',tl.capfor(Dm));
            % 'Q': sized from the support of Dm itself, not like the
            % generalized storage T'*Q
%           if strcmp(form,'Q'),    o = rmfield(o,'like');  end             % MMP, 11/21/2026 (was)
            if strcmp(form,'Q') || dq,  o = rmfield(o,'like');  end         % MMP, 11/21/2026
            [prog,Nm,li] = lpi_ineq(prog,-Dm,o);
        end
%       aux = struct('obj',[],'Pm',Pm,'Dm',Dm,'Nm',Nm,'li',li,'form',form); % MMP, 11/21/2026 (was)
        aux = struct('obj',[],'Pm',Pm,'Dm',Dm,'Nm',Nm,'li',li,'form',form, ...
                     'Npm',Npm,'lip',lip);  % the slack of the margin of P  % MMP, 11/21/2026
    end
end
