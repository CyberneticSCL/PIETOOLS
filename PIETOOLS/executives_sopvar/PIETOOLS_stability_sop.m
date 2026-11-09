function [prog,P,info] = PIETOOLS_stability_sop(PIE,settings,form)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,INFO] = PIETOOLS_STABILITY_SOP(PIE,SETTINGS,FORM) tests stability
% of the PIE
%
%       T x_t = A x
%
% by searching for a Lyapunov certificate, posed as an LPI on 'cdopvar'
% operators and solved as an SDP. FORM selects the LPI; with c = 1 in 1-D
% and c = 2 in two or more spatial variables,
%
%   'P'   P >= eppos I,          T'*P*A + A'*P*T + c*epneg*T'*P*T <= 0;
%   'Pd'  P >= eppos I,          T*P*A' + A*P*T' + c*epneg*T*P*T' <= 0;
%   'Q'   P >= eppos2*T'*T,  T'*Q = P,   A'*Q + Q'*A + epneg*P <= 0;
%   'Qd'  P >= eppos2*T*T',  T*Q = P,    A*Q + Q'*A' + epneg*P <= 0;
%   'FP'  P >= 0,                T'*P*A + A'*P*T + c*epneg*T'*T <= 0;
%   'FQ'  P >= 0,            T'*Q = P,   A'*Q + Q'*A + epneg*T'*T <= 0;
%
% with P a positive operator, Q a free operator, and eppos I the identity
% scaled by the margin of each space. The forms 'FP', 'FQ' test finite
% energy, int_0^inf ||T x(t)||^2 dt <= C ||T x(0)||^2 for some C > 0, and
% need epneg > 0. The 'Q' forms exist in 1-D only; for a PIE in two or
% more variables, 'Q', 'Qd' and 'FQ' are posed as 'P', 'Pd' and 'FP'. While
% the solver
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
%               settings_2d                                 2-D degrees of P;
%               settings_nd    degrees and margin of P in 3 or more
%                              variables (EXEC_TOOLS_SOP);
%               eppos, eppos2 (1-D), settings_2d.eppos (2-D)  margins of P;
%               epneg   the decay rate, >= 0;
%               sos_opts  solver options;
%               sop     container options (EXEC_SOP_SETTINGS), among them
%                       dw_feas, max_raise_feas and eppos_stab.
%             Default: in 1-D lpisettings_sop('heavy') with eppos 1e-4,
%             eppos2 1e-6, epneg 0; in 2-D and higher lpisettings_sop('light')
%             with sos_opts.simplify 1, settings_2d.eppos 1e-2 on every
%             space, epneg 0; epneg 1e-5 for the forms 'FP', 'FQ';
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
        settings.settings_2d.eppos = 1e-2*ones(4,1);    settings.settings_2d.epneg = 0;
        if fe,  settings.settings_2d.epneg = 1e-5;  end % the stock FE default % MMP, 10/28/2026
    else
%       settings = lpisettings('heavy');                                    % MMP, 11/08/2026 (was)
        settings = lpisettings_sop('heavy');                                % MMP, 11/08/2026
        settings.eppos = 1e-4;  settings.eppos2 = 1e-6;     settings.epneg = 0;
        if fe,  settings.epneg = 1e-5;  end             % the stock FE default % MMP, 10/28/2026
    end
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
if nv>=2 && strcmp(form,'Q'),   form = 'P';     end     % no Q forms in N-D % MMP, 10/29/2026
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
eppos = tl.eppos;   epneg = tl.epneg;
if ~isempty(st.sop.eppos_stab),     eppos(:) = st.sop.eppos_stab;   end
if ~isempty(st.sop.max_raise_feas)
%   st.sop.max_raise = st.sop.max_raise_feas;                               % MMP, 10/29/2026 (was)
    st.sop.max_raise = st.sop.max_raise_feas(min(nv,end));  % per dimension % MMP, 10/29/2026
end
% epfac = 1;  if nv==2,   epfac = 2;  end           % the 2-D LPI carries 2 epneg % MMP, 10/29/2026 (was)
epfac = 1;  if nv>=2,   epfac = 2;  end     % 2 epneg, as the 2-D LPI       % MMP, 10/29/2026
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
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        if st.sop.keepdual,     prog.sopeq = {};    end
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
                [prog,Pm] = tl.storage(prog,Tm,'in');
                Pm = Pm + eppos(2)*(Tm'*Tm);
                % Q free, from the output space of T to its input space
                [spO,dmO] = copvar_space_list(Tm,'out');    [spI,dmI] = copvar_space_list(Tm,'in');
                [prog,Qm] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs_sop(Pm));
                prog = lpi_eq_sop(prog,Tm'*Qm - Pm);
                Dm = Am'*Qm + Qm'*Am;
                if epneg~=0,    Dm = Dm + epneg*Pm;     end
            case 'Qd'
                [prog,Pm] = tl.storage(prog,Tm,'out');
                Pm = Pm + eppos(2)*(Tm*Tm');
                [spO,dmO] = copvar_space_list(Tm,'in');     [spI,dmI] = copvar_space_list(Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs_sop(Pm));
                prog = lpi_eq_sop(prog,Tm*Qm - Pm);
                Dm = Am*Qm + Qm'*Am';
                if epneg~=0,    Dm = Dm + epneg*Pm;     end
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
                Dm = Am'*Qm + Qm'*Am + epneg*(Tm'*Tm);                      % MMP, 10/28/2026
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
            o = struct('like',Pm,'dw',dwbase+dw,'margin',st.sop.margin,'cap',tl.capfor(Dm)); % MMP, 10/29/2026
            [prog,Nm,li] = lpi_ineq_sop(prog,-Dm,o);
        end
        aux = struct('obj',[],'Pm',Pm,'Dm',Dm,'Nm',Nm,'li',li,'form',form);
    end
end
