function [prog,P,R,omega,info] = PIETOOLS_well_posedness_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,R,OMEGA,INFO] = PIETOOLS_WELL_POSEDNESS_SOP(PIE,SETTINGS) tests
% well-posedness of the 1-D PIE
%
%       T x_t = A x
%
% by the LPI, posed on 'cdopvar' operators and solved as an SDP,
%
%   P >= eppos I,  R >= I,
%   T'*P*A + A'*P*T - 2*omega*T'*P*T <= 0,
%   (T - A)*R*(T - A)' - eppos I >= 0,
%
% with P and R positive operators and omega = -settings.epneg. While the
% solver returns UNKNOWN or INFEASIBLE, the program is rebuilt with the
% degree of both slacks raised by one, up to settings.sop.max_raise_feas
% (its 1-D entry; or max_raise) times.
%
% INPUT
% - PIE:      'pie_struct' object of a 1-D PIE, with T and A 'opvar' objects,
%             or 'copvar' objects (a 1-D container PIE, EXEC_PIE_SOP); its
%             other operators are not used; a PIE in two or more variables
%             is an error;
% - settings: (optional) 'struct' of LPI settings, as returned by lpisettings
%             or lpisettings_sop:
%               dd1, dd12, options1, options12, override1   degrees of P, R;
%               eppos, eppos2   margins of P;
%               epneg   minus omega;
%               sos_opts  solver options;
%               sop     container options (EXEC_SOP_SETTINGS), among them
%                       dw_feas, max_raise_feas and eppos_stab.
%             Default lpisettings('heavy') with sos_opts.simplify 1, eppos =
%             eppos2 = 1e-2, epneg 0.
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - P, R:     'opvar' objects, the solved P and R ([] when the solver
%             returned no solution);
% - omega:    double, the omega of the LPI;
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded' or 'unknown';
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               omega, eppos   as posed;
%               D1, D2     'cdopvar' objects, the two constrained operators;
%               maxre      double, with settings.sop.witness the largest
%                          real part of the eigenvalues of the Legendre
%                          discretization of (A,T), else NaN;
%               witness    'struct' of that discretization ([] if none);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          as returned by EXEC_INFO_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_well_posedness_sop
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
% MMP, 10/28/2026: settings.sop.dw_feas (its 1-D entry), max_raise_feas and
%                eppos_stab replace dw_P, max_raise and eppos when set: the
%                1-D feasibility LPIs certify at dw 0 with 2.0-3.7x fewer
%                variables (measured 10/16/2026). Unset, nothing changes.
% MMP, 10/28/2026: Header to the header protocol, step comments in the body;
%                code unchanged. The deleted header text compared the file
%                with the stock executive and named the helpers it calls.
% MMP, 10/29/2026: settings.sop.max_raise_feas may hold one entry per
%                dimension (PIETOOLS_stability_sop reads entry min(nv,end));
%                this 1-D executive reads its 1-D entry, as for dw_feas.
% MMP, 10/29/2026: The PIE goes through EXEC_PIE_SOP (INITIALIZE for an
%                opvar PIE, as before; a 1-D container PIE accepted), and
%                every PIE in two or more variables is rejected, not only a
%                2-D one.
% MMP, 10/29/2026: Each slack takes opts.cap = capfor(target) of
%                EXEC_TOOLS_SOP, the per-direction bounds of
%                settings.sop.cap_delta; empty by default, so every program
%                is unchanged unless cap_delta is set.
% MMP, 10/29/2026: st.sop takes tl.sop of EXEC_TOOLS_SOP: each per-dimension
%                vector (dw_feas, dw_P, dw_Q, max_raise, max_raise_feas,
%                cap_delta) reduced to its entry for the dimension of the
%                PIE, so the settings tiers may hold [1-D 2-D 3-D] values;
%                scalar settings are unchanged.

% Default settings, as the stock executive; check the PIE
if nargin<2 || isempty(settings)
    settings = lpisettings('heavy');
    settings.sos_opts.simplify = 1;
    settings.eppos = 1e-2;  settings.eppos2 = 1e-2;   settings.epneg = 0;
end
label = 'well_posedness_sop';
st = exec_sop_settings(settings);
if ~isa(PIE,'pie_struct'),  error('PIETOOLS_well_posedness_sop:pie','The PIE should be a ''pie_struct''.');  end
% PIE = initialize(PIE);                                                    % MMP, 10/29/2026 (was)
% if PIE.dim==2,  error('PIETOOLS_well_posedness_sop:dim','Well-posedness analysis of 2D PDEs is currently not supported.');   end % MMP, 10/29/2026 (was)
PIE = exec_pie_sop(PIE);                    % opvar or container PIE        % MMP, 10/29/2026
if PIE.dim>=2,  error('PIETOOLS_well_posedness_sop:dim','Well-posedness analysis of PDEs in two or more variables is currently not supported.');   end % MMP, 10/29/2026

% Tools of the PIE's dimension (storage, stock slack, margins); omega, and
% the degree raise and raise count of a feasibility LPI
tl = exec_tools_sop(PIE,st);
st.sop = tl.sop;    % the options of this dimension (EXEC_TOOLS_SOP)        % MMP, 10/29/2026
eppos = tl.eppos;   omega = -tl.epneg;
if ~isempty(st.sop.eppos_stab),     eppos(:) = st.sop.eppos_stab;   end     % MMP, 10/28/2026
dwbase = st.sop.dw_P;                                                       % MMP, 10/28/2026
if ~isempty(st.sop.dw_feas),    dwbase = st.sop.dw_feas(1);     end         % MMP, 10/28/2026
if ~isempty(st.sop.max_raise_feas)                                          % MMP, 10/28/2026
%   st.sop.max_raise = st.sop.max_raise_feas;                               % MMP, 10/29/2026 (was)
    st.sop.max_raise = st.sop.max_raise_feas(1);    % the 1-D entry         % MMP, 10/29/2026
end                                                                         % MMP, 10/28/2026

% The PIE operators as containers
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);

% Build and solve, raising the slack degrees while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;

% The solved P and R, converted back to opvars
P = [];     R = [];
try
    P = tl.toop(getsol_lpivar_sop(prog,aux.Pm));
    R = tl.toop(getsol_lpivar_sop(prog,aux.Rm));
catch ME
    if st.sop.verbose,  fprintf('  [%s] solution extraction: %s\n',label,ME.message);  end
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info_sop(prog,aux,st,label);
info.omega = omega;     info.eppos = eppos;
info.D1 = aux.Dm1;      info.D2 = aux.Dm2;

% Numerical witness: the spectrum of the discretized pencil (A,T)
info.witness = [];      info.maxre = NaN;
if st.sop.witness
    Wt = pie_witness_sop(PIE,'rate',st.sop);
    info.witness = Wt;
    if Wt.ok,   info.maxre = Wt.maxre;  end
end
if st.sop.verbose
    fprintf('  [%s] %s: certificate at omega %g: %d; numerical spectrum max Re %.6g\n', ...
            label,cl.status,omega,info.certified,info.maxre);
end

    function [prog,aux] = assemble(dw)
        % [PROG,AUX] = ASSEMBLE(DW) builds the well-posedness LPI with both
        % slack degrees raised by DWBASE + DW.
        % INPUT  dw:    nonnegative integer, the raise of this pass.
        % OUTPUT prog:  'struct', the LPI program;
        %        aux:   'struct' with fields obj ([]: a feasibility LPI), Pm,
        %               Rm ('cdopvar' P and R), Dm1, Dm2 (the two
        %               constrained operators), li (slack info of
        %               LPI_INEQ_SOP for Dm1, with degrees2 those for Dm2).
        % New program; record the equality rows when the dual is kept
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        if st.sop.keepdual,     prog.sopeq = {};    end
        % P >= eppos I on the output space of T, R >= I on its input space
        [prog,Pm] = tl.storage(prog,Tm,'out');
        [prog,Rm] = tl.storage(prog,Tm,'in');
        I1 = tl.eye(1,eppos);
        Pm = Pm + I1;
        Rm = Rm + tl.eye(2,[1;1]);
        % Dissipativity Dm1 <= 0 and surjectivity Dm2 >= 0
        Dm1 = Tm'*Pm*Am + Am'*Pm*Tm;
        if omega~=0,    Dm1 = Dm1 - 2*omega*(Tm'*Pm*Tm);    end
        TA = Tm - Am;
        Dm2 = TA*Rm*TA' - I1;
        % The stock slacks and equalities, or the slacks sized from P and R
        if strcmp(st.sop.slack,'stock')
            [prog,N1] = tl.slack(prog,Dm1);
            [prog,tag] = lpi_eq_sop(prog,Dm1+N1,'symmetric');
            [prog,N2] = tl.slack(prog,Dm2);
            prog = lpi_eq_sop(prog,Dm2-N2,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
%           o1 = struct('like',Pm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin); % MMP, 10/28/2026 (was)
%           o1 = struct('like',Pm,'dw',dwbase+dw,'margin',st.sop.margin);   % MMP, 10/29/2026 (was)
            o1 = struct('like',Pm,'dw',dwbase+dw,'margin',st.sop.margin,'cap',tl.capfor(Dm1)); % MMP, 10/29/2026
            [prog,~,li] = lpi_ineq_sop(prog,-Dm1,o1);
%           o2 = struct('like',Rm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin); % MMP, 10/28/2026 (was)
%           o2 = struct('like',Rm,'dw',dwbase+dw,'margin',st.sop.margin);   % MMP, 10/29/2026 (was)
            o2 = struct('like',Rm,'dw',dwbase+dw,'margin',st.sop.margin,'cap',tl.capfor(Dm2)); % MMP, 10/29/2026
            [prog,~,li2] = lpi_ineq_sop(prog,Dm2,o2);
            li.degrees2 = li2.degrees;
        end
        aux = struct('obj',[],'Pm',Pm,'Rm',Rm,'Dm1',Dm1,'Dm2',Dm2,'li',li);
    end
end
