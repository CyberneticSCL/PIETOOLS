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
% - PIE:      'pie_struct' object of a PIE in one spatial variable,
%             initialized when its T is not yet a 'copvar' object; T and A
%             'copvar' objects; its other operators are not used; K empty
%             (a constraint K x = 0 on the state, e.g. periodic BCs, is an
%             error); a PIE in two or more variables is an error;
% - settings: (optional) 'struct' of LPI settings:
%               lpi.pos   the storage row of P and R;
%               eppos     the margin of P and of the surjectivity condition
%                         (eppos2, when set, the margin on L2);
%               epneg     minus omega;
%               sos_opts  solver options;
%               sop       the options of the executives, among them dw_feas
%                         (the degree raise of the slacks; else dw_P),
%                         max_raise_feas and eppos_stab (the margin of P,
%                         in place of eppos and eppos2).
%             Default: the 'heavy' tier with sos_opts.simplify 1,
%             eppos = eppos2 = 1e-2, epneg 0.
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - P, R:     'copvar' objects, the solved P and R ([] when the solver
%             returned no solution);
% - omega:    double, the omega of the LPI;
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded', 'unknown' or 'residual' (cuADMM,
%                          which decides no feasibility: never certified);
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               omega, eppos   as posed;
%               D1, D2     'cdopvar' objects, the two constrained operators;
%               maxre      double, with settings.sop.witness the largest
%                          real part of the eigenvalues of the Legendre
%                          discretization of (A,T), else NaN;
%               witness    'struct' of that discretization ([] if none);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          the last raise and the loop history, the slack
%                          mode, the degrees and terms of the slack of the
%                          dissipativity constraint, the size of the SDP,
%                          and the dual kernel of that constraint with its
%                          1-norm.
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

% Default settings; complete the options of the executives; check the PIE
if nargin<2 || isempty(settings)
    settings = lpisettings_sop('heavy');
    settings.sos_opts.simplify = 1;
    settings.eppos = 1e-2;  settings.eppos2 = 1e-2;   settings.epneg = 0;
end
label = 'well_posedness_sop';
st = exec_sop_settings(settings);
if ~isa(PIE,'pie_struct'),  error('PIETOOLS_well_posedness_sop:pie','The PIE should be a ''pie_struct''.');  end
PIE = exec_pie_sop(PIE);                    % every operator a 'copvar' or 0
if PIE.dim>=2,  error('PIETOOLS_well_posedness_sop:dim','Well-posedness analysis of PDEs in two or more variables is currently not supported.');   end
% With a constraint K x = 0 (auxiliary states) the range of T - A is a
% proper subspace: the surjectivity condition below cannot hold
if isa(PIE.K,'copvar') && size(PIE.K,1)>0
    error('PIETOOLS_well_posedness_sop:K','PIEs with an algebraic constraint K are not supported.');
end

% Tools of the dimension of the PIE (storage, 'stock' slack, margins);
% omega, and the degree raise and raise count of a feasibility LPI
tl = exec_tools_sop(PIE,st);
st.sop = tl.sop;    % the options of this dimension
eppos = tl.eppos;   omega = -tl.epneg;
if ~isempty(st.sop.eppos_stab),     eppos(:) = st.sop.eppos_stab;   end
dwbase = st.sop.dw_P;
if ~isempty(st.sop.dw_feas),    dwbase = st.sop.dw_feas(1);     end
if ~isempty(st.sop.max_raise_feas)
    st.sop.max_raise = st.sop.max_raise_feas(1);    % the 1-D entry
end

% The PIE operators, 'copvar' objects
Tm = PIE.T;     Am = PIE.A;

% Build and solve, raising the slack degrees while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;

% The solved P and R, as 'copvar' objects
P = [];     R = [];
try
    P = getsol_lpivar_sop(prog,aux.Pm);
    R = getsol_lpivar_sop(prog,aux.Rm);
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
        %               constrained operators), li (the info of the slack of
        %               Dm1, with degrees2 the degrees of the slack of Dm2).
        % New program; record the equality rows when the dual is kept
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);
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
        % The 'stock' slacks and equalities, or the slacks sized from P and R
        if strcmp(st.sop.slack,'stock')
            [prog,N1] = tl.slack(prog,Dm1);
            [prog,tag] = lpi_eq_sop(prog,Dm1+N1,'symmetric');
            [prog,N2] = tl.slack(prog,Dm2);
            prog = lpi_eq_sop(prog,Dm2-N2,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
            % The slack of Dm1 sized from P, that of Dm2 from R
            o1 = struct('like',Pm,'dw',dwbase+dw,'margin',st.sop.margin,'cap',tl.capfor(Dm1));
            [prog,~,li] = lpi_ineq_sop(prog,-Dm1,o1);
            o2 = struct('like',Rm,'dw',dwbase+dw,'margin',st.sop.margin,'cap',tl.capfor(Dm2));
            [prog,~,li2] = lpi_ineq_sop(prog,Dm2,o2);
            li.degrees2 = li2.degrees;
        end
        aux = struct('obj',[],'Pm',Pm,'Rm',Rm,'Dm1',Dm1,'Dm2',Dm2,'li',li);
    end
end
