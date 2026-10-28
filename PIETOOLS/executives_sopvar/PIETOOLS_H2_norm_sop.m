function [prog,W,gam,R,Q,info] = PIETOOLS_H2_norm_sop(PIE,settings,form)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,W,GAM,R,Q,INFO] = PIETOOLS_H2_NORM_SOP(PIE,SETTINGS,FORM) computes
% an upper bound on the H2 norm from w to z of the PIE
%
%       T x_t = A x + B1 w,      z = C1 x,
%
% by an LPI on 'cdopvar' operators, solved as an SDP that minimizes gam.
% FORM selects the LPI; R and W are positive operators, Q a free operator:
%
%   'c'    R >= 0, T*Q = R, W >= 0 on the space of z, gam >= trace(W),
%          [-gam I, B1'; B1, Q'*A' + A*Q] <= 0,  [W, C1*Q; Q'*C1', R] >= 0;
%   'o'    R >= 0, T'*Q = R, W >= 0 on the space of w, gam >= trace(W),
%          [-gam I, C1; C1', Q'*A + A'*Q] <= 0,  [W, B1'*Q; Q'*B1, R] >= 0;
%   'cco'  W >= 0 (2-D: W >= eppos I on R^n and L2[s1,s2]),
%          (A*W)*T' + T*(W*A') + B1*B1' <= 0,  gam >= trace(C1*W*C1');
%   'oco'  W >= 0 (2-D: as 'cco'),
%          (A'*W)*T + T'*(W*A) + C1'*C1 <= 0,  gam >= trace(B1'*W*B1);
%
% the trace taken over the finite-dimensional block. The bound GAM is gam
% in the forms 'c', 'o' and sqrt(gam) in 'cco', 'oco'. While the solver
% returns UNKNOWN or INFEASIBLE, the program is rebuilt with the degree of
% the slacks raised by one, up to settings.sop.max_raise times.
%
% INPUT
% - PIE:      'pie_struct' object of a 1-D or 2-D PIE, with T, A, B1 and C1
%             'opvar' (1-D) or 'opvar2d' (2-D) objects; Tw zero (no
%             disturbance at the boundary);
% - settings: (optional) 'struct' of LPI settings, as returned by lpisettings
%             or lpisettings_sop:
%               dd1, dd12, options1, options12, override1   1-D degrees of R
%                       (of W in 'cco', 'oco');
%               settings_2d   2-D degrees;
%               settings_2d.eppos (2-D)  margins of W in 'cco', 'oco';
%               sos_opts  solver options;
%               sop     container options (EXEC_SOP_SETTINGS): dw_Q for 'c'
%                       and 'o', dw_P for 'cco' and 'oco', max_raise.
%             Default lpisettings('light');
% - form:     (optional) char, 'c' (default), 'o', 'cco' or 'oco'.
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - W:        'opvar' (1-D) or 'opvar2d' (2-D) object, the solved W ([] when
%             the solver returned no solution);
% - gam:      double, the certified bound on the H2 norm; NaN when the
%             solver returned no value;
% - R, Q:     'opvar' / 'opvar2d' objects, the solved R and Q in the forms
%             'c' and 'o'; [] in 'cco' and 'oco';
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded' or 'unknown';
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               form, gam, gam_lpi, nv  the LPI posed, the bound, the LPI
%                          objective, the dimension;
%               W, R, Q    as the outputs;
%               h2_num     double, in 1-D with settings.sop.witness the H2
%                          norm of the Legendre discretization of the PIE,
%                          else NaN; gap  gam/h2_num - 1;
%               witness    'struct' of the discretization ([] if none);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          as returned by EXEC_INFO_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_H2_norm_sop
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
% Initial coding MMP, 10/28/2026: one executive for the eight H2 norm
%                executives of 10/08/2026 (PIETOOLS_H2_norm_c_sop, _o_sop,
%                _c_coercive_sop, _o_coercive_sop and their 2-D versions
%                _2D_c_non_coercive_sop, _2D_o_non_coercive_sop,
%                _2D_c_sop, _2D_o_sop), which called private/h2_build_sop
%                (MMP, 10/08/2026, 10/09/2026); their LPIs are the forms
%                here, the program of each unchanged.

% Default settings and form; check the PIE
if nargin<2 || isempty(settings),   settings = lpisettings('light');    end
if nargin<3 || isempty(form),       form = 'c';     end
if ~isa(PIE,'pie_struct')
    error('PIETOOLS_H2_norm_sop:pie','The PIE should be a ''pie_struct'' object.');
end
PIE = initialize(PIE);
label = ['H2_norm_sop ' form];

% Complete the container options and read the tools of the PIE's dimension
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
if ~(PIE.Tw==0),    error('PIETOOLS_H2_norm_sop:Tw','The H2 norm LPI cannot be solved with disturbances at the boundary.');  end
if ~any(strcmp(form,{'c','o','cco','oco'}))
    error('PIETOOLS_H2_norm_sop:form','FORM should be ''c'', ''o'', ''cco'' or ''oco''.');
end
qform = any(strcmp(form,{'c','o'}));

% Degree raise of the slacks and the margins
dwbase = st.sop.dw_P;   if qform,    dwbase = st.sop.dw_Q;   end
eppos = tl.eppos;

% The PIE operators as containers, and the identities on w and z
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);
Bw = op2copvar_sop(PIE.B1); Cz = op2copvar_sop(PIE.C1);
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE; the bound is sqrt(gam) in the coercive forms
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
gam_lpi = cl.obj;
gam = gam_lpi;
if any(strcmp(form,{'cco','oco'})),     gam = sqrt(gam_lpi);    end

% The solved W (and R, Q), converted back to opvar / opvar2d
W = [];     R = [];     Q = [];
try
    W = tl.toop(getsol_lpivar_sop(prog,aux.Wm));
    if qform
        R = tl.toop(getsol_lpivar_sop(prog,aux.Rm));
        Q = tl.toop(getsol_lpivar_sop(prog,aux.Qm));
    end
catch ME
    if st.sop.verbose,  fprintf('  [%s] output extraction: %s\n',label,ME.message);  end
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info_sop(prog,aux,st,label);
info.form = form;   info.gam = gam;  info.gam_lpi = gam_lpi;    info.nv = nv;
info.W = W;     info.R = R;     info.Q = Q;

% Numerical witness (1-D): the H2 norm of the discretized PIE
info.witness = [];  info.h2_num = NaN;  info.gap = NaN;
if st.sop.witness && nv==1
    Wt = pie_witness_sop(PIE,'h2',st.sop);
    info.witness = Wt;
    if Wt.ok
        info.h2_num = Wt.h2;
        if ~isnan(gam) && Wt.h2>0,  info.gap = gam/Wt.h2 - 1;   end
    end
end
if st.sop.verbose
    fprintf('  [%s] %s: value %.8g (LPI objective %.8g); numerical H2 norm %.8g; gap %.2e\n', ...
            label,cl.status,gam,gam_lpi,info.h2_num,info.gap);
end

    function [prog,aux] = assemble(dw)
        % [PROG,AUX] = ASSEMBLE(DW) builds the H2 LPI of FORM with the slack
        % degrees raised by DWBASE + DW.
        % INPUT  dw:    nonnegative integer, the raise of this pass.
        % OUTPUT prog:  'struct', the LPI program;
        %        aux:   'struct' with fields obj (the 'dpvar' gam), Rm, Qm,
        %               Wm ('cdopvar' R, Q, W; Rm, Qm [] in 'cco', 'oco'), li
        %               (slack info of LPI_INEQ_SOP for the negativity
        %               constraint), form.
        % New program; record the equality rows when the dual is kept
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        if st.sop.keepdual,     prog.sopeq = {};    end
        % gam as the objective, gam >= 0 except in 'o' (as the stock 1-D files)
        [prog,gamv] = lpidecvar(prog,'gam');
        if ~strcmp(form,'o'),   prog = lpi_ineq(prog,gamv);     end
        prog = lpisetobj(prog,gamv);
        % The operators, the negativity (Dneg) and positivity (Dpos)
        % constraints of FORM, the support of the slacks and the trace
        Rm = [];    Qm = [];    Wm = [];
        switch form
            case 'c'
                [sp,dm] = copvar_space_list(Tm,'out');
                [prog,Rm] = tl.storage(prog,Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,dm,sp,tl.dom,get_lpivar_degs_sop(Rm));
                prog = lpi_eq_sop(prog,Tm*Qm - Rm);
                [prog,Wm] = wdecl(prog,Cz,'out');
                Dneg = [-(gamv*Iw),  Bw';
                        Bw,          Qm'*Am' + Am*Qm];
                Dpos = [Wm,         Cz*Qm;
                        Qm'*Cz',    Rm];
                like = Rm;  tr = trace_rn_sop(Wm);
            case 'o'
                [sp,dm] = copvar_space_list(Tm,'out');
                [prog,Rm] = tl.storage(prog,Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,dm,sp,tl.dom,get_lpivar_degs_sop(Rm));
                prog = lpi_eq_sop(prog,Tm'*Qm - Rm);
                [prog,Wm] = wdecl(prog,Bw,'in');
                Dneg = [-(gamv*Iz),  Cz;
                        Cz',         Qm'*Am + Am'*Qm];
                Dpos = [Wm,         Bw'*Qm;
                        Qm'*Bw,     Rm];
                like = Rm;  tr = trace_rn_sop(Wm);
            case 'cco'
                [prog,Wm] = tl.storage(prog,Tm,'out');
                Wm = add_eppos2d(Wm);
                Dneg = (Am*Wm)*Tm' + Tm*(Wm*Am') + Bw*Bw';
                Dpos = [];  like = Wm;  tr = trace_rn_sop((Cz*Wm)*Cz');
            case 'oco'
                [prog,Wm] = tl.storage(prog,Tm,'out');
                Wm = add_eppos2d(Wm);
                Dneg = (Am'*Wm)*Tm + Tm'*(Wm*Am) + Cz'*Cz;
                Dpos = [];  like = Wm;  tr = trace_rn_sop((Bw'*Wm)*Bw);
        end
        % Dneg <= 0 and Dpos >= 0: the stock slacks and equalities, or the
        % slacks sized from the storage
        if strcmp(st.sop.slack,'stock')
            [prog,Nneg] = tl.slack(prog,Dneg);
            [prog,tag] = lpi_eq_sop(prog,Nneg + Dneg,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
            if ~isempty(Dpos)
                [prog,Npos] = tl.slack(prog,Dpos);
                prog = lpi_eq_sop(prog,Npos - Dpos,'symmetric');
            end
        else
            o = struct('like',like,'dw',dwbase+dw,'margin',st.sop.margin);
            [prog,~,li] = lpi_ineq_sop(prog,-Dneg,o);
            if ~isempty(Dpos)
                [prog,~,~] = lpi_ineq_sop(prog,Dpos,o);
            end
        end
        % gam >= trace
        prog = lpi_ineq(prog,gamv - tr);
        aux = struct('obj',gamv,'Rm',Rm,'Qm',Qm,'Pm',[],'Zm',[],'Wm',Wm,'li',li,'form',form);
    end

    function [prog,Wm] = wdecl(prog,X,side)
        % [PROG,WM] = WDECL(PROG,X,SIDE) declares W >= 0 on the SIDE space of
        % X in the forms 'c' and 'o'.
        % INPUT  prog:  'struct', the LPI program;
        %        X:     'copvar' object (C1 for 'c', B1 for 'o');
        %        side:  char, 'out' or 'in'.
        % OUTPUT prog:  'struct', the program with the Gram matrices of W;
        %        Wm:    'cdopvar' object, W.
        % 1-D: poslpivar's default degrees; 2-D: the storage degrees on a
        % space with variables, the identity basis on R^n
        if nv==1
            [prog,Wm] = poslpivar_sop(prog,X,[],struct(),side,PIE.dom);
        elseif any(any(X.(['space_' side])))
            [prog,Wm] = tl.storage(prog,X,side);
        else
            [prog,Wm] = poslpivar_sop(prog,X,[],struct(),side);
        end
    end

    function Wm = add_eppos2d(Wm)
        % WM = ADD_EPPOS2D(WM) adds eppos I on R^n and L2[s1,s2] to W in 2-D
        % (blkdiag(eppos(1) I, 0, 0, eppos(4) I)); returns W unchanged in 1-D.
        % INPUT/OUTPUT  Wm: 'cdopvar' object, W.
        if nv==2
            e = eppos;  e(2:3) = 0;
            if ~all(e==0),  Wm = Wm + tl.eye(1,e);  end
        end
    end
end
