function [prog,W,gam,R,Q,info] = PIETOOLS_H2_norm(PIE,settings,form)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,W,GAM,R,Q,INFO] = PIETOOLS_H2_NORM(PIE,SETTINGS,FORM) computes
% an upper bound on the H2 norm from w to z of the PIE
%
%       T x_t = A x + B1 w,      z = C1 x,
%
% by an LPI on 'dopvar' operators, solved as an SDP that minimizes gam.
% FORM selects the LPI; R and W are positive operators, Q a free operator:
%
%   'c'    R >= 0, T*Q = R, W >= 0 on the space of z, gam >= trace(W),
%          [-gam I, B1'; B1, Q'*A' + A*Q] <= 0,  [W, C1*Q; Q'*C1', R] >= 0;
%   'o'    R >= 0, T'*Q = R, W >= 0 on the space of w, gam >= trace(W),
%          [-gam I, C1; C1', Q'*A + A'*Q] <= 0,  [W, B1'*Q; Q'*B1, R] >= 0;
%   'cco'  W >= eppos I on R^n and on the space in every variable,
%          (A*W)*T' + T*(W*A') + B1*B1' <= 0,  gam >= trace(C1*W*C1');
%   'oco'  W with the margin of 'cco',
%          (A'*W)*T + T'*(W*A) + C1'*C1 <= 0,  gam >= trace(B1'*W*B1);
%
% the trace taken over the finite-dimensional block. The bound GAM is gam
% in the forms 'c', 'o' and sqrt(gam) in 'cco', 'oco'. While the solver
% returns UNKNOWN or INFEASIBLE, the program is rebuilt with the degree of
% the slacks raised by one, up to settings.sop.max_raise times.
%
% PIEs with an algebraic constraint (K, Kw or Ku with rows, e.g. from
% periodic boundary conditions) are not supported.
%
% INPUT
% - PIE:      'pie_struct' object of a PIE in any number of spatial
%             variables, initialized when its T is not yet an 'opvar'
%             object; T, A, B1 and C1 'opvar' objects; Tw zero or absent
%             (no disturbance at the boundary); K, Kw, Ku empty or without
%             rows; its other operators are not used;
% - settings: (optional) 'struct' of LPI settings:
%               lpi.pos   the storage row of R (of W in 'cco', 'oco'), in
%                         every dimension; settings_2d (optional) the
%                         16-basis 2-D degrees, in place of lpi.pos in two
%                         variables;
%               eppos     the margin of W in 'cco' and 'oco', in every
%                         dimension;
%               sos_opts  solver options;
%               sop       the options of the executives, among them dw_Q
%                         (the degree raise of the slacks in 'c' and 'o'),
%                         dw_P (in 'cco' and 'oco') and max_raise.
%             Default: the 'light' tier;
% - form:     (optional) char, 'c' (default), 'o', 'cco' or 'oco'.
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - W:        'opvar' object, the solved W ([] when
%             the solver returned no solution);
% - gam:      double, the certified bound on the H2 norm; NaN when the
%             solver returned no value;
% - R, Q:     'opvar' objects, the solved R and Q in the forms
%             'c' and 'o'; [] in 'cco' and 'oco';
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded', 'unknown' or 'residual' (cuADMM,
%                          which decides no feasibility: never certified);
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
%                          the last raise and the loop history, the slack
%                          mode, its degrees and terms, the size of the SDP,
%                          and the dual kernel of the negativity constraint
%                          with its 1-norm;
%               aux        'struct', the operators of the last pass: obj    % MMP, 11/21/2026
%                          (gam), Rm, Qm, Pm, Zm, Wm, li, form.             % MMP, 11/21/2026
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_H2_norm
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
% MMP, 11/21/2026: with K, 'c' and 'cco' are errors, where they were
%                replaced by 'o' and 'oco' without a message: each form
%                certifies its own bound (maintainer, 11/21: "there should
%                NEVER be a swap"). The programs of the admitted forms are
%                unchanged. Help text before: "which vanishes on ker K; 'c'
%                and 'cco' are posed as 'o' and 'oco'."
% DJ, 11/21/2026: Error for a PIE with an algebraic constraint (K, Kw, Ku
%                with rows);
% MMP, 11/21/2026: info.aux returns the operators of the last pass, which an
%                a posteriori check of the solve needs (maintainer, 11/21: in
%                the toolbox).

% Default settings and form; check the PIE
if nargin<2 || isempty(settings),  settings = lpisettings('light'); end
if nargin<3 || isempty(form),       form = 'c';     end
if ~isa(PIE,'pie_struct')
    error('PIETOOLS_H2_norm:pie','The PIE should be a ''pie_struct'' object.');
end
PIE = exec_pie(PIE);                        % every operator an 'opvar' or 0
label = ['H2_norm ' form];

% Complete the options of the executives and read the tools of the
% dimension of the PIE (storage, margins); boundary disturbances and an
% unknown FORM are errors
st = exec_sop_settings(settings);
tl = exec_tools(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension
if ~(PIE.Tw==0),    error('PIETOOLS_H2_norm:Tw','The H2 norm LPI cannot be solved with disturbances at the boundary.');  end
if ~any(strcmp(form,{'c','o','cco','oco'}))
    error('PIETOOLS_H2_norm:form','FORM should be ''c'', ''o'', ''cco'' or ''oco''.');
end
% The constraint K x = 0 on the state, if any (empty for most PIEs). 'c'
% % and 'cco' have no term for it: posed as 'o' and 'oco'                   % MMP, 11/21/2026 (was)
% and 'cco' have no term for it: an error                                   % MMP, 11/21/2026
Kcon = PIE.K;   hasK = isa(Kcon,'opvar') && size(Kcon,1)>0;
% PIEs with an algebraic constraint (K, Kw, Ku with rows) are not supported
if any(cellfun(@(n) isa(PIE.(n),'opvar') && size(PIE.(n),1)>0,{'K','Kw','Ku'})) % DJ, 11/21/2026
    error('PIETOOLS_H2_norm:algebraicConstraint','The H2 norm LPI is not supported for PIEs with algebraic constraints.');
end
% A constraint involving the input (Kw nonzero) is not supported; u is zero
if hasK && isa(PIE.Kw,'opvar') && any(cellfun(@(B) ~isempty(B) && any(cellfun(@nnz,B.params(:))),PIE.Kw.C(:)))
    error('PIETOOLS_H2_norm:constraintInputs','Algebraic constraints involving inputs are currently not supported.');
end
% if hasK && strcmp(form,'c'),    form = 'o';     end                       % MMP, 11/21/2026 (was)
% if hasK && strcmp(form,'cco'),  form = 'oco';   end                       % MMP, 11/21/2026 (was)
if hasK && any(strcmp(form,{'c','cco'}))                                    % MMP, 11/21/2026
    error('PIETOOLS_H2_norm:formK',['PIE.K is nonempty (the state satisfies K x = 0), ', ...
        'for which form ''%s'' has no term; use ''o'' or ''oco''.'],form);  % MMP, 11/21/2026
end                                                                         % MMP, 11/21/2026
qform = any(strcmp(form,{'c','o'}));

% Degree raise of the slacks and the margins
dwbase = st.sop.dw_P;   if qform,    dwbase = st.sop.dw_Q;   end
eppos = tl.eppos;

% The PIE operators, 'opvar' objects, and the identities on w and z
Tm = PIE.T;     Am = PIE.A;
Bw = PIE.B1;    Cz = PIE.C1;
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE; the bound is sqrt(gam) in the coercive forms
[prog,aux] = lpi_solve_loop(@assemble,st,label);
cl = aux.classify;
gam_lpi = cl.obj;
gam = gam_lpi;
if any(strcmp(form,{'cco','oco'})),     gam = sqrt(gam_lpi);    end

% The solved W (and R, Q), as 'opvar' objects
W = [];     R = [];     Q = [];
try
    W = getsol_lpivar(prog,aux.Wm);
    if qform
        R = getsol_lpivar(prog,aux.Rm);
        Q = getsol_lpivar(prog,aux.Qm);
    end
catch ME
    if st.sop.verbose,  fprintf('  [%s] output extraction: %s\n',label,ME.message);  end
end

% Output struct: solver status, slack, SDP shape, dual kernel, and this LPI
info = exec_info(prog,aux,st,label);
info.form = form;   info.gam = gam;  info.gam_lpi = gam_lpi;    info.nv = nv;
info.W = W;     info.R = R;     info.Q = Q;
info.aux = aux;     % the internals of the last pass, for a posteriori checks % MMP, 11/21/2026

% Numerical witness (1-D): the H2 norm of the discretized PIE
info.witness = [];  info.h2_num = NaN;  info.gap = NaN;
if st.sop.witness && nv==1
    Wt = pie_witness(PIE,'h2',st.sop);
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
        %               Wm ('dopvar' R, Q, W; Rm, Qm [] in 'cco', 'oco'), li
        %               (the info of the slack of the negativity
        %               constraint), form.
        % New program; record the equality rows when the dual is kept
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);
        if st.sop.keepdual,     prog.sopeq = {};    end
        % gam as the objective, gam >= 0 except in 'o'
        [prog,gamv] = lpidecvar(prog,'gam');
        if ~strcmp(form,'o'),   prog = lpi_ineq(prog,gamv);     end
        prog = lpisetobj(prog,gamv);
        % The operators, the negativity (Dneg) and positivity (Dpos)
        % constraints of FORM, the support of the slacks and the trace
        Rm = [];    Qm = [];    Wm = [];
        switch form
            case 'c'
                % R = T*Q >= 0, W >= 0 on the space of z
                [sp,dm] = opvar_space_list(Tm,'out');
                [prog,Rm] = tl.storage(prog,Tm,'out');
                [prog,Qm] = lpivar(prog,dm,sp,tl.dom,get_lpivar_degs(Rm));
                prog = lpi_eq(prog,Tm*Qm - Rm);
                [prog,Wm] = wdecl(prog,Cz,'out');
                Dneg = [-(gamv*Iw),  Bw';
                        Bw,          Qm'*Am' + Am*Qm];
                Dpos = [Wm,         Cz*Qm;
                        Qm'*Cz',    Rm];
                like = Rm;  tr = trace_rn(Wm);
            case 'o'
                % R = T'*Q >= 0, W >= 0 on the space of w
                [sp,dm] = opvar_space_list(Tm,'out');
                [prog,Rm] = tl.storage(prog,Tm,'out');
                [prog,Qm] = lpivar(prog,dm,sp,tl.dom,get_lpivar_degs(Rm));
                prog = lpi_eq(prog,Tm'*Qm - Rm);
                [prog,Wm] = wdecl(prog,Bw,'in');
                % The state block, with K'*Z + Z'*K if K is given
                X22 = Qm'*Am + Am'*Qm;
                if hasK,    [prog,X22] = exec_kmult(prog,X22,Kcon,Tm,tl.dom,Rm);    end
                Dneg = [-(gamv*Iz),  Cz;
                        Cz',         X22];
                Dpos = [Wm,         Bw'*Qm;
                        Qm'*Bw,     Rm];
                like = Rm;  tr = trace_rn(Wm);
            case 'cco'
                % W with its margin; no positivity constraint (Dpos empty)
                [prog,Wm] = tl.storage(prog,Tm,'out');
                Wm = add_eppos2d(Wm);
                Dneg = (Am*Wm)*Tm' + Tm*(Wm*Am') + Bw*Bw';
                Dpos = [];  like = Wm;  tr = trace_rn((Cz*Wm)*Cz');
            case 'oco'
                [prog,Wm] = tl.storage(prog,Tm,'out');
                Wm = add_eppos2d(Wm);
                Dneg = (Am'*Wm)*Tm + Tm'*(Wm*Am) + Cz'*Cz;
                % Acts on the state: K'*Z + Z'*K if K is given
                if hasK,    [prog,Dneg] = exec_kmult(prog,Dneg,Kcon,Tm,tl.dom,Wm);  end
                Dpos = [];  like = Wm;  tr = trace_rn((Bw'*Wm)*Bw);
        end
        % Dneg <= 0 and Dpos >= 0: the 'stock' slacks and equalities, or the
        % slacks sized from the storage
        if strcmp(st.sop.slack,'stock')
            [prog,Nneg] = tl.slack(prog,Dneg);
            [prog,tag] = lpi_eq(prog,Nneg + Dneg,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
            if ~isempty(Dpos)
                [prog,Npos] = tl.slack(prog,Dpos);
                prog = lpi_eq(prog,Npos - Dpos,'symmetric');
            end
        else
            o = struct('like',like,'dw',dwbase+dw,'margin',st.sop.margin,'cap',tl.capfor(Dneg));
            [prog,~,li] = lpi_ineq(prog,-Dneg,o);
            if ~isempty(Dpos)
                op = o;     op.cap = tl.capfor(Dpos);   % its own variables
                [prog,~,~] = lpi_ineq(prog,Dpos,op);
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
        %        X:     'opvar' object (C1 for 'c', B1 for 'o');
        %        side:  char, 'out' or 'in'.
        % OUTPUT prog:  'struct', the program with the Gram matrices of W;
        %        Wm:    'dopvar' object, W.
        % 1-D: the degrees {1,[1 1 1],[1 1 1]}; two or more variables: the
        % storage when a space has variables, else the identity basis on R^n
        if nv==1
            [prog,Wm] = poslpivar_legacy(prog,X,[],struct(),side,PIE.dom);
        elseif any(any(X.(['space_' side])))
            [prog,Wm] = tl.storage(prog,X,side);
        else
            [prog,Wm] = poslpivar_legacy(prog,X,[],struct(),side);
        end
    end

    function Wm = add_eppos2d(Wm)
        % WM = ADD_EPPOS2D(WM) adds eppos I on R^n and on the space in every
        % variable to W, in every dimension (1-D: both spaces; 2-D:
        % blkdiag(eppos(1) I, 0, 0, eppos(4) I)).
        % INPUT/OUTPUT  Wm: 'dopvar' object, W.
            % One margin per variable subset, zero on the partial ones, of
            % which 1-D has none
            e = eppos(:);   if isscalar(e),     e = e*ones(2^nv,1);     end
            e(2:end-1) = 0;
            if ~all(e==0),  Wm = Wm + tl.eye(1,e);  end
    end
end
