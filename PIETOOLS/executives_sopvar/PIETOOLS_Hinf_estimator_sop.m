function [prog,L,gam,P,Z,info] = PIETOOLS_Hinf_estimator_sop(PIE,settings,gain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,L,GAM,P,Z,INFO] = PIETOOLS_HINF_ESTIMATOR_SOP(PIE,SETTINGS,GAIN)
% computes an observer gain L for the PIE
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
% with P a positive operator, Z a free operator, and, in one or two spatial
% variables, L the gain of bounded degree that minimizes ||P*L - Z||^2, an
% approximation of P^{-1}*Z (in two variables, should that computation
% fail, L = P^{-1}*Z for a separable P). While the solver returns UNKNOWN or
% INFEASIBLE, the program is rebuilt with the degree of the negativity
% slack raised by one, up to settings.sop.max_raise times. PIEs with an
% algebraic constraint (K, Kw or Ku with rows, e.g. from periodic boundary
% conditions) are not supported.
%
% INPUT
% - PIE:      'pie_struct' object of a PIE in any number of spatial
%             variables, initialized when its T is not yet a 'copvar'
%             object; T, Tw, A, B1, C1, D11, C2 and D21 'copvar' objects (Tw
%             may be absent), and K, Kw, Ku empty or without rows; its
%             other operators are not used; in two or more variables, Tw
%             and B1 not both nonzero; C2 with at least one row (a sensed
%             output), else an error;
% - settings: (optional) 'struct' of LPI settings:
%               lpi.pos   the storage row of P, in every dimension;
%               lpi.var   field deg, a 1x3 'cell': deg{min(nv,3)} the
%                         degrees of Z for a PIE in nv variables;
%               settings_2d  (optional) the 16-basis 2-D degrees of P and of
%                         Z (the largest entry of Zop_deg), in place of
%                         lpi.pos and lpi.var in two variables;
%               eppos     the margin of P on every space, in every
%                         dimension;
%               sos_opts  solver options;
%               sop       the options of the executives, among them dw_P,
%                         the degree raise of the slack, max_raise and
%                         gam_fixed (as GAIN).
%             Default: the 'light' tier;
% - gain:     (optional) positive double: solve the feasibility LPI at
%             gam = gain instead of minimizing gam (0 or [] minimizes).
%
% OUTPUT
% - prog:     'struct', the solved LPI program;
% - L:        'copvar' object, the observer gain from
%             y to the state ([] when the solver returned no solution, and
%             in 3 or more variables, where no gain extraction exists: the
%             LPI is solved and P, Z are returned);
% - gam:      double, the certified bound; with GAIN, GAIN when certified and
%             Inf otherwise; NaN when the solver returned no value;
% - P, Z:     'copvar' objects, the solved P and Z;
% - info:     'struct' with fields
%               status     char: 'optimal', 'inaccurate', 'infeasible',
%                          'unbounded', 'unknown' or 'residual' (cuADMM,
%                          which decides no feasibility: never certified);
%               certified  logical, true when status is 'optimal' or
%                          'inaccurate';
%               form ('est'), gam, gam_fixed, nv;
%               P, Z       as the outputs; Psol, Zsol the same objects;
%               K          'cdopvar' object, the LPI operator;
%               G          'copvar' object, the gain L (one or two
%                          variables; [] otherwise); inv  'struct', the
%                          diagnostics of the computation of L (output INFO
%                          of GETOBSERVER_SOP: for the method 'qp' the field
%                          qp, no inverse of P being formed; for the 2-D
%                          fallback 'inv' the data of the inverse);
%               dw, hist, slack, degrees, terms, shape, dual, dual_norm1
%                          the last raise and the loop history, the slack
%                          mode, its degrees and terms, the size of the SDP,
%                          and the dual kernel of the negativity constraint
%                          with its 1-norm.
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
% DJ, 11/21/2026: Error for a PIE with an algebraic constraint (K, Kw, Ku
%                with rows);
% MMP, 11/24/2026: Error for a PIE without sensed output (C2 without rows):
%                no estimation is possible; the LPI was built and solved,
%                then the gain extraction failed ('Index exceeds array
%                bounds'), and L was returned empty. Help: C2 needs a row.

% Default settings and fixed gain; check the PIE
if nargin<2 || isempty(settings),  settings = lpisettings_sop('light'); end
if nargin>=3 && ~isempty(gain) && gain~=0,  settings.sop.gam_fixed = gain;  end
if ~isa(PIE,'pie_struct')
    error('PIETOOLS_Hinf_estimator_sop:pie','The PIE should be a ''pie_struct'' object.');
end
PIE = exec_pie_sop(PIE);                    % every operator a 'copvar' or 0
label = 'Hinf_estimator_sop';

% Complete the options of the executives and read the tools of the
% dimension of the PIE (storage of P, margins, degrees of Z); in two or
% more variables, w and dw/dt together are an error
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
st.sop = tl.sop;    % the options of this dimension
Tw0 = (PIE.Tw==0);
if nv>=2 && ~Tw0 && ~(PIE.B1==0)
    error('PIETOOLS_Hinf_estimator_sop:TwB','The PIE takes both the input w and its derivative dw/dt; LPI based estimator synthesis is currently not supported.');
end
gfix = [];  if isfield(st.sop,'gam_fixed') && ~isempty(st.sop.gam_fixed) && st.sop.gam_fixed~=0,    gfix = st.sop.gam_fixed;    end
eppos = tl.eppos;

% The PIE operators, 'copvar' objects, and the identities on w and z
Tm = PIE.T;     Am = PIE.A;
% The constraint K x = 0 on the state, if any (empty for most PIEs)
Kcon = PIE.K;   hasK = isa(Kcon,'copvar') && size(Kcon,1)>0;
% PIEs with an algebraic constraint (K, Kw, Ku with rows) are not supported
if any(cellfun(@(n) isa(PIE.(n),'copvar') && size(PIE.(n),1)>0,{'K','Kw','Ku'})) % DJ, 11/21/2026
    error('PIETOOLS_Hinf_estimator_sop:algebraicConstraint','The H-infinity estimator synthesis LPI is not supported for PIEs with algebraic constraints.');
end
% A constraint involving the input (Kw nonzero) is not supported; u is zero
if hasK && isa(PIE.Kw,'copvar') && any(cellfun(@(B) ~isempty(B) && any(cellfun(@nnz,B.params(:))),PIE.Kw.C(:)))
    error('PIETOOLS_Hinf_estimator_sop:constraintInputs','Algebraic constraints involving inputs are currently not supported.');
end
Bw = PIE.B1;    Cz = PIE.C1;    Dzw = PIE.D11;
Cy = PIE.C2;    Dyw = PIE.D21;
% No sensed output y (C2 without rows): no estimation is possible           % MMP, 11/24/2026
if size(Cy,1)==0                                                            % MMP, 11/24/2026
    error('PIETOOLS_Hinf_estimator_sop:noOutput','The PIE has no sensed outputs (C2 has no rows); estimation is not possible.'); % MMP, 11/24/2026
end                                                                         % MMP, 11/24/2026
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);
Twm = [];   if ~Tw0,    Twm = PIE.Tw;    end

% Build and solve, raising the slack degree while the solver returns
% UNKNOWN or INFEASIBLE
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
certified = any(strcmp(cl.status,{'optimal','inaccurate'}));
if isempty(gfix),   gam = cl.obj;
elseif certified,   gam = gfix;
else,               gam = Inf;
end

% The solved P and Z, as 'copvar' objects, and the gain L, an approximation
% of P^{-1}*Z, in one or two variables
L = [];     P = [];     Z = [];     Psol = [];  Zsol = [];  Gc = [];    iinv = [];
try
    Psol = getsol_lpivar_sop(prog,aux.Pm);  Zsol = getsol_lpivar_sop(prog,aux.Zm);
    % L of bounded degree minimizing ||P*L - Z||^2
    P = Psol;      Z = Zsol;
    if nv==1,   [Gc,~,iinv] = getObserver_sop(Psol,Zsol);   L = Gc;
    elseif nv==2
        % 2-D: the gain of the default method 'qp'; after an error, the
        % other method, 'inv' (the inverse of a separable P)
        try
            [Gc,~,iinv] = getObserver_sop(Psol,Zsol);
        catch ME2
            if st.sop.verbose,  fprintf('  [%s] gain extraction: %s; inverse (inv).\n',label,ME2.message);  end
            [Gc,~,iinv] = getObserver_sop(Psol,Zsol,struct('method','inv'));
        end
        L = Gc;
    elseif st.sop.verbose                   % no N-D gain extraction
        fprintf('  [%s] gain extraction: none in %d variables; P and Z returned.\n',label,nv);
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
        %               operator), Nm (the slack), li (the info of the
        %               slack), form.
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
        % The state block, with K'*Z_K + Z_K'*K if K is given
        X33 = PA'*Tm + Tm'*PA;
        if hasK,    [prog,X33] = exec_kmult_sop(prog,X33,Kcon,Tm,tl.dom,Pm);    end
        Km = [K11,      -Dzw',       K13;
              -Dzw,     -(gamv*Iz),  Cz;
              K31,      Cz',         X33];
        % Km <= 0: the 'stock' slack and an equality, or the slack sized from P
        if strcmp(st.sop.slack,'stock')
            [prog,Nm] = tl.slack(prog,Km);
            [prog,tag] = lpi_eq_sop(prog,Nm+Km,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
            o = struct('like',Pm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin,'cap',tl.capfor(Km));
            [prog,Nm,li] = lpi_ineq_sop(prog,-Km,o);
        end
        % The objective gam, none at a fixed gain
        obj = gamv;     if ~isempty(gfix),   obj = [];  end
        aux = struct('obj',obj,'Pm',Pm,'Zm',Zm,'Km',Km,'Nm',Nm,'li',li,'form','est');
    end
end
