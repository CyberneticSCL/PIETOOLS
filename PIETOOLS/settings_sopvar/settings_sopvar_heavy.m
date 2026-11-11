function st = settings_sopvar_heavy()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = SETTINGS_SOPVAR_HEAVY() returns the settings of the 'heavy' tier of
% LPISETTINGS_SOP: the degrees and options from which the container
% executives build the positive operator P of an LPI, the margins of P, and
% the options of the executives. The degrees are those of the stock heavy
% level in one spatial variable and of the stock stripped level in two.
%
% INPUT
% none.
%
% OUTPUT
% - st:   'struct' with fields
%           override1, options1, options12, dd1, dd12   P on a PIE in one
%                         spatial variable;
%           ddZ           1x3 'double', the free operator of the synthesis
%                         LPIs in one spatial variable;
%           sosineq_on, override2, options2, dd2, options3, dd3 (or opts)
%                         the slack of the 'stock' slack mode in one
%                         spatial variable;
%           settings_2d   'struct', the same in two spatial variables:
%                         is2D, LF_deg, LF_opts, LF_use_psatz,
%                         LF_opts_psatz, LF_deg_psatz (P); Zop_deg,
%                         Zop_opts (free operator); use_sosineq and
%                         ineq_opts or eq_deg, eq_opts, eq_use_psatz,
%                         eq_opts_psatz, eq_deg_psatz ('stock' slack);
%           eppos         'double', the margin of P on every space, in every
%                         dimension;
%           epneg         'double', the decay rate, in every dimension;
%           sos_opts      'struct' with field simplify (0 or 1);
%           settings_nd   'struct' with fields LF_int, LF_mult: P on a PIE in
%                         three or more spatial variables;
%           sop           'struct', the options of the container executives:
%                         slack, dw_feas, dw_P, dw_Q, max_raise_feas,
%                         max_raise, cap_delta, margin, keepdual,
%                         witness, N_cheb, nfreq;
%           lpi           'struct' with field ineq, the degree options of
%                         LPI_INEQ_SOP: dD, dw, wmin.
%         The body gives the value of each field and what it does.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - settings_sopvar_heavy
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
% MMP, 10/29/2026: 3-D entries: settings_nd, and a third entry of dw_feas,
%                max_raise_feas, max_raise and cap_delta: weight 2, no
%                raise, the per-direction bound Jmax(T) + 1.
% MMP, 11/09/2026: The storage values are in this file and the body is
%                reorganized by the code that reads each field (maintainer:
%                "Files should not be simply containers for other files",
%                "Re-organize them so it is clear what is going on"). The
%                values of the stock heavy (1-D) and stripped (2-D) levels
%                (M. Peet, S. Shivakumar, D. Jagt, 2021-2022), called from
%                settings/ until 11/08/2026 and from copies in
%                settings_sopvar/private (removed) on 11/08/2026, are
%                written out. Deleted: the commented-out st.sop of
%                10/29/2026 and the commented-out calls of
%                settings_PIETOOLS_heavy and settings_PIETOOLS_stripped_2D
%                of 11/08/2026. The returned structure is byte-identical,
%                with and without sosineq_on or use_sosineq in the base
%                workspace. The 11/08/2026 entry is removed.
% MMP, 11/09/2026: One margin eppos and one decay rate epneg for every
%                dimension (maintainer: "I see no reason to treat the eppos
%                and epneg options differently for 1D,2D and ND"): eppos2,
%                settings_2d.eppos and epneg, settings_nd.eppos and
%                sop.eppos_stab are removed, and EXEC_SOP_SETTINGS sets the
%                fields each dimension reads from eppos and epneg. eppos is
%                1e-4, the former 1-D value: the margin in two or more
%                variables, and in stability and well-posedness, falls from
%                1e-2 (README, measurements of 11/09/2026).
% MMP, 11/10/2026: st.lpi.ineq, the degree options of LPI_INEQ_SOP (dD 1,
%                dw 1, wmin [2 2 0]; maintainer, 11/10/2026: "settings files
%                should have the knobs used in lpi_programming, especially
%                lpi_ineq").

% BEGIN MMP, 11/09/2026: body reorganized, one margin and decay rate (header)

%% 1. PIEs in one spatial variable s in [a,b]
% P = P1 + P2, each Pi = Zi'*Qi*Zi with Qi >= 0 (POSLPIVAR_SETTINGS_SOP,
% role 'lf'). Degrees dd = {d1, d2, d3}: d1 of the multiplier basis; d2, d3
% of the lower- and upper-integral bases, each [degree in the integration
% variable, in s, joint]. Options: sep 1 takes one integral over [a,b] for
% the lower and upper ones; exclude (1x4) drops the R^n, multiplier, lower-
% and upper-integral bases; psatz 1 weights Q by (s-a)(b-s).
st.override1 = 1;                                   % 1: P = P1, dd12 and options12 unread
st.options1  = struct('sep',0,'exclude',[0 0 0 0]);
st.options12 = struct('sep',0,'exclude',[0 0 0 0],'psatz',1);
st.dd1  = {2, [1 1 2], [1 1 2]};
st.dd12 = {2, [1 1 2], [1 1 2]};
% Free operator of the synthesis LPIs (LPIVAR_CDOPVAR): [d1 d2 d3], d1 the
% degree of the multiplier and of the maps to and from R^n, d3 and d2 the
% degrees of the integral kernels in s and in the integration variable.
st.ddZ = [4 2 2];
% Slack N = N1 + N2 of the 'stock' slack mode (role 'slack'), read only
% with sop.slack = 'stock': N1 of dd2 and options2, N2 of dd3 and options3,
% N2 omitted when override2 = 1. sosineq_on is read from the base workspace
% when it exists there, as in the stock settings files: 1 selects the stock
% sosineq, which the container executives reject, and sets only opts, which
% none of them reads.
st.sosineq_on = 0;
if evalin('base','exist(''sosineq_on'',''var'')')
    st.sosineq_on = evalin('base','sosineq_on');
end
if st.sosineq_on
    st.opts = struct('psatz',1,'pure',1);
else
    st.override2 = 0;                               % 0: N = N1 + N2
    st.options2  = struct('psatz',0,'exclude',[0 0 0 0]);
    st.dd2  = {3, [2 2 3], [2 2 3]};
    st.options3  = struct('psatz',1,'exclude',[0 0 0 0]);
    st.dd3  = {2, [1 1 2], [1 1 2]};
end

%% 2. PIEs in two spatial variables (x,y): st.settings_2d
% P = P0 + P1 + ..., each term Z'*Q*Z with Q >= 0
% (POSLPIVAR_SETTINGS_2D_SOP, role 'lf'), Z built from the 16 bases
% [R^n; Zx^o, Zx^a, Zx^b; Zy^o, Zy^a, Zy^b; Z2^oo, Z2^ao, Z2^bo, Z2^oa,
% Z2^ob, Z2^aa, Z2^ba, Z2^ab, Z2^bb]: o the multiplier, a the lower and b
% the upper integral, in x for the first index and in y for the second;
% tt and nu the integration variables in x and y.
% Degrees: dx on L2[x], dy on L2[y], d2 on L2[x,y]. The caps of Z2^{ij}
% have rows [-; x] (i = o) or [-; x; tt; x*tt] (i = a, b) and columns
% [-, y] (j = o) or [-, y, nu, y*nu] (j = a, b), each entry the cap on the
% joint degree in its row and column variables.
dx = {1; [1;1;1]; [1;1;1]};                         % Zx^o: degree in x; Zx^a, Zx^b: caps [x; tt; x*tt]
dy = {1, [1 1 1], [1 1 1]};                         % the same in y and nu
d2 = cell(3,3);
d2{1,1} = [0 1;1 2];                                % Z2^oo
d2{1,2} = [0 1 1 1;1 2 2 2];                        % Z2^oa
d2{1,3} = [0 1 1 1;1 2 2 2];                        % Z2^ob
d2{2,1} = [0 1;1 2;1 2;1 2];                        % Z2^ao
d2{2,2} = [0 1 1 1;1 2 2 2;1 2 2 2;1 2 2 2];        % Z2^aa
d2{2,3} = [0 1 1 1;1 2 2 2;1 2 2 2;1 2 2 2];        % Z2^ab
d2{3,1} = [0 1;1 2;1 2;1 2];                        % Z2^bo
d2{3,2} = [0 1 1 1;1 2 2 2;1 2 2 2;1 2 2 2];        % Z2^ba
d2{3,3} = [0 1 1 1;1 2 2 2;1 2 2 2;1 2 2 2];        % Z2^bb
% Options: exclude (1x16) drops the bases in the order above; sep (6x1)
% equates the lower- and upper-integral parameters of P: Rxx{2} = Rxx{3}
% (1), Ryy{2} = Ryy{3} (2), R22{2,1} = R22{3,1} (3), R22{1,2} = R22{1,3}
% (4), R22{2,j} = R22{3,j} for j = 2,3 (5), R22{i,2} = R22{i,3} for
% i = 2,3 (6), with R22{i,j} the parameter of Z2^{ij}; psatz 1 weights Q
% by the domain, 2 by a disk (rejected by POSCOPVAR). Each nonzero entry
% of LF_use_psatz adds a term of the options LF_opts_psatz{j} and the
% degrees LF_deg_psatz{j}.
s2.is2D = 1;
s2.LF_deg  = struct('dx',{dx},'dy',{dy},'d2',{d2});
s2.LF_opts = struct('sep',[1;1;1;1;1;1],'exclude',[0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0],'psatz',0);
s2.LF_use_psatz  = 0;                               % no added term: the two fields below unread
s2.LF_opts_psatz = {struct('psatz',1,'sep',s2.LF_opts.sep,'exclude',s2.LF_opts.exclude), ...
                    struct('psatz',2,'sep',s2.LF_opts.sep,'exclude',s2.LF_opts.exclude)};
s2.LF_deg_psatz  = {s2.LF_deg, s2.LF_deg};
% Free operator of the synthesis LPIs: twice the degrees of P, of which only
% the largest entry is read (EXEC_TOOLS_SOP); Zop_opts is read by no
% container executive.
twice = @(c) cellfun(@(d) 2*d, c, 'UniformOutput', false);
s2.Zop_deg  = struct('dx',{twice(dx)},'dy',{twice(dy)},'d2',{twice(d2)});
s2.Zop_opts = struct('sep',[0;0;0;0;0],'isscalar',0,'ismultiplier',1);
% Slack of the 'stock' slack mode (role 'slack'), read only with sop.slack
% = 'stock'. use_sosineq is read from the base workspace (sosineq_on, else
% use_sosineq) when it exists there, as in the stock settings files: 1
% selects the stock sosineq, which the container executives reject, and
% sets only ineq_opts, which none of them reads. With 0, the value here, the
% slack has the degrees of P raised by 2 in dx and dy and by 3 in d2, and
% one term per nonzero entry of eq_use_psatz (none here) with those raised
% by 1 and by 2.
s2.use_sosineq = 0;
if evalin('base','exist(''sosineq_on'',''var'')')
    s2.use_sosineq = evalin('base','sosineq_on');
elseif evalin('base','exist(''use_sosineq'',''var'')')
    s2.use_sosineq = evalin('base','use_sosineq');
end
if s2.use_sosineq
    s2.ineq_opts = struct('psatz',0);
else
    up = @(c,k) cellfun(@(d) d+k, c, 'UniformOutput', false);
    s2.eq_deg  = struct('dx',{up(dx,2)},'dy',{up(dy,2)},'d2',{up(d2,3)});
    s2.eq_opts = struct('psatz',0,'sep',zeros(1,6),'exclude',zeros(1,16));
    s2.eq_use_psatz = [0;0];
    for j = 1:numel(s2.eq_use_psatz)
        s2.eq_opts_psatz{j} = struct('psatz',s2.eq_use_psatz(j),'sep',zeros(1,6),'exclude',zeros(1,16));
        s2.eq_deg_psatz{j}  = struct('dx',{up(dx,1)},'dy',{up(dy,1)},'d2',{up(d2,2)});
    end
end
st.settings_2d = s2;

%% 3. Margin, decay rate and solver, the same in every dimension
% eppos: the margin of P, P >= eppos*I on every space (R^n and each L2), in
% each LPI that imposes one, stability and well-posedness included; epneg:
% the decay rate. EXEC_SOP_SETTINGS sets from these the fields each
% dimension reads (eppos2, settings_2d.eppos and epneg, settings_nd.eppos);
% such a field, set here, takes precedence. LPISETTINGS_SOP overwrites
% epneg and simplify with its inputs, and the solver when it is given one.
st.eppos = 1e-4;    st.epneg = 0;
st.sos_opts = struct('simplify',0);

%% 4. PIEs in three or more spatial variables: st.settings_nd
% P on each space (EXEC_TOOLS_SOP): a multiplier basis of degree LF_mult
% and one full-integral basis, separable in every direction, of degree
% LF_int in each variable and 2*LF_int in total; the margin and decay rate
% of section 3.
st.settings_nd = struct('LF_int',1,'LF_mult',0);

%% 5. Options of the container executives: st.sop (EXEC_SOP_SETTINGS)
% A vector holds the entries for 1, 2, and 3 or more spatial variables; a
% scalar serves every dimension. The slack of each negativity constraint
% -D >= 0 has a weight degree, raised by dw over the rule of GET_LIFT_DEGS,
% and the degree loop LPI_SOLVE_LOOP_SOP rebuilds the LPI with the weight
% one higher after UNKNOWN or INFEASIBLE, up to max_raise times.
st.sop = struct( ...
    'slack','like', ...             'like': slack sized from P (LPI_INEQ_SOP); 'stock': the slacks above
    'dw_feas',[1 2 2], ...          weight raise in stability and well-posedness
    'dw_P',2, ...                   weight raise in coercive gains and norms, and synthesis
    'dw_Q',1, ...                   weight raise in the Q (non-coercive) forms
    'max_raise_feas',[2 2 0], ...   raises in stability and well-posedness
    'max_raise',[1 1 0], ...        raises in every other LPI
    'cap_delta',[NaN NaN 1], ...    slack bases bounded per direction by Jmax(T) + cap_delta; NaN: no bound
    'margin',0, ...                 eps > 0 certifies -D + eps*e_D >= 0 in place of -D >= 0
    'keepdual',true, ...            assemble the dual kernel of the negativity constraint
    'witness',true, ...             one variable: gain, decay rate or H2 norm of a Legendre model
    'N_cheb',24, ...                its Legendre degree
    'nfreq',240);                   % its number of frequencies
% END MMP, 11/09/2026

%% 6. Options of LPI_INEQ_SOP: st.lpi.ineq                                  % MMP, 11/10/2026
% The certificate of P >= 0 in a hand-written LPI: lift degree Dmin + dD    % MMP, 11/10/2026
% and weight degree reach + dw over the support of P (GET_LIFT_DEGS), at    % MMP, 11/10/2026
% least wmin when read from P itself (no opts.like). wmin 2: with dw alone  % MMP, 11/10/2026
% the Poincare LPI of DEMO3 is 22% above 1/pi, at w = 2 within 2e-5         % MMP, 11/10/2026
% (README). A vector holds the entries for 1, 2, and 3 or more variables;   % MMP, 11/10/2026
% a call that sets dw and not wmin gets no floor (LPI_INEQ_SOP).            % MMP, 11/10/2026
st.lpi.ineq = struct('dD',1,'dw',1,'wmin',[2 2 0]);                         % MMP, 11/10/2026
end
