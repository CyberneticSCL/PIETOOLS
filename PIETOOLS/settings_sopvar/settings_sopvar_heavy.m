function st = settings_sopvar_heavy()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = SETTINGS_SOPVAR_HEAVY() returns the settings of the 'heavy' tier of
% LPISETTINGS_SOP: the degrees and options from which the container
% executives build the positive operator P of an LPI, the margins of P, and
% the options of the executives. The storage is a row (POSCOPVAR, DEG       % MMP, 11/17/2026
% omitted), in one variable the stock heavy storage.                        % MMP, 11/17/2026
%
% INPUT
% none.
%
% OUTPUT
% - st:   'struct' with fields (settings report 11/16/2026, Sec. 5, 6):     % MMP, 11/17/2026
%           lpi       'struct' with fields pos, the storage row that        % MMP, 11/17/2026
%                     POSCOPVAR reads with DEG omitted (len_R, len_inc,     % MMP, 11/17/2026
%                     len_com, len_sos, patterns, cap, len_cross, terms,    % MMP, 11/17/2026
%                     offsets); ineq, the slack options of LPI_INEQ_SOP     % MMP, 11/17/2026
%                     (dD, dw, wmin, wcross, psatz, psatz_offset); var,     % MMP, 11/17/2026
%                     field deg, the free operator of LPIVAR_CDOPVAR;       % MMP, 11/17/2026
%           sop       'struct', the options of the container executives:    % MMP, 11/17/2026
%                     slack, dw_feas, dw_P, dw_Q, max_raise_feas,           % MMP, 11/17/2026
%                     max_raise, cap_delta, margin;                         % MMP, 11/17/2026
%           eppos     'double', the margin of P on every space, in every    % MMP, 11/17/2026
%                     dimension;                                            % MMP, 11/17/2026
%           epneg     'double', the decay rate, in every dimension;         % MMP, 11/17/2026
%           sos_opts  'struct' with fields simplify (0 or 1) and solver.    % MMP, 11/17/2026
%         The body gives the value of each field and what it does.          % MMP, 11/17/2026
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
% MMP, 11/12/2026: The positive and free operators in native form, lpi.pos
%                and lpi.var (maintainer, 11/12/2026: phase 2 of the knob
%                reduction): they replace the stock-format fields of
%                section 1 (override1, options1, options12, dd1, dd12, ddZ)
%                and settings_nd (section 4), with the same programs. The
%                'stock' slack (sosineq_on, override2, options2, dd2,
%                options3, dd3) leaves for LPISETTINGS_STOCK_SOP
%                (maintainer: "retire it, but create a utility function
%                which creates an ad hoc settings file which reproduces the
%                baselines"); the slack comment of sop names it. Sections 5
%                and 6 are now 4 and 5.
% MMP, 11/17/2026: The structure of the settings report 11/16/2026, Sec. 5
%                (A tier layout) and Sec. 6 (Tables of the storage rows
%                and of the slack and raise knobs; work package 5): lpi.pos
%                is the storage row that POSCOPVAR reads with DEG omitted
%                in every dimension; lpi.ineq gains wcross 0, psatz 'auto'
%                and psatz_offset [], and sos_opts the field solver ''.
%                lpi.var.deg{2} 2 (was 1): the direct default test certified
%                S17 at 0.5 r2* with degree 2 and not 1 (report, Sec. 6).
%                Deleted, described in the report, Sec. 7 (the present
%                tiers): settings_2d (the stock 2-D storage LF_*, free
%                operator Zop_* and 'stock' slack); the lpi.pos fields mult,
%                int, joint, sep; the sop fields keepdual, witness, N_cheb,
%                nfreq, which govern the output and not the program (now
%                defaults of EXEC_SOP_SETTINGS and LPISETTINGS_SOP). The
%                body of the entries of 11/09/2026, 11/10/2026 and
%                11/12/2026 is replaced; their values stand only where the
%                report keeps them.
% MMP, 11/17/2026: lpi.ineq.psatz {'auto','faces1','auto'} (was 'auto';
%                maintainer decision, work package 5): in two variables
%                the plain term and the four faces one weight degree below
%                it whatever w ('faces1', GET_LIFT_DEGS); 'auto' in one and
%                in three or more. With wcross 0 and 'auto' (faces at
%                offset 0 when w = 1) S17 at 0.5 r2* in the direct form is
%                UNKNOWN (676,410 variables after one raise; terms_results.md,
%                "The cross length"), where the faces at offset 1 keep the
%                certificate; on every 2-D case measured on 11/17/2026
%                (gap_results.md G1-G4; S1, S3, S4, S8, PDE 33) the offset-1
%                faces certify whatever offset 0 certifies at about a
%                quarter of the variables (S3 95,930 against 383,270; S8
%                123,168 against 492,048; PDE 33 52,660 against 190,900).
% MMP, 11/18/2026: Three fields of section 4 for the direct stability form
%                'Q' (PIETOOLS_STABILITY_SOP; maintainer: "make the
%                separation"): deg_Q, the degrees of its Q (this tier's
%                lpi.var.deg, which stays the degree of the synthesis free
%                operators), max_raise_Q, its raises (this tier's
%                max_raise_feas), and raise_Q 'slack', what a raise adds
%                (the slack weight, as before; 'deg' the first entry of
%                deg_Q; 'both'). Values of 11/17/2026: programs unchanged.
% MMP, 11/18/2026: lpi.ineq gains Dmax false and cap_inc [] (GET_LIFT_DEGS,
%                same date): the lift degree at the support ceiling and the
%                triangular cap of Shivakumar's bounds (11/17/2026), both
%                off; programs unchanged.

% BEGIN MMP, 11/17/2026: body replaced by the tier layout of the settings
% report 11/16/2026, Sec. 5, values of Sec. 6 (header). Deleted, stamped
% MMP 11/09/2026, 11/10/2026, 11/12/2026: sections 1 (lpi.pos mult, int,
% joint, sep; lpi.var), 2 (settings_2d), 3 (eppos, epneg, sos_opts), 4
% (st.sop with keepdual, witness, N_cheb, nfreq) and 5 (lpi.ineq dD, dw,
% wmin); the values of 1, 2 and 4 are in the report, Sec. 7.
% Heavy tier. Entry k (cell) or column k (row): a space or inequality in
% k variables, the last entry for every higher k.
% 1. The storage (POSCOPVAR, DEG omitted): lengths of monomial vectors,
% 0 = absent; d1, d2, d3 of the report, Table of the storage rows.
st.lpi.pos = struct( ...
   'len_R',    1, ...                                   % R^q: the identity
   'len_inc',  {{2,       [0 0],         [0 0 0]}}, ...             % d1
   'len_com',  {{2,       [2 2],         [2 2 2]}}, ...             % d2
   'len_sos',  {{[3;2;2], [2 2;2 2;0 0], [1 1 1;2 2 2;0 0 0]}}, ... % d3
   'patterns', {{'all',   'separable',   'uniform'}}, ...           % Pi
   'cap',      {{Inf,     [Inf Inf],     [Inf Inf Inf]}}, ...       % c
   'len_cross', {{1,      1,             1}}, ...        % non-owned, 1: none
   'terms',    {{'none',  'none',        'none'}}, ...              % F
   'offsets',  {{[],      [],            []}});                     % o
% 2. The slack (LPI_INEQ_SOP, GET_LIFT_DEGS); degrees.
st.lpi.ineq = struct( ...
   'dD',           1, ...         % lift raise: D = Dmin + dD
   'dw',           1, ...         % weight raise, calls without dw
   'wmin',         [2 2 0], ...   % weight floor, target the storage
   'wcross',       0, ...         % weight, non-owned dirs; 0: none
   'psatz',        {{'auto','faces1','auto'}}, ... % F: 2-D faces1
   'Dmax',         false, ...     % lift at the support ceiling; off        % MMP, 11/18/2026
   'cap_inc',      [], ...        % triangular cap; 'auto' = Dmax-1; off    % MMP, 11/18/2026
   'psatz_offset', []);           % o: read with a row of codes only        % MMP, 11/17/2026
% st.lpi.ineq.psatz = 'auto';  the rule in every dimension                  % MMP, 11/17/2026 (was)
% 3. The free operator (LPIVAR_CDOPVAR); degrees.
st.lpi.var = struct('deg', {{[4 2 2], 2, 1}});  % [d1 d2 d3]; caps
% 4. The executives (EXEC_SOP_SETTINGS).
st.sop = struct( ...
   'slack',          'like', ...      % slack sized from the storage
   'dw_feas',        [1 2 2], ...     % weight raise: stability, w.-p.
   'dw_P',           2, ...           % weight raise: gains, synthesis
   'dw_Q',           1, ...           % weight raise: the Q forms
   'max_raise_feas', [2 2 0], ...     % raises: stability, w.-p.
   'max_raise',      [1 1 0], ...     % raises: every other LPI
   'cap_delta',      [NaN NaN 1], ... % Jmax(T) + cap_delta; NaN: none
   'deg_Q',          {{[4 2 2], 2, 1}}, ... % Q, direct stability form      % MMP, 11/18/2026
   'max_raise_Q',    [2 2 0], ...     % its raises                          % MMP, 11/18/2026
   'raise_Q',        'slack', ...     % a raise adds: slack weight          % MMP, 11/18/2026
   'margin',         0);           % eps of -D + eps*e_D >= 0
% 5. Tolerances and solver.
st.eppos = 1e-4;                      % P >= eppos*I on every space
st.epneg = 0;                         % decay rate
st.sos_opts = struct( ...
   'simplify', 0, ...                 % 1: simplify the SDP
   'solver',   '');                   % '': MOSEK if found, else SeDuMi
% END MMP, 11/17/2026
end
