function prog = lpiprogram_sop(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPIPROGRAM_SOP(VARTAB,DUMVARTAB,DOM,DECVARTAB,FREEVARTAB) declares
% an LPI program structure in any number of spatial variables. It is
% 'lpiprogram' (which since 10/01/2026 takes any number of variables) with
% two more ways to name the variables, so that one program serves the
% legacy operator classes and the containers of the 'sopvar' family:
%
%   PROG = LPIPROGRAM_SOP(VARTAB,DUMVARTAB,DOM[,DECVARTAB[,FREEVARTAB]])
%   PROG = LPIPROGRAM_SOP(VARTAB,DOM[,DECVARTAB[,FREEVARTAB]])
%          every form 'lpiprogram' accepts, with the same result;
%   PROG = LPIPROGRAM_SOP(NAMES,DOM,...), NAMES a cellstr (the registry
%          P.vars of a container, with DOM = P.dom);
%   PROG = LPIPROGRAM_SOP(P[,DECVARTAB[,FREEVARTAB]]), P a 'copvar',
%          'cdopvar', 'sopvar' or 'sdopvar': the variables and domains of
%          P, in the order of its sorted registry;
%   PROG = LPIPROGRAM_SOP(S,...), S a struct with fields 'vars' (cellstr)
%          and 'dom' (the struct form of DOM in 'lpivar_cdopvar').
%
% INPUT
% - vartab:     n x 1 'polynomial' array of the spatial variables, or a
%               cellstr of their names;
% - dumvartab:  (optional) n x 1 'polynomial' array (or cellstr) of dummy
%               variables, one per spatial variable; [] gives s_dum for
%               variable s, as 'lpiprogram' does;
% - dom:        n x 2 'double', row i the interval of vartab(i); one 1 x 2
%               row applies to every variable;
% - decvartab:  (optional) q x 1 'dpvar' of decision variables;
% - freevartab: (optional) m x 1 'polynomial' of non-spatial variables.
%
% OUTPUT
% - prog:       SOSTOOLS program structure with the field 'dom' added, the
%               fields and values 'lpiprogram' builds: vartable =
%               [vartab; dumvartab; freevartab] ('polynomial'), dom n x 2,
%               decvartable a cellstr. It is the output of 'lpiprogram' on
%               the same variables.
%
% NOTES
% The program is built by 'lpiprogram'; this routine only turns the new
% input forms into its arguments. The container routines (lpivar_cdopvar, poscopvar, lpi_eq_cdopvar,
% lpi_eq_sdopvar, getsol_lpivar_sop) read neither prog.vartable nor
% prog.dom; 'lpidecvar', 'lpisetobj' and 'lpisolve' read neither either.
% The legacy 'lpi_eq' and 'lpi_ineq' read the dummy variables at
% prog.vartable(n+1:2n), which is where they are here, as in 'lpiprogram'.
% The legacy 1-D and 2-D operator routines (lpivar, poslpivar) still
% support at most 2 variables and check that themselves.
%
% Cost: that of 'lpiprogram', plus O(n) for the name conversion.
%
% See also LPIPROGRAM, LPIDECVAR, LPI_EQ_SOP, LPIGETSOL_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpiprogram_sop
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
% Initial coding MMP, 09/29/2026. 'lpiprogram' (DJ, 10/13/2024; 11/30/2024;
%                01/23/2025) refuses more than 2 spatial variables, and the
%                N-D container scripts (heatNd_lpi, heatNd_poincare,
%                test_copquadvar_faces) rebuilt its output by hand for
%                N > 2. This follows lpiprogram's parsing and construction
%                in full, without the cap, so one code path serves every N
%                and is checked against lpiprogram where lpiprogram runs.
%                Deviations, all on inputs lpiprogram accepts wrongly or
%                not at all: the dummy variables' ispvar check tests the
%                dummy variables (lpiprogram tests vartab again); names
%                may be given as a cellstr or a container; no spatial
%                variable gives dom = zeros(0,2), where lpiprogram errors.
% MMP, 10/01/2026: Reduced to an adapter: 'lpiprogram' itself lost its cap,
%                and gained the two other deviations above (the dummy check
%                and zero variables), so the copy of its parsing and
%                construction is deleted (BEGIN/END below) and the new
%                input forms are converted and handed to it. The entry
%                above therefore describes deleted code. Outputs and error
%                messages are unchanged on every form tested
%                (tests/test_lpiprogram_sop and the HEAD comparison).

% BEGIN MMP, 10/01/2026: 'lpiprogram' itself now takes any number of        % MMP, 10/01/2026
% spatial variables, checks the dummy variables and accepts none, so this   % MMP, 10/01/2026
% keeps only what lpiprogram does not read - a container or block, a        % MMP, 10/01/2026
% struct of 'vars' and 'dom', names as a cellstr - and calls it. Deleted:   % MMP, 10/01/2026
% the copy of lpiprogram's parsing and construction (Initial coding         % MMP, 10/01/2026
% 09/29/2026), whose three deviations lpiprogram now shares.                % MMP, 10/01/2026
a = varargin;                                                               % MMP, 10/01/2026
if isempty(a)                                                               % MMP, 10/01/2026
    error("No domain has been specified for the spatial variables.")        % MMP, 10/01/2026
end                                                                         % MMP, 10/01/2026
% A container, a block or a struct carries its own domains: (names,[],dom,...). % MMP, 10/01/2026
if is_sop_class(a{1}) || (isstruct(a{1}) && isfield(a{1},'vars') && isfield(a{1},'dom')) % MMP, 10/01/2026
    [nm,dm] = registry_of(a{1});                                            % MMP, 10/01/2026
    a = [{nm,[],dm}, a(2:end)];                                             % MMP, 10/01/2026
end                                                                         % MMP, 10/01/2026
% Names given as cellstr (the registry 'copvar.vars'), for the variables    % MMP, 10/01/2026
% and for the dummies.                                                      % MMP, 10/01/2026
if iscellstr(a{1}) || isstring(a{1}),   a{1} = names2pvar(a{1});    end %#ok<ISCLSTR> % MMP, 10/01/2026
if numel(a)>=2 && (iscellstr(a{2}) || isstring(a{2})), a{2} = names2pvar(a{2}); end %#ok<ISCLSTR> % MMP, 10/01/2026
prog = lpiprogram(a{:});                                                    % MMP, 10/01/2026
% END MMP, 10/01/2026

end



% ========================================================================
function tf = is_sop_class(x)
tf = isa(x,'copvar') || isa(x,'cdopvar') || isa(x,'sopvar') || isa(x,'sdopvar');
end


% ========================================================================
function [nm,dm] = registry_of(X)
% Names (1 x nv cellstr, sorted) and domains (nv x 2) of the variables of
% X. A container stores both as its registry (vars, dom); a block stores
% them per side (vars.out/in, dom.out/in), which agree on a shared
% variable (sopvar class invariant); a struct gives them directly.
if isa(X,'copvar') || isa(X,'cdopvar')
    nm = X.vars;    dm = X.dom;
elseif isstruct(X)
    nm = cellstr(X.vars);    dm = X.dom;
else
    vo = cellstr(X.vars.out);   vi = cellstr(X.vars.in);
    [nm,io] = unique([vo(:); vi(:)]);
    D = [reshape(X.dom.out,[],2); reshape(X.dom.in,[],2)];
    dm = D(io,:);
end
nm = reshape(cellstr(nm),[],1);
if isempty(nm),     dm = zeros(0,2);    end
end


% ========================================================================
function p = names2pvar(names)
% Column 'polynomial' of the named variables, in the given order.
names = cellstr(names);
if isempty(names),  p = polynomial(zeros(0,1));     return,     end
p = polynomial(names(:));
end
