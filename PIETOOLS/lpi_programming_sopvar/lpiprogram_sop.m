function prog = lpiprogram_sop(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPIPROGRAM_SOP(VARTAB,DUMVARTAB,DOM,DECVARTAB,FREEVARTAB) declares
% an LPI program in the n spatial variables VARTAB on the domains DOM, with
% the dummy variables DUMVARTAB, the decision variables DECVARTAB and the
% independent variables FREEVARTAB, which have no domain: it converts the
% names given as a 'cellstr' or string array to 'polynomial' variables and
% builds the program by LPIPROGRAM. The variables and domains may instead
% be read off an operator P or a struct S:
%
%   PROG = LPIPROGRAM_SOP(VARTAB,DUMVARTAB,DOM,DECVARTAB,FREEVARTAB)
%   PROG = LPIPROGRAM_SOP(VARTAB,DOM,DECVARTAB,FREEVARTAB)
%   PROG = LPIPROGRAM_SOP(P,DECVARTAB,FREEVARTAB)
%   PROG = LPIPROGRAM_SOP(S,DECVARTAB,FREEVARTAB)
%
% with DECVARTAB and FREEVARTAB optional and in either order. Without
% DUMVARTAB, the dummy variable of a spatial variable s is s_dum.
%
% INPUT
% - vartab:     n x 1 'polynomial' of the spatial variables, or a 'cellstr'
%               or string array of their names; in the form without
%               DUMVARTAB also an n x 2 'polynomial' whose second column
%               holds the dummy variables, n = size(DOM,1);
% - dumvartab:  (optional) n x 1 'polynomial', 'cellstr' or string array of
%               the dummy variables, entry i that of VARTAB(i); [] gives
%               the default names;
% - dom:        n x 2 'double', row i the interval [a,b], a < b, of
%               VARTAB(i); one 1 x 2 row applies to every variable;
% - decvartab:  (optional) q x 1 'dpvar' of decision variables;
% - freevartab: (optional) m x 1 'polynomial' of independent variables
%               without a domain;
% - P:          'copvar', 'cdopvar', 'sopvar' or 'sdopvar' object, whose
%               spatial variables, sorted by name, and their domains are
%               taken;
% - S:          'struct' with fields vars ('cellstr' or string array of n
%               names, taken in the order given) and dom (n x 2 or 1 x 2
%               'double', the domains).
%
% OUTPUT
% - prog:       'struct', an SOS program structure (SOSPROGRAM) with the
%               fields
%                 vartable     (2n+m) x 1 'polynomial', [VARTAB; DUMVARTAB;
%                              FREEVARTAB];
%                 dom          n x 2 'double', the domains;
%                 decvartable  'cellstr', the names of DECVARTAB, sorted ({}
%                              without DECVARTAB);
%               and the other fields of a new SOS program (var, expr,
%               extravar, objective, solinfo, varmat), with no constraint.
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
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body and headers
%                of the local functions; code unchanged. The deleted
%                header text was history, the comparison with lpiprogram,
%                the routines that read prog.vartable and prog.dom, and a
%                cost note.

% BEGIN MMP, 10/01/2026: 'lpiprogram' itself now takes any number of        % MMP, 10/01/2026
% spatial variables, checks the dummy variables and accepts none, so this   % MMP, 10/01/2026
% keeps only what lpiprogram does not read - a container or block, a        % MMP, 10/01/2026
% struct of 'vars' and 'dom', names as a cellstr - and calls it. Deleted:   % MMP, 10/01/2026
% the copy of lpiprogram's parsing and construction (Initial coding         % MMP, 10/01/2026
% 09/29/2026), whose three deviations lpiprogram now shares.                % MMP, 10/01/2026
% The arguments; with none, no domain is given
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
% The program, built by lpiprogram from the converted arguments
prog = lpiprogram(a{:});                                                    % MMP, 10/01/2026
% END MMP, 10/01/2026

end



% ========================================================================
function tf = is_sop_class(x)
% TF = IS_SOP_CLASS(X) returns true when X is a 'copvar', 'cdopvar',
% 'sopvar' or 'sdopvar' object.
% INPUT  x:  any value.
% OUTPUT tf: logical.
tf = isa(x,'copvar') || isa(x,'cdopvar') || isa(x,'sopvar') || isa(x,'sdopvar');
end


% ========================================================================
function [nm,dm] = registry_of(X)
% [NM,DM] = REGISTRY_OF(X) returns the names and the domains of the spatial
% variables of X.
% INPUT  X:  'copvar' or 'cdopvar' object; 'sopvar' or 'sdopvar' object;
%            or 'struct' with fields vars ('cellstr' or string array) and
%            dom ('double').
% OUTPUT nm: nv x 1 cellstr, sorted for an object, in the given order for
%            a struct;
%        dm: nv x 2 double, row k the domain of nm{k} (a struct's dom as
%            given); zeros(0,2) when nv = 0.
if isa(X,'copvar') || isa(X,'cdopvar')
    % A container stores both as its registry (vars, dom)
    nm = X.vars;    dm = X.dom;
elseif isstruct(X)
    % A struct gives them directly
    nm = cellstr(X.vars);    dm = X.dom;
else
    % A block stores them per side (vars.out/in, dom.out/in), which agree
    % on a shared variable (sopvar class invariant): the sorted union
    vo = cellstr(X.vars.out);   vi = cellstr(X.vars.in);
    [nm,io] = unique([vo(:); vi(:)]);
    D = [reshape(X.dom.out,[],2); reshape(X.dom.in,[],2)];
    dm = D(io,:);
end
% Names as a column; a 0 x 2 domain for no variable
nm = reshape(cellstr(nm),[],1);
if isempty(nm),     dm = zeros(0,2);    end
end


% ========================================================================
function p = names2pvar(names)
% P = NAMES2PVAR(NAMES) returns the column 'polynomial' of the variables
% named NAMES, in the given order.
% INPUT  names: 'cellstr' or string array.
% OUTPUT p:     numel(NAMES) x 1 'polynomial'; 0 x 1 for no name.
names = cellstr(names);
if isempty(names),  p = polynomial(zeros(0,1));     return,     end
p = polynomial(names(:));
end
