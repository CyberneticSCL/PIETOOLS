function varargout = merge_copvar_registry(cls,opname,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [A,B,...] = MERGE_COPVAR_REGISTRY(CLS,OPNAME,A,B,...) restates containers
% over the union of their variable registries, so that 'plus' and 'mtimes'
% of 'copvar' and 'cdopvar' accept operands built over different ones.
%
% INPUTS
% - cls:    'copvar' or 'cdopvar', the class named in the error identifier
%           cls:domConflict, as for concatenation;
% - opname: the operation named in the error message, e.g. 'plus';
% - A,B,..: 'copvar' or 'cdopvar' containers, in operand order;
%
% OUTPUTS
% - A,B,..: the same operators, each over the sorted union registry, with
%           dom and both space masks remapped by variable NAME. Blocks,
%           dimensions and decision variable lists are untouched;
%
% NOTES
% The registry is container metadata (Sec. 8; see 'copvar'): every block
% carries its own vars and dom, and a mask over a larger registry, false in
% the added columns, names the same space. Restating therefore changes no
% operator. It is needed because an R^n -> R^n operator (-gam*Iw, Dzu, Dyw)
% has an EMPTY registry, and a product such as PB'*Tw keeps its left
% factor's registry although it maps R^n -> R^n, so operands that are
% compatible space by space can still disagree on the registry.
%
% The merge rule is the one concatenation applies in 'cat_copvar_grid':
% sorted union, and a variable given two different domains is an error,
% cls:domConflict, with the same message. 'cat_copvar_grid' keeps its own
% copy of that loop (09/25/2026); it could call this routine instead.
%
% When all registries already agree nothing is built, and an operand that
% is already on the union is returned as it is.
%
% Cost: O(K*nv) for the merge, plus O((M+N)*nv) logical copies per restated
% operand, nv the number of spatial variables (one per direction). No
% coefficient matrix or decision variable list is read or copied, so the
% cost is flat in the number of decision variables. Measured on 2 x 2
% 'poscopvar' containers, q = 1e3 to 3.9e5, 1-D and 2-D, union nv = 1 to 4:
% 3 us when the registries agree, 20-40 us when one is restated, no memory
% held; P*A with A restated onto {s1,...,s4} runs as fast as P*A.
%
% See also CAT_COPVAR_GRID, PLUS, MTIMES, COPVAR, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - merge_copvar_registry
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
% Initial coding MMP, 09/26/2026

K = numel(varargin);
varargout = varargin;

% % % One registry already: the common case, nothing to build.
same = true;
for k = 2:K
    if ~isequal(varargin{k}.vars,varargin{1}.vars) || ...
            ~isequal(varargin{k}.dom,varargin{1}.dom)
        same = false;
        break
    end
end
if same,    return,     end

% % % The union, sorted as the registry must be, and one domain per name.
allv = cell(1,0);
for k = 1:K
    allv = [allv, reshape(varargin{k}.vars,1,[])];                 %#ok<AGROW>
end
vars = reshape(unique(allv),1,[]);
dom = nan(numel(vars),2);
for k = 1:K
    [~,loc] = ismember(varargin{k}.vars,vars);
    for r = 1:numel(loc)
        d = varargin{k}.dom(r,:);
        if all(isnan(dom(loc(r),:)))
            dom(loc(r),:) = d;
        elseif any(dom(loc(r),:)~=d)
            error([cls ':domConflict'],['Variable ''%s'' is on [%g,%g] in '...
                'operand %d of %s and on [%g,%g] in an earlier operand.'],...
                vars{loc(r)},d(1),d(2),k,opname,dom(loc(r),1),dom(loc(r),2))
        end
    end
end

% % % Remap each operand's masks onto the union, column by variable name.
% Assigned in place, so the class, the blocks and Zd carry over as they are.
nv = numel(vars);
for k = 1:K
    X = varargin{k};
    if isequal(X.vars,vars) && isequal(X.dom,dom),  continue,   end
    [~,loc] = ismember(X.vars,vars);
    so = false(size(X.space_out,1),nv);     so(:,loc) = X.space_out;
    si = false(size(X.space_in,1),nv);      si(:,loc) = X.space_in;
    X.vars = vars;          X.dom = dom;
    X.space_out = so;       X.space_in = si;
    varargout{k} = X;
end

end
