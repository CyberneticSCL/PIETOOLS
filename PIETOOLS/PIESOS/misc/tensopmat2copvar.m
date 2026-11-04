function Pc = tensopmat2copvar(Pt,dims,vars,depmat)
% PC = TENSOPMAT2COPVAR(PT,DIMS,VARS,DEPMAT) converts the linear operator
% Pt, a 'tensopmat' object, to a 'copvar' object, with a block for each
% pair of a space in the output and input space, composed of the blocks of
% Pt between the components of these spaces.
%
% INPUT
% - Pt:     'tensopmat' object representing a linear operator, with each
%           block either empty, a 'tensopvar' object involving a single
%           'sopvar' object, a 'sopvar' object, or a constant matrix of
%           type 'double' or 'polynomial' mapping from R;
% - dims:   (optional) struct with fields
%       out:    nspaces_out x 1 array with the sizes of the spaces in the
%               output space of the copvar;
%       in:     nspaces_in x 1 array with the sizes of the spaces in the
%               input space of the copvar;
% - vars:   (optional) struct with fields
%       out:    1 x N cellstr with the names of the variables in the
%               columns of depmat.out;
%       in:     1 x N cellstr with the names of the variables in the
%               columns of depmat.in;
%           Default to Pt.vars.out and Pt.vars.in;
% - depmat: (optional) struct with fields
%       out:    nspaces_out x N logical array, with element (i,j) true if
%               space i of the output space depends on variable
%               vars.out{j};
%       in:     nspaces_in x N logical array, similarly for the input
%               space;
%
% OUTPUT
% - Pc:     'copvar' object representing the same operator as Pt, with a
%           row for each space in the output space, and a column for each
%           space in the input space. An empty 'copvar' if there are no
%           components in the output or the input space of Pt.
%
% NOTES
% - If dims.out (dims.in) is omitted or [], each component (row or column
%   of blocks) of Pt is a separate space of Pc, with the same size and
%   dependence on the variables, so that e.g. a 'tensopmat' mapping to
%   L2^2[s] x L2^3[s] gives a 'copvar' mapping to the same space.
%   Otherwise, vars.out and depmat.out (vars.in and depmat.in) describe the
%   spaces of Pc, to which the components of Pt are assigned in order, so
%   that e.g. dims.out=5 and depmat.out=1 combine these components in
%   L2^5[s].
% - Each component of Pt must be a part of a specified space with the same
%   variables, and the components must fill the spaces.
% - The blocks of Pt cannot involve tensor products of multiple operators
%   or functionals ('intop' objects), nor operators depending on variables
%   to a degree higher than 1, as these are not linear PI operators.
% - See also COPVAR2TENSOPMAT.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - tensopmat2copvar
%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 11/03/2026: Initial coding;
% DJ (CC), 11/03/2026: Reuse metadata of zero operator for the output;

if ~isa(Pt,'tensopmat')
    error("Operator must be of type 'tensopmat'.")
end
if nargin<2
    dims = [];
end
if nargin<3
    vars = [];
end
if nargin<4
    depmat = [];
end
% A field is only used if it is specified, and not equal to [].
isgiven = @(S,sd) isstruct(S) && isfield(S,sd) && ~isequal(S.(sd),[]);
if any(Pt.depmat.in(:)>1) || any(Pt.depmat.out(:)>1)
    error("Only linear operators are supported, which depend on each variable to degree at most 1.")
end

% Declare the names and domain of all variables in the copvar registry.
allnames = [Pt.vars.out,Pt.vars.in];
alldom = [Pt.dom.out;Pt.dom.in];
[vnames,idx] = unique(allnames);
dom = alldom(idx,:);
for k=1:numel(vnames)
    is_k = strcmp(allnames,vnames{k});
    if any(any(alldom(is_k,:)~=dom(k,:)))
        error(['The variable ',vnames{k},' is defined on different domains in the input and output.'])
    end
end
meta = struct('vars',{vnames},'dom',dom);

% Determine the components of Pt and the spaces of Pc, and assign each of
% the components to a space and a range of indices in this space.
sides = {'out','in'};
cdims = cell(1,2);      cnames = cell(1,2);     
gdims = cell(1,2);      gnames = cell(1,2);
crow = cell(1,2);       cstart = cell(1,2);
for j=1:2
    sd = sides{j};
    cdims{j} = Pt.dims.(sd)(:);
    cnames{j} = arrayfun(@(k) reshape(sort(Pt.vars.(sd)(Pt.depmat.(sd)(k,:)>0)),1,[]),1:numel(cdims{j}),'UniformOutput',false);
    if ~isgiven(dims,sd)
        % Each component is a separate space.
        gdims{j} = cdims{j};
        gnames{j} = cnames{j};
    else
        gdims{j} = dims.(sd)(:);
        Pvars = Pt.vars.(sd);
        if isgiven(vars,sd)
            Pvars = vars.(sd);
        end
        if ~isgiven(depmat,sd) || size(depmat.(sd),1)~=numel(gdims{j}) || size(depmat.(sd),2)~=numel(Pvars)
            error(['depmat.',sd,' should be a logical array with a row for each space in dims.',sd,' and a column for each variable in vars.',sd,'.'])
        end
        gnames{j} = arrayfun(@(k) reshape(sort(Pvars(logical(depmat.(sd)(k,:)))),1,[]),1:numel(gdims{j}),'UniformOutput',false);
    end
    used = zeros(numel(gdims{j}),1);
    crow{j} = zeros(numel(cdims{j}),1);     cstart{j} = zeros(numel(cdims{j}),1);
    for k=1:numel(cdims{j})
        r = find(cellfun(@(n) isequal(n,cnames{j}{k}),gnames{j}) & gdims{j}(:)'-used(:)'>=cdims{j}(k),1);
        if isempty(r)
            error("The components of the "+sd+"put space do not match the specified spaces.")
        end
        crow{j}(k) = r;
        cstart{j}(k) = used(r);
        used(r) = used(r) + cdims{j}(k);
    end
    if any(used~=gdims{j}(:))
        error("The components of the "+sd+"put space do not fill the specified spaces.")
    end
end
if isempty(gdims{1}) || isempty(gdims{2})
    Pc = copvar();
    return
end

% Declare the block of each pair of spaces, from the blocks of the
% components in these spaces, starting from the zero operator.
Zop = zeros_copvar_sop(struct('out',gdims{1},'in',gdims{2}),...
        struct('out',{gnames{1}},'in',{gnames{2}}),meta);
C = Zop.C;
for r=1:numel(gdims{1})
    for c=1:numel(gdims{2})
        % Components in the spaces, in the order in which they fill them.
        [~,order_r] = sort(cstart{1}(crow{1}==r));      comps_r = find(crow{1}==r);     comps_r = comps_r(order_r);
        [~,order_c] = sort(cstart{2}(crow{2}==c));      comps_c = find(crow{2}==c);     comps_c = comps_c(order_c);
        Bcell = cell(numel(comps_r),numel(comps_c));
        is_zero = true;
        for a=1:numel(comps_r)
            for b=1:numel(comps_c)
                i = comps_r(a);     j = comps_c(b);
                spaces_ij = struct('out',{{cnames{1}{i}}},'in',{{cnames{2}{j}}});
                Bij = Pt.ops{i,j};
                if isempty(Bij)
                    Zij = zeros_copvar_sop(struct('out',cdims{1}(i),'in',cdims{2}(j)),spaces_ij,meta);
                    Bcell{a,b} = Zij.C{1,1};
                    continue
                end
                is_zero = false;
                if isa(Bij,'tensopvar')
                    if ~isequal(size(Bij.ops),[1,1]) || ~isa(Bij.ops{1,1},'sopvar')
                        error("Only blocks involving a single 'sopvar' object are supported.")
                    end
                    Bcell{a,b} = Bij.ops{1,1};
                elseif isa(Bij,'sopvar')
                    Bcell{a,b} = Bij;
                elseif isa(Bij,'intop') || isa(Bij,'dpvar')
                    error("Blocks of type '"+class(Bij)+"' are not supported.")
                elseif isnumeric(Bij) || isa(Bij,'polynomial')
                    % Matrix acting on R, mapping to a constant function.
                    Mij = polynomial(Bij);
                    if ~isempty(Mij.varname)
                        error("Matrix blocks must be constant, they cannot depend on spatial variables.")
                    elseif ~isempty(cnames{2}{j})
                        error("Matrix blocks must map from R.")
                    end
                    Mij = mat2copvar_sop(double(Mij),struct('out',cdims{1}(i),'in',cdims{2}(j)),spaces_ij,meta);
                    Bcell{a,b} = Mij.C{1,1};
                else
                    error("Blocks of type '"+class(Bij)+"' are not supported.")
                end
                % The operator must map between the component spaces.
                Bdims = Bcell{a,b}.dims;
                if ~isequal(Bdims(:)',[cdims{1}(i),cdims{2}(j)]) || ~isequal(reshape(sort(Bcell{a,b}.vars.out),1,[]),cnames{1}{i}) ...
                        || ~isequal(reshape(sort(Bcell{a,b}.vars.in),1,[]),cnames{2}{j})
                    error("The block of the operator does not map between the spaces of its components.")
                end
            end
        end
        if is_zero
            continue
        end
        % Concatenate the blocks of the components into a single block.
        Rows = arrayfun(@(a) [Bcell{a,:}],1:size(Bcell,1),'UniformOutput',false);
        C{r,c} = vertcat(Rows{:});
    end
end
% Every block maps between the spaces of the zero operator (checked above),
% so reuse its metadata to reduce computation time.
Pc = copvar(C,metadata(Zop));

end
