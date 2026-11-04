function Pt = copvar2tensopmat(P,dims,vars,depmat)
% PT = COPVAR2TENSOPMAT(P,DIMS,VARS,DEPMAT) converts the 'copvar'
% operator P to a 'tensopmat' object, with a block for each pair of a
% component of the output space and a component of the input space, taken
% from the 'sopvar' blocks of P (a zero operator for an empty block).
%
% INPUT
% - P:      'copvar' object, with grid of spaces
%               L2^{p_1}[s^1] x ... x L2^{p_N}[s^N]
%                   --> L2^{q_1}[t^1] x ... x L2^{q_M}[t^M];
% - dims:   (optional) struct with fields
%       out:    ncomps_out x 1 array with the sizes of the components of
%               the output space;
%       in:     ncomps_in x 1 array with the sizes of the components of the
%               input space;
% - vars:   (optional) struct with fields
%       out:    1 x N cellstr with the names of the variables in the
%               columns of depmat.out;
%       in:     1 x N cellstr with the names of the variables in the
%               columns of depmat.in;
%           Default to P.vars;
% - depmat: (optional) struct with fields
%       out:    ncomps_out x N logical array, with element (i,j) true if
%               component i of the output space depends on variable
%               vars.out{j};
%       in:     ncomps_in x N logical array, similarly for the input space;
%
% OUTPUT
% - Pt:     'tensopmat' object representing the same operator as P, with a
%           row for each component of the output space, and a column for
%           each component of the input space. [] if there are no
%           components in the output or the input space.
%
% NOTES
% - If dims.out (dims.in) is omitted or [], each space (row of the grid) of
%   the output (input) space of P with nonzero size is taken as a single
%   component, and vars.out and depmat.out (vars.in and depmat.in) are
%   ignored.
% - A component of the specified output or input space must be a part of a
%   space of P, with the same variables, and the components in the same
%   space must be ordered as in P, so that e.g. two scalar components in
%   L2[s] correspond to the first and second row of a block in the
%   space L2^2[s] of P.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - copvar2tensopmat
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

if ~isa(P,'copvar')
    error("Operator must be of type 'copvar'.")
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

% Determine the components of the output and input spaces, and the space
% (row of the grid of P) and indices of each of the components.
sides = {'out','in'};
cdims = cell(1,2);      cnames = cell(1,2);     crow = cell(1,2);   cidx = cell(1,2);
for s=1:2
    sd = sides{s};
    dims_g = P.(['dim_',sd])(:);
    mask_g = P.(['space_',sd]);
    names_g = arrayfun(@(r) reshape(sort(P.vars(mask_g(r,:))),1,[]),1:numel(dims_g),'UniformOutput',false);
    if ~isgiven(dims,sd)
        % Each (nonzero) space of P is a component.
        keep = dims_g>0;
        cdims{s} = dims_g(keep);        cnames{s} = names_g(keep);
    else
        cdims{s} = dims.(sd)(:);
        Pvars = P.vars;
        if isgiven(vars,sd)    
            Pvars = vars.(sd);
        end
        if ~isgiven(depmat,sd) || size(depmat.(sd),1)~=numel(cdims{s}) || size(depmat.(sd),2)~=numel(Pvars)
            error(['depmat.',sd,' should be a logical array with a row for each component in dims.',sd,' and a column for each variable in vars.',sd,'.'])
        end
        cnames{s} = arrayfun(@(k) reshape(sort(Pvars(logical(depmat.(sd)(k,:)))),1,[]),1:numel(cdims{s}),'UniformOutput',false);
    end
    % Match each component with a space of P, and indices in that space.
    used = zeros(numel(dims_g),1);
    crow{s} = zeros(numel(cdims{s}),1);     cidx{s} = cell(numel(cdims{s}),1);
    for k=1:numel(cdims{s})
        r = find(cellfun(@(n) isequal(n,cnames{s}{k}),names_g) & dims_g(:)'-used(:)'>=cdims{s}(k),1);
        if isempty(r)
            error("The components of the "+sd+"put space do not match the spaces of the operator.")
        end
        crow{s}(k) = r;     
        cidx{s}{k} = used(r) + (1:cdims{s}(k));     
        used(r) = used(r) + cdims{s}(k);
    end
end
if isempty(cdims{1}) || isempty(cdims{2})
    Pt = [];
    return
end

% Declare the block of each pair of components.
Cell = cell(numel(cdims{1}),numel(cdims{2}));
for i=1:size(Cell,1)
    for j=1:size(Cell,2)
        Pij = P.C{crow{1}(i),crow{2}(j)};
        if isempty(Pij)
            Zij = zeros_copvar_sop(struct('out',cdims{1}(i),'in',cdims{2}(j)),...
                    struct('out',{{cnames{1}{i}}},'in',{{cnames{2}{j}}}),struct('vars',{P.vars},'dom',P.dom));
            Pij = Zij.C{1,1};
        else
            Pij = Pij(cidx{1}{i},cidx{2}{j});
        end
        Cell{i,j} = tensopmat(tensopvar(Pij));
    end
end
Rows = arrayfun(@(i) [Cell{i,:}],1:size(Cell,1),'UniformOutput',false);
Pt = vertcat(Rows{:});

end
