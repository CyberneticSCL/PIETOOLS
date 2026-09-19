function meta = derive_mopvar_meta(C,cls,okclasses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% META = DERIVE_MOPVAR_META(C,CLS,OKCLASSES) derives and validates the
% container metadata of an M x N cell of PI operator blocks, for the
% 'mopvar' and 'mdopvar' constructors.
%
% INPUTS
% - C:          M x N cell of blocks; a cell holding [] is a structurally
%               zero block;
% - cls:        name of the calling container class, for error identifiers
%               and messages ('mopvar' or 'mdopvar');
% - okclasses:  cellstr of block classes this container accepts, e.g.
%               {'sopvar'} or {'sopvar','sdopvar'};
%
% OUTPUTS
% - meta:       struct with fields vars, dom, space_out, space_in, dim_out,
%               dim_in. The decision variable list is NOT handled here; a
%               container that has one reconciles it itself.
%
% NOTES
% 'mopvar' and 'mdopvar' duplicate their METHOD files, following how
% @sopvar and @sdopvar are structured, but they share this routine rather
% than duplicating it: it is the whole of the validation and the part that
% would silently diverge, and the only thing that differs between the two
% containers is which block classes are admitted, which is an argument.
%
% Cost: O(M*N) block visits. The spatial set operations run on the variable
% registry, which has one entry per spatial direction, so they are not on a
% decision-variable-scaled axis.
%
% See also MOPVAR, MDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - derive_mopvar_meta
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
% Initial coding MMP, 09/17/2026

[M,N] = size(C);
occ = ~cellfun(@isempty,C);       % which cells hold an operator

% % % An all-empty row leaves s^i and q_i undetermined: nothing in the grid
% % % says what they are. Put an explicitly zero-valued block there instead.
if M>0 && N>0
    if any(~any(occ,2))
        error([cls ':emptyRow'],['Row(s) %s have no populated block, so the '...
            'output space and dimension are undetermined.'],mat2str(find(~any(occ,2))'))
    end
    if any(~any(occ,1))
        error([cls ':emptyColumn'],['Column(s) %s have no populated block, so '...
            'the input space and dimension are undetermined.'],mat2str(find(~any(occ,1))))
    end
end

% % % Block classes and matrix dimensions. 'dims' must be 1x2: 'sopvar' also
% % % accepts a block-dims matrix (@sopvar/size sums it, @sdopvar/size
% % % returns it raw), which would duplicate this container's role for blocks
% % % that share a variable set. Reject it rather than support both.
for i = 1:M
    for j = 1:N
        if ~occ(i,j),   continue,   end
        if ~any(cellfun(@(k) isa(C{i,j},k),okclasses))
            error([cls ':badBlockClass'],...
                'Block (%d,%d) is a ''%s''; %s blocks must be %s.',...
                i,j,class(C{i,j}),cls,strjoin(okclasses,' or '))
        end
        if numel(C{i,j}.dims)~=2
            error([cls ':blockDimsNotScalar'],['Block (%d,%d) has %d matrix '...
                'dimensions; a %s block needs a 1x2 ''dims''. Use the '...
                'container grid for block structure, not sopvar''s dims matrix.'],...
                i,j,numel(C{i,j}.dims),cls)
        end
    end
end

% % % Variable registry and domains. This runs over spatial directions, of
% % % which there are a handful, so the set operations are not on a hot axis.
names = cell(1,0);  lo = zeros(1,0);  hi = zeros(1,0);
for i = 1:M
    for j = 1:N
        if ~occ(i,j),   continue,   end
        [names,lo,hi] = register(names,lo,hi,C{i,j}.vars.out,C{i,j}.dom.out,i,j,cls);
        [names,lo,hi] = register(names,lo,hi,C{i,j}.vars.in ,C{i,j}.dom.in ,i,j,cls);
    end
end
[vars,ord] = sort(names);
dom = [lo(ord)',hi(ord)'];

% % % Row consistency: every populated block in row i maps into the same
% % % space and with the same output dimension.
space_out = false(M,numel(vars));       dim_out = zeros(M,1);
for i = 1:M
    j0 = find(occ(i,:),1);
    dim_out(i) = C{i,j0}.dims(1);
    space_out(i,:) = ismember(vars,C{i,j0}.vars.out);
    s_i = sort(C{i,j0}.vars.out(:))';
    for j = find(occ(i,:))
        if C{i,j}.dims(1)~=dim_out(i)
            error([cls ':rowDimMismatch'],['Blocks (%d,%d) and (%d,%d) are in '...
                'one row but have output dimensions %d and %d.'],...
                i,j0,i,j,dim_out(i),C{i,j}.dims(1))
        end
        % SET equality: blocks order vars.out as [S2,S3] with S2 = s^i/s^j,
        % so the ORDER legitimately differs along a row.
        if ~isequal(sort(C{i,j}.vars.out(:))',s_i)
            error([cls ':rowSpaceMismatch'],['Blocks (%d,%d) and (%d,%d) are in '...
                'one row but map into {%s} and {%s}.'],...
                i,j0,i,j,strjoin(s_i,','),strjoin(sort(C{i,j}.vars.out(:))',','))
        end
    end
end

% % % Column consistency: the same for the input side.
space_in = false(N,numel(vars));        dim_in = zeros(N,1);
for j = 1:N
    i0 = find(occ(:,j),1);
    dim_in(j) = C{i0,j}.dims(2);
    space_in(j,:) = ismember(vars,C{i0,j}.vars.in);
    s_j = sort(C{i0,j}.vars.in(:))';
    for i = find(occ(:,j))'
        if C{i,j}.dims(2)~=dim_in(j)
            error([cls ':colDimMismatch'],['Blocks (%d,%d) and (%d,%d) are in '...
                'one column but have input dimensions %d and %d.'],...
                i0,j,i,j,dim_in(j),C{i,j}.dims(2))
        end
        if ~isequal(sort(C{i,j}.vars.in(:))',s_j)
            error([cls ':colSpaceMismatch'],['Blocks (%d,%d) and (%d,%d) are in '...
                'one column but map out of {%s} and {%s}.'],...
                i0,j,i,j,strjoin(s_j,','),strjoin(sort(C{i,j}.vars.in(:))',','))
        end
    end
end

meta = struct('vars',{vars},'dom',dom,'space_out',space_out,...
    'space_in',space_in,'dim_out',dim_out,'dim_in',dim_in);

end


% ========================================================================
function [names,lo,hi] = register(names,lo,hi,v,d,i,j,cls)
% Add the variables v with domains d to the registry, rejecting a variable
% that two blocks, or the two sides of one block, place on different
% domains. This is the check that cannot be made one block at a time.

v = v(:)';
if size(d,1)~=numel(v)
    error([cls ':domSizeMismatch'],...
        'Block (%d,%d) lists %d variables against %d domain rows.',i,j,numel(v),size(d,1))
end
for k = 1:numel(v)
    idx = find(strcmp(names,v{k}),1);
    if isempty(idx)
        names{end+1} = v{k};    lo(end+1) = d(k,1);     hi(end+1) = d(k,2); %#ok<AGROW>
    elseif lo(idx)~=d(k,1) || hi(idx)~=d(k,2)
        error([cls ':domConflict'],['Variable ''%s'' is on [%g,%g] in block '...
            '(%d,%d) and [%g,%g] elsewhere in the container.'],...
            v{k},d(k,1),d(k,2),i,j,lo(idx),hi(idx))
    end
end
end
