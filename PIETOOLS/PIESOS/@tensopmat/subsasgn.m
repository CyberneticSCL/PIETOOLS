function out = subsasgn(obj,ref,RHS)
% OUT = SUBSASGN(OBJ,REF,RHS) assigns the value RHS to the field or the
% block(s) REF of the 'tensopmat' object OBJ. This function is
% automatically called when calling e.g. C{i,j} = RHS or C.ops = RHS for a
% 'tensopmat' object C.
%
% INPUTS
% - obj:    'tensopmat' object of which to assign elements;
% - ref:    'struct' specifying which elements of obj to adjust;
% - RHS:    object specifying the new value(s) of the elements specified
%           by ref. For ref(1).type='{}', this is one of:
%             - if ref(1).subs is a single (scalar) row,column pair: the
%               'tensopvar' or 'intop' object, or 'double', 'polynomial'
%               or 'dpvar' matrix, to store directly as obj.ops{i,j};
%             - if ref(1).subs is a pair of index vectors: a 'tensopmat'
%               object whose .ops matches the selected shape, substituted
%               into obj.ops(idcs1,idcs2) as a block;
%
% OUTPUTS
% - out:    'tensopmat' object with the desired component value;
%
% NOTES
% - Only '{}'-type subsasgn (replacing row/column-block(s) in obj.ops) and
%   '.'-type subsasgn (setting a field) are supported; '()'-type subsasgn
%   is not supported.
% - After '{}'-type subsasgn, the object is passed through 'initialize',
%   which recomputes dims/vars/dom/depmat from obj.ops and errors if the
%   new block(s) are inconsistent in matrix dimension, variables, or
%   domain with the rest of the row(s)/column(s) they belong to.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subsasgn
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 10/15/2026: Initial coding;

switch ref(1).type
    case '.'
        out = builtin('subsasgn',obj,ref,RHS);
    case '()'
        error("'()'-type subsassign is not supported for tensopmat objects.")
    case '{}'
        % Use '{}'-type subsindexing to set values of obj.ops
        [ncomps,ntrms] = size(obj.ops);
        if isscalar(ref(1).subs)
            % Support only a single linear index
            if ~isnumeric(ref(1).subs{1}) || ~isscalar(ref(1).subs{1})
                error("Linear indexing of 'tensopmat' objects is not supported.")
            end
            [indr,indc] = ind2sub([ncomps,ntrms],ref(1).subs{1});
        elseif numel(ref(1).subs)==2
            % Extract row and column numbers
            indr = ref(1).subs{1};
            indc = ref(1).subs{2};
        else
            error("At most two indices are supported for 'tensopmat' objects.")
        end
        if any(indr<=0) || any(indc<=0)
            error("Indices must be positive integers.")
        elseif any(indr>ncomps)
            error("Row indices cannot exceed number of row-blocks.")
        elseif any(indc>ntrms)
            error("Column indices cannot exceed number of column-blocks.")
        end
        if numel(ref)>=2
            % e.g. C{i,j}.field=... -- would require reaching into the
            % extracted block, reassigning it, and only then validating.
            error("Chained assignment into a 'tensopmat' block (e.g. C{i,j}.field=...) is not supported.")
        end

        ref(1).subs = {indr,indc};
        if isscalar(indr) && isscalar(indc)
            % Setting a single block C{i,j}: store the operator directly
            % as obj.ops{i,j}.
            if ~isa(RHS,'tensopvar') && ~isa(RHS,'intop') && ~is_native_op(RHS)
                error("For setting a single block of a 'tensopmat' object, the new value must be specified as 'tensopvar' or 'intop' object, or as matrix of type 'double', 'polynomial' or 'dpvar'.")
            end
            ref(1).type = '()';
            new_ops = builtin('subsasgn',obj.ops,ref(1),{RHS});
        else
            % Setting multiple blocks C{idcs1,idcs2}: new value must be a
            % 'tensopmat' whose .ops matches the selected shape.
            if ~isa(RHS,'tensopmat')
                error("For setting multiple blocks of a 'tensopmat' object, the new value must be specified as 'tensopmat' object.")
            elseif ~isequal(size(RHS.ops),[numel(indr),numel(indc)])
                error("Size of the 'tensopmat' object on the right-hand side must match the number of indices on the left-hand side.")
            end
            ref(1).type = '()';
            new_ops = builtin('subsasgn',obj.ops,ref(1),RHS.ops);
        end

        % Recompute dims/vars/dom/depmat from the updated obj.ops, erroring
        % if the new block(s) are inconsistent (dims, vars, or domains)
        % with the rest of their row/column.
        out = obj;
        out.ops = new_ops;
        out = initialize(out);
end

end
