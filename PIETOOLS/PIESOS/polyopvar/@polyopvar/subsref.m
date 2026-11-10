function out = subsref(obj,ref)
% OUT = SUBSREF(OBJ,REF) extracts the elements of the 'polyopvar' object
% OBJ specified by REF. This function is automatically called when calling
% e.g. f(i,j), f{i,j}, or f.C for a 'polyopvar' object f.
%
% INPUTS
% - obj:    'polyopvar' object of which to extract certain elements;
% - ref:    struct array indicating which elements/fields of obj to
%           extract, with fields 'type' ('.', '()', or '{}') and 'subs';
%
% OUTPUTS
% - out:    object representing the elements of obj specified by ref:
%           - if ref(1).type = '.', out = obj.(ref(1).subs);
%           - if ref(1).type = '()', a 'polyopvar' object representing
%             the element(s) (i,j) of the matrix-valued polynomial obj,
%             with coefficients C(i,:), for C = obj.C;
%           - if ref(1).type = '{}', a 'polyopvar' object representing
%             the block(s) (i,j) of the matrix-valued polynomial obj,
%             i.e. the rows of blocks C{i,:} of the coefficients, paired
%             with the monomials Z(k,:) of its nonempty terms C.ops{i,k};
%           If ref has more than one element, the remaining elements
%           ref(2:end) are applied to the resulting object;
%
% NOTES
% - Only vector-valued polynomials are supported: obj is a single column
%   (see 'size'), so that the column index j of obj(i,j) and obj{i,j} can
%   only be 1 or ':'. With a single index k, obj(k) is obj(k,1) and obj{k}
%   is obj{k,1}.
% - For ref(1).type = '()', i counts the scalar rows of obj, i.e. of the
%   'tensopmat' obj.C, see 'subsref' for 'tensopmat'. All monomials are
%   retained. If all blocks are scalar-valued, obj(i,1) represents the
%   same polynomial as obj{i,1}.
% - For ref(1).type = '{}', i counts the blocks of rows of obj.C.ops
%   (components, of row dimension obj.dims.out); ':' selects all of them.
%   Monomials for which all extracted blocks are empty are dropped from
%   the output (so f{i} returns the nonempty terms of component i), and an
%   error is thrown if no nonempty block is extracted. A single term of f
%   is extracted by f.C{':',k} and f.Z(k,:), the blocks and the monomial k.
% - The keyword 'end' cannot be used: it evaluates to 1 for a 'polyopvar'
%   (MATLAB's size of the scalar object), so f{end} is f{1} and f(end) is
%   f(1). Use f{size(f.C.ops,1)} or f(size(f,1)) instead.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subsref
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
% DJ, 10/15/2026: Initial coding;
% DJ, 10/23/2026: Support ':' in f{i,j}, and drop the monomials without a
%                   nonempty block; several blocks via C{I,J} of tensopmat;
% DJ, 11/06/2026: Output space of f(i,j), f{i,j} follows from C (see set.C);
% DJ, 11/07/2026: Reject adjoints for f(i,j) and f{i,j}, see 'ctranspose';
% DJ, 11/08/2026: Support f(i,1) to extract element (i,1) of f, and f{i,1}
%                   to extract block (i,1) of (vector-valued) elements;

if isempty(ref)
    out = obj;
    return
end

switch ref(1).type
    case '.'
        % Just use the built-in function
        out = builtin('subsref',obj,ref);
        
    case '()'
        if obj.is_adjoint
            error("polyopvar:adjoint_not_supported","'subsref' is not supported for adjoint 'polyopvar' objects (u'), see 'ctranspose'.")
        end
        % Extract elements of the matrix-valued polynomial
        out = obj;
        % Use 'tensopvar' subsref to extract desired elements
        C = obj.C;
        % obj is vector-valued, a single column: obj(k) is obj(k,1), and obj(i,1) is C(i,:), with all
        % monomials, since the columns of C are monomials and not columns of obj
        if isscalar(ref(1).subs)
            ref(1).subs{2} = ':';
        elseif numel(ref(1).subs)==2 && ~strcmp(ref(1).subs{2},':')
            if ~isequal(ref(1).subs{2},1)
                error("The column index must be 1 or ':', only vector-valued 'polyopvar' objects are supported.")
            end
            ref(1).subs{2} = ':';
        end
        C = subsref(C,ref(1));
        % The output space of out is set by the assignment of out.C, see set.C
        out.C = C;
        out.chkval = false;
        out = initialize(out);

        % Perform remaining subsref
        if numel(ref)>=2
            out = subsref(out,ref(2:end));
        end
        
    case '{}'
        if obj.is_adjoint
            error("polyopvar:adjoint_not_supported","'subsref' is not supported for adjoint 'polyopvar' objects (u'), see 'ctranspose'.")
        end
        % When calling f{i,j}, extract block (i,j) of the matrix-valued 
        % polynomial, for j=1 (f is vector-valued).
        [ncomps,ntrms] = size(obj.C.ops);
        if isscalar(ref(1).subs)
            % Support only a single linear index
            if ~isnumeric(ref(1).subs{1}) || ~isscalar(ref(1).subs{1})
                error("Linear indexing of 'polyopvar' objects is not supported.")
            end
            % Linear index along the single column of blocks, f{k} is f{k,1}
            indr = ref(1).subs{1};
            indc = 1;
            ref(1).subs{1} = indr;
            ref(1).subs{2} = indc;
        elseif numel(ref(1).subs)==2
            % Extract row and column indices
            indr = ref(1).subs{1};
            indc = ref(1).subs{2};
        else
            error("At most two indices are supported for 'polyopvar' objects.")
        end
        % ':' selects all components (else compared as a char code below)
        if strcmp(indr,':'),    indr = 1:ncomps;    end
        % Only vector-valued polynomials: a single column of blocks, which contains all monomials
        if strcmp(indc,':'),    indc = 1;           end
        if ~isequal(indc,1)
            error("The column index must be 1 or ':', only vector-valued 'polyopvar' objects are supported.")
        end
        % All monomials, those without nonempty block in the extracted rows are dropped below
        indc = 1:ntrms;
        if any(indr>ncomps)
            error("Row indices cannot exceed number of components.")
        end

        % Build the extracted block(s) into their own 'tensopmat', reusing
        % its cell constructor (which calls its own 'initialize').
        % Drop the monomials without a nonempty block in extracted rows.
        ops = obj.C.ops(indr,indc);
        keep = any(~cellfun('isempty',ops),1);
        if ~any(keep)
            error("Extracted component(s) contain no nonempty term.")
        end
        indc = indc(keep);
        if isscalar(indr) && isscalar(indc)
            C = tensopmat(ops(:,keep));
        else
            % Slices dims/depmat from obj.C, so a row without a block keeps its dimensions
            C = obj.C{indr,indc};
        end

        % Retain only the corresponding monomials; the output side follows from
        % the now-reduced coefficient operator (set.C), and 'initialize'
        % validates the result.
        out = obj;
        out.C = C;
        out.Z = obj.Z(indc,:);
        out.chkval = false;
        out = initialize(out);

        % Perform remaining subsref
        if numel(ref)>=2
            out = subsref(out,ref(2:end));
        end

end
end
