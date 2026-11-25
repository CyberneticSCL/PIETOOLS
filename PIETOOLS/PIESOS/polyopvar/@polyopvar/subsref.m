function out = subsref(obj,ref)
% OUT = SUBSREF(OBJ,REF) extracts the elements or fields of the
% 'polyopvar' object OBJ specified by REF. Called for e.g. f(i,j), f{i,j}
% or f.C.
%
% INPUTS
% - obj:    'polyopvar' object, not an adjoint for '()' and '{}';
% - ref:    'struct' array with fields 'type' ('.', '()' or '{}') and
%           'subs', the elements or fields of obj to extract;
%
% OUTPUTS
% - out:    - for '.', the field obj.(ref(1).subs);
%           - for '()', a 'polyopvar' object of the rows i of obj, with
%             coefficients C(i,:) for C = obj.C, and all monomials;
%           - for '{}', a 'polyopvar' object of the row-blocks i of obj,
%             C{i,:}, with only the monomials on which a nonempty block
%             of these acts (error if there is none);
%           ref(2:end) is then applied to the result;
%
% NOTES
% - obj is a single column, so the column index j of obj(i,j) and
%   obj{i,j} can only be 1 or ':'; obj(k) is obj(k,1), obj{k} is obj{k,1}.
% - For '()', i counts the scalar rows of obj; for '{}', the row-blocks
%   of obj.C (of dimension obj.dims.out). If all row-blocks are scalar,
%   obj(i) and obj{i} represent the same polynomial.
% - A single term of f is f.C{':',k} acting on the monomial f.Z(k,:).
% - The keyword 'end' evaluates to 1 for a 'polyopvar' object: use
%   f(size(f,1)) or f{size(f.C.ops,1)} instead.
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 10/15/2026: Initial coding;

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
        % The output space of out is set by the assignment of out.C
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
