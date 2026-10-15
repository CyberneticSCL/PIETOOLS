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
%           - if ref(1).type = '()', a 'polyopvar' object with
%             coefficients C(ref(1).subs{:}), for C = obj.C;
%           - if ref(1).type = '{}', a 'polyopvar' object representing
%             the term(s) C.ops{i,j} paired with the corresponding
%             monomial(s) Z(j,:), for i,j specified by ref(1).subs;
%           If ref has more than one element, the remaining elements
%           ref(2:end) are applied to the resulting object;
%
% NOTES
% - For ref(1).type = '()', ref(1).subs index the rows and columns of the
%   'tensopmat' obj.C, see 'subsref' for 'tensopmat'.
% - For ref(1).type = '{}', either two indices or a single scalar linear
%   index must be specified.
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

if isempty(ref)
    out = obj;
    return
end

switch ref(1).type
    case '.'
        % Just use the built-in function
        out = builtin('subsref',obj,ref);
        
    case '()'
        % Extract elements of the matrix-valued polynomial
        out = obj;
        % Use 'tensopvar' subsref to extract desired elements
        C = obj.C;
        C = subsref(C,ref(1));
        out.C = C;
        out.dims.out = C.dims.out;
        out.vars.out = C.vars.out;
        out.dom.out = C.dom.out;
        out.depmat.out = C.depmat.out;
        out.chkval = false;
        out = initialize(out);

        % Perform remaining subsref
        if numel(ref)>=2
            out = subsref(out,ref(2:end));
        end
        
    case '{}'
        % When calling f{i,j}, extract the term associated to the jth
        % monomial in the ith function
        [ncomps,ntrms] = size(obj.C);
        if isscalar(ref(1).subs)
            % Support only a single linear index
            if ~isnumeric(ref(1).subs{1}) || ~isscalar(ref(1).subs{1})
                error("Linear indexing of 'polyopvar' objects is not supported.")
            end
            [indr,indc] = ind2sub([ncomps,ntrms],ref(1).subs{1});
            ref(1).subs{1} = indr;
            ref(1).subs{2} = indc;
        elseif numel(ref(1).subs)==2
            % Extract row and column indices
            indr = ref(1).subs{1};
            indc = ref(1).subs{2};
        else
            error("At most two indices are supported for 'polyopvar' objects.")
        end
        if any(indr>ncomps)
            error("Row indices cannot exceed number of components.")
        elseif any(indc>ntrms)
            error("Column indices cannot exceed number of monomials.")
        end

        % Build the extracted block(s) into their own 'tensopmat', reusing
        % its cell constructor (which calls its own 'initialize') rather
        % than @tensopmat/subsref.m's '{}' (unwraps to a bare factor) or
        % '()' (re-slices element-wise within each block) cases.
        C = tensopmat(obj.C.ops(indr,indc));

        % Retain only the corresponding monomials, and let 'initialize'
        % recompute the output side (and validate the result) from the
        % now-reduced coefficient operator.
        out = obj;
        out.C = C;
        out.Z = obj.Z(indc,:);
        out.dims.out = [];
        out.vars.out = cell(1,0);
        out.dom.out = zeros(0,2);
        out.depmat.out = [];
        out.chkval = false;
        out = initialize(out);

        % Perform remaining subsref
        if numel(ref)>=2
            out = subsref(out,ref(2:end));
        end

end
end
