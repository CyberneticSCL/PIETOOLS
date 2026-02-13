function out = subsref(obj,ref)
% OUT = SUBSREF(OBJ,REF) extracts the elements of 'polyopvar' object OBJ
% specified by REF. This function is automatically called when calling e.g.
% f(i,j) or f.C for a 'polyopvar' object f.
% 
% INPUT
% - obj:    'polyopvar' object of which to extract certain elements;
% - ref:    struct indicating which element/fields of obj to extract;
% 
% OUTPUT
% - out:    object representing the elements of 'obj' specified by 'ref'. 
%           If ref.type='.', this will be out=obj.(ref.subs). If 
%           ref.type = '()', the output will be a polyopvar object 
%           including the elements in the indices specified by ref.subs;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIEOOLS - subsref
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
% DJ, 02/13/2026: Initial coding

if isempty(ref)
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

        % Perform remaining subsref
        if numel(ref)>=2
            out = subsref(out,ref(2:end));
        end
        
    case '{}'
        % When calling f{i,j}, extract the term associated to the jth
        % monomial in the ith function
        [ncomps,ntrms] = size(obj.C.ops);
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
            error("At most two indices are supported for 'tensopvar' objects.")
        end
        if any(indr>ncomps)
            error("Row indices cannot exceed number of components.")
        elseif any(indc>ntrms)
            error("Column indices cannot exceed number of monomials.")
        end

        % Extract the desired coefficients
        C = subsref(obj.C,ref(1));

        % Retain only the corresponding monomials
        out = obj;
        out.C = C;
        out.degmat = obj.degmat(indc,:);

        % Perform remaining subsref
        if numel(ref)>=2
            out = subsref(out,ref(2:end));
        end
        
end
end
