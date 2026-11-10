function out = subsref(obj,ref)
% OUT = SUBSREF(OBJ,REF) extracts the elements of 'tensopvar' object OBJ
% specified by REF. This function is automatically called when calling e.g.
% C(i,j), C{i,j} or C.ops for a 'tensopvar' object C.
%
% INPUTS
% - obj:    'tensopvar' object of which to extract certain elements;
% - ref:    struct indicating which element/fields of obj to extract;
%
% OUTPUTS
% - out:    object representing the elements of 'obj' specified by 'ref'.
%           If ref(1).type='.', this will be out=obj.(ref(1).subs). If
%           ref(1).type='()', the output will be a 'tensopvar' object
%           including the elements in the indices specified by
%           ref(1).subs. If ref(1).type='{}', the output will be the
%           element(s) of obj.ops specified by ref(1).subs, as for a cell
%           array;
%
% NOTES
% - '()'-type subsref is not supported for operators consisting of
%   multiple terms, or with a Kronecker product (otype==2) of multiple
%   factors along a dimension of size larger than 1.
% - C.dims.in and C.dims.out return C.dims(:,2) and C.dims(:,1), the
%   dimensions of the inputs and outputs of the factors of C.
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
% DJ, 10/22/2026: Allow '()' subsref for single-factor Kronecker products.

if isempty(ref)
    out = obj;
    return
end

switch ref(1).type
    case '.'
        if strcmp(ref(1).subs,'dims') && numel(ref)>=2 && strcmp(ref(2).type,'.')...
                && any(strcmp(ref(2).subs,{'in','out'}))
            % C.dims.in and C.dims.out return the input and output dimensions
            % of the factors, the second and first column of C.dims
            if strcmp(ref(2).subs,'in')
                out = obj.dims(:,2);
            else
                out = obj.dims(:,1);
            end
            % Perform remaining subsref
            if numel(ref)>=3
                out = subsref(out,ref(3:end));
            end
        else
            % Just use the built-in function
            out = builtin('subsref',obj,ref);
        end
        
    case '()'
        % When calling C(i,j), just extract C.ops(i,j);
        % With a single factor, otype is immaterial
        if size(obj.ops,2)>1 && ((obj.otype(1)==2 && size(obj,1)>1) || (obj.otype(2)==2 && size(obj,2)>1))
            error("'()'-type subsref is currently not supported for Kronecker products in 'tensopvar' objects.")
        end
        matdim = obj.matdim;
        if isscalar(ref(1).subs)
            % % Distinguish case of linear indexing
            error("Linear indexing of 'tensopvar' objects is not supported.")
        
        elseif length(ref(1).subs)==2
            % % Consider case of two indices
            indr = ref(1).subs{1};
            indc = ref(1).subs{2};
            
            % Allow for ':' for extracting all rows/columns
            if strcmp(indr,':')
                indr = 1:matdim(1);                                        
            end
            if strcmp(indc,':')
                indc = 1:matdim(2);                                        
            end

            % Error checking
            if any(indr<=0) || any(indc<=0)
                error('Indices must be positive integers');
            elseif any(indr>matdim(1))||any(indc>matdim(2))
                error('Indices exceed matrix dimensions');
            end
            if size(obj.ops,1)>1
                error("Subsindexing for tensor-PI operators consisting of multiple components is not supported.")
            end

            % Extract elements of the PI or integral operators
            % corresponding to the desired row and column numbers
            out = obj;
            for j=1:numel(out.ops)
                op_j = obj.ops{j};
                op_j = op_j(indr,indc);
                out.ops{j} = op_j;
            end
            out.dims(:,1) = numel(indr);             
            out.dims(:,2) = numel(indc);

            % Perform remaining subsref
            if numel(ref)>=2
                out =subsref(out,ref(2:end));
            end
            
        else
            % In case someone specifies more than two arguments
            error("At most two indices are supported for 'tensopvar' class objects.")
        end
        
    case '{}'
        % When calling C{i,j}, just extract C.ops{i,j};
        [ncomps,ntrms] = size(obj.ops);
        if isscalar(ref(1).subs)
            % Support only a single linear index
            if ~isnumeric(ref(1).subs{1}) || ~isscalar(ref(1).subs{1})
                error("Linear indexing of 'tensopvar' objects is not supported.")
            end
            [indr,indc] = ind2sub([ncomps,ntrms],ref(1).subs{1});
        elseif numel(ref(1).subs)==2
            % Extract row and column numbers
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
        % Return obj.ops{indr,indc} directly, matching cell '{}' semantics
        ref(1).subs = {indr,indc};
        out = builtin('subsref',obj.ops,ref(1));

        % Perform remaining subsref, e.g. for C{i,j}.property
        if numel(ref)>=2 
            out = subsref(out,ref(2:end)); 
        end
end
end
