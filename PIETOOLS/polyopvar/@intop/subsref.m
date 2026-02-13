function out = subsref(obj,ref)
% OUT = SUBSREF(OBJ,S) extracts the elements of 'intop' object OBJ
% specified by REF. This function is automatically called when calling e.g.
% K(i,j) or K.params for an 'intop' object K.
% 
% INPUT
% - obj:    'intop' object of which to extract certain elements;
% - ref:    struct indicating which element/fields of obj to extract;
% 
% OUTPUT
% - out:    object representing the elements of 'obj' specified by 'ref'.
%           If ref.type='.', this will be out=obj.(ref.subs). If 
%           ref.type='()', the output will be a intop object including the 
%           elements in the indices specified by ref.subs;
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
        % When calling C(ridcs,cidcs), extract elements ridcs, cidcs from 
        % the matrix-valued operator;
        matdim = obj.matdim;
        n_trms = size(obj.omat,1);
        if isscalar(ref(1).subs)
            % % Distinguish case of linear indexing
            error("Linear indexing of 'intop' objects is not supported.")
        
        elseif length(ref(1).subs)==2
            % % Consider case of two indices
            indr = ref(1).subs{1};
            indc = ref(1).subs{2};
            
            % Allow for ':' for extracting all rows
            if strcmp(indr,':')
                indr = 1:matdim(1);
            end
            % Determine columns in parameter array associated with the
            % desired columns of the matrix-valued operator
            if strcmp(indc,':')
                indc = 1:size(matdim(2),2);
            elseif length(indc)>=2
                indc = repelem(matdim(2)*(0:n_trms-1)',numel(indc),1) + repmat(indc(:),n_trms,1);
            elseif isscalar(indc)
                indc = (indc(1)-1)*n_trms+1:(indc(1))*n_trms;
            else
                indc = [];
            end
            
            % Error checking
            if any(indr<=0) || any(indc<=0)
                error('Indices must be positive integers');
            elseif any(indr>matdim(1))||any(indc>matdim(2)*n_trms)
                error('Indices exceed matrix dimensions');
            end
            
            % Extract parameters corresponding to desired elements of
            % matrix-valued operator
            Kparam = obj.params;
            Kparam = Kparam(indr,indc);
                        
            % Build the new operator
            out = intop(Kparam,obj.omat,obj.pvarname,obj.dom);
            
        else
            % In case someone specifies more than two arguments            
            if all(cell2mat(ref(1).subs(3:end))==1)
                ref(1).subs = ref(1).subs(1:2);
                out = subsref(obj,ref);
            else            
                error("At most two indices are supported for 'intop' class objects.")
            end 
        end
        
    case '{}'
        error("'{}' type indexing of 'tensopvar' objects is not supported.")
end
end
