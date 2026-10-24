function out = subsref(obj,ref)
% OUT = SUBSREF(OBJ,REF) extracts the elements of 'tensopmat' object OBJ
% specified by REF. This function is automatically called when calling e.g.
% C(i,j), C{i,j} or C.ops for a 'tensopmat' object C.
%
% INPUTS
% - obj:    'tensopmat' object of which to extract certain elements;
% - ref:    struct indicating which element/fields of obj to extract;
%
% OUTPUTS
% - out:    object representing the elements of 'obj' specified by 'ref'.
%           If ref(1).type='.', this will be out=obj.(ref(1).subs). If
%           ref(1).type='()', the output will be a 'tensopmat' object
%           including the rows and columns with the indices specified by
%           ref(1).subs. If ref(1).type='{}', the output will be the
%           block of obj.ops with the indices specified by ref(1).subs. If
%           ref(1).subs are not both scalar (index vectors, or ':' for all
%           blocks), the output will be the 'tensopmat' object including
%           these row- and column-blocks (like subsasgn for C{I,J}=RHS);
%
% NOTES
% - The indices of C{i,j} count the row- and column-blocks of obj.ops, not
%   the scalar rows and columns that C(i,j) indexes. The keyword 'end'
%   cannot be used: it evaluates to 1 for a 'tensopmat', so C{i,end} is
%   C{i,1}. Use size(C.ops,2) instead.
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
% DJ, 10/23/2026: Support ':' and index vectors in C{I,J}, returning a
%                   'tensopmat' of the selected blocks;

if isempty(ref)
    out = obj;
    return
end

switch ref(1).type
    case '.'
        % Just use the built-in function
        out = builtin('subsref',obj,ref);
        
    case '()'
        % When calling C(i,j), just extract C.ops(i,j);
        matdim = size(obj);
        if isscalar(ref(1).subs)
            % % Distinguish case of linear indexing
            error("Linear indexing of 'tensopmat' objects is not supported.")

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

            % Establish to which operator each row/column index belongs
            nr_arr = obj.dims.out;      nc_arr = obj.dims.in;
            nR = numel(nr_arr);         nC = numel(nc_arr);
            nnr_arr = cumsum([0;nr_arr(:)]);   
            nnc_arr = cumsum([0;nc_arr(:)]);
            rblcks = (indr(:)>nnr_arr(1:end-1)' & indr(:)<=nnr_arr(2:end)')*(1:nR)';
            cblcks = (indc(:)>nnc_arr(1:end-1)' & indc(:)<=nnc_arr(2:end)')*(1:nC)';
            
            % Split the indices per block in the new operator
            r_strt = find([true; rblcks(2:end)~=rblcks(1:end-1); true]);
            c_strt = find([true; cblcks(2:end)~=cblcks(1:end-1); true]);
            blck_nums_r = rblcks(r_strt(1:end-1));
            blck_nums_c = cblcks(c_strt(1:end-1));
            nR_new = numel(blck_nums_r);
            nC_new = numel(blck_nums_c);

            % Extract elements of the tensor-PI or integral operators
            % corresponding to the desired row and column numbers
            out = obj;
            out.ops = cell(nR_new,nC_new);
            out.depmat.out = out.depmat.out(blck_nums_r,:);
            out.depmat.in = out.depmat.in(blck_nums_c,:);
            out.dims.out = reshape(r_strt(2:end)-r_strt(1:end-1),[],1);
            out.dims.in = reshape(c_strt(2:end)-c_strt(1:end-1),[],1);
            for i=1:nR_new
                rblck_i = rblcks(r_strt(i));
                ridcs_i = indr(r_strt(i):r_strt(i+1)-1)-nnr_arr(rblck_i);
                for j=1:nC_new
                    cblck_j = cblcks(c_strt(j));
                    if ~isempty(obj.ops{rblck_i,cblck_j})
                        out.ops{i,j} = zeros(0,0);
                        cidcs_j = indc(c_strt(j):c_strt(j+1)-1)-nnc_arr(cblck_j);
                        out.ops{i,j} = obj.ops{rblck_i,cblck_j}(ridcs_i,cidcs_j);
                    end
                end
            end
            out.chkval = true;

            % Perform remaining subsref
            if numel(ref)>=2
                out =subsref(out,ref(2:end));
            end
            
        else
            % In case someone specifies more than two arguments
            error("At most two indices are supported for 'tensopmat' class objects.")
        end
        
    case '{}'
        % When calling C{i,j}, just extract C.ops{i,j};
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
        % A single block, C{i,j}, is returned as such, anything else as 'tensopmat'
        is_blck = isnumeric(indr) && isscalar(indr);
        is_blck = is_blck && isnumeric(indc) && isscalar(indc);
        % ':' selects all row/column-blocks (else compared as a char code below)
        if strcmp(indr,':'),    indr = 1:ncomps;    end
        if strcmp(indc,':'),    indc = 1:ntrms;     end
        if any(indr>ncomps)
            error("Row indices cannot exceed number of components.")
        elseif any(indc>ntrms)
            error("Column indices cannot exceed number of monomials.")
        end
        % Several blocks, C{I,J} or with ':': the 'tensopmat' of these blocks. Slicing obj.ops, and
        % dims/depmat as for '()', keeps the dimensions of blocks in an empty row or column
        if ~is_blck
            out = obj;
            out.ops = obj.ops(indr,indc);
            out.dims.out = reshape(obj.dims.out(indr),[],1);
            out.dims.in = reshape(obj.dims.in(indc),[],1);
            out.depmat.out = obj.depmat.out(indr,:);
            out.depmat.in = obj.depmat.in(indc,:);
            out.chkval = true;
            if numel(ref)>=2
                out = subsref(out,ref(2:end));
            end
            return
        end
        % Extract desired coefficient operators
        ref(1).type = '{}';
        ref(1).subs = {indr,indc};
        ops = subsref(obj.ops,ref);
        out = ops;     
end
end
