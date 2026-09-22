function [A,B] = unpack_sheets(C,mrow,ncol,nsheet)
% Inverse of 'pack_sheets' applied to the output of 'int_semisep': the
% column index of the result is (sheet,ncol,NR) with the semiseparable index
% fastest, so each sheet still occupies a contiguous block of ncol columns.
% Sheet 1 is returned as the vectorized A, the remaining nsheet sheets as
% the rows of B.
%
% The whole array is scanned once and then partitioned by column, rather
% than sliced first. Slicing a sparse matrix whose column count carries the
% decision variable dimension, as C(:,ncol+1:end) does, copies that block
% before anything useful happens, and the copy dominates once there are many
% decision variables.

if size(C,1)~=mrow || size(C,2)~=(1+nsheet)*ncol
    error("Internal error: unexpected coefficient dimensions.")
end

[irow,icol,val] = find(C);
is_A = icol<=ncol;

A = sparse((icol(is_A)-1)*mrow+irow(is_A),1,val(is_A),mrow*ncol,1);

jcol = icol(~is_A)-ncol;
sg = floor((jcol-1)/ncol)+1;
cc = jcol - (sg-1)*ncol;
B = sparse(sg,(cc-1)*mrow+irow(~is_A),val(~is_A),nsheet,mrow*ncol);

end


