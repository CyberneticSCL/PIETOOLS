function CG = pack_sheets(V,mrow,ncol)
% Reshape each row of V into an mrow x ncol matrix and lay the results out
% side by side, so that sheet sg occupies columns (sg-1)*ncol+1 : sg*ncol.

nsheet = size(V,1);
[sg,lin,val] = find(V);
irow = mod(lin-1,mrow)+1;
icol = floor((lin-1)/mrow)+1;
CG = sparse(irow,(sg-1)*ncol+icol,val,mrow,nsheet*ncol);

end


