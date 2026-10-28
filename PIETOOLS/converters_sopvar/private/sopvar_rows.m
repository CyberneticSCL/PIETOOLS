function Q = sopvar_rows(P,rows)
% Q = SOPVAR_ROWS(P,ROWS) returns the 'sopvar' operator that maps to only
% the components ROWS of the output of the 'sopvar' P.
nz = prod(cellfun(@numel,P.ZL));
if isempty(P.ZL),   nz = 1;     end
rowsel = reshape((rows(:)'-1)*nz + (1:nz)',[],1);
params = P.params;
for m=1:numel(params)
    params{m} = params{m}(rowsel,:);
end
Q = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,[numel(rows),P.dims(2)]);
end
