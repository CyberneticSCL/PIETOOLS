function out = randOpvar(vars_in,vars_out,dim,degree,density)
vars = union(vars_out, vars_in);
n= numel(vars);
celldim = repmat(3,1,n);
[~,idx] = ismember(vars, intersect(vars_in,vars_out));
celldim(~idx) = 1;
if n == 0
    params = repmat({zeros(dim(1),dim(2))},1,1);
else
    params = repmat({zeros(dim(1),dim(2))},celldim);
end

nMons = [10,10];

for i=1:numel(params)
Aidx = cell(1,numel(celldim));
[Aidx{:}] = ind2sub(celldim,i);

var_s = vars_out;
var_t = vars_in;

for j=1:n
    if Aidx{j}==1 && celldim(j)==3
        var_t = setdiff(var_t,vars{j});
    end
end

var_t = strrep(var_t,'s','t');

params{i} = quadPoly.randquadPoly(dim,nMons,var_s,var_t,[degree,degree],density);
end
out = sopvar(vars_in,vars_out,dim,params);
end