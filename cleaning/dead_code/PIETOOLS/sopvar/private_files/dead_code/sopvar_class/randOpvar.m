function out = randOpvar(vars_in,vars_out,dim,degree,density)
% This creates a randOpvar object with specified parameters. 
% out: L_2^dim(2)[vars_in] --> L_2^dim(1)[vars_out]
% where each parameter in P is a sparse polynomial object, quadPoly(),
% with density of coefficients as specified and degree in main variables si
% is degree(1) and degree in dummy variables (si_dum) is degree(2)
%
% Each parameter has the form 
% (I_{dim(1)} \otimes Z_{degree(1)}(s)^T) C (I_{dim(2)} \otimes Z_{degree(2)}(s_dum)) 

% create normalized domains
dom_in = repmat([0,1],numel(vars_in),1);
dom_out = repmat([0,1],numel(vars_out),1);
% find unique variables
vars = union(vars_out, vars_in);
n= numel(vars);

% cell dimensions for Pout.params. Initialize an empty cell
celldim = repmat(3,1,n);
[~,idx] = ismember(vars, intersect(vars_in,vars_out));
celldim(~idx) = 1;
if n == 0
    params = repmat({zeros(dim(1),dim(2))},1,1);
else
    params = repmat({zeros(dim(1),dim(2))},celldim);
end

% number of monomials on each side in quadPoly. Arbitrarily fixed at 10
nMons = [10,10];

% for each parameter, repeat
for i=1:numel(params)
%linear to multi index
Aidx = cell(1,numel(celldim));
[Aidx{:}] = ind2sub(celldim,i);

var_s = vars_out;   % left monomials must only have si that appear in range space
var_t = vars_in;    % right monomials must have "at most" the si_dum corresponding to
                    % si that appear in domain space.

% if sj appears in both range and domain and the term is a multiplier, then
% there cannot be any sj_dum in the parameter, so remove that from var_t
for j=1:n          
    if Aidx{j}==1 && celldim(j)==3
        var_t = setdiff(var_t,vars{j});
    end
end

var_t = strrep(var_t,'s','t');

% crreate a randquadPoly with the monomials chosen; ideally this should
% create parameters leading to well-defined opvars. For example, no
% variables corresponding to sj if there is no sj in output, and if sj is
% in input then only sj_dum appears making sure that sj is integrated out,
% etc.
params{i} = reduce(quadPoly.randquadPoly(dim,nMons,var_s,var_t,[degree,degree],density));
end
out = sopvar(vars_in,vars_out,dim,dom_in, dom_out,params);
end