function out = randsopvar(vars_S1,vars_S2,vars_S3,dim,degree,density)
% This creates a randOpvar object with specified parameters.
% out: L_2^dim(2)[vars_in] --> L_2^dim(1)[vars_out]
% where each parameter in P is a sparse polynomial object, quadPoly(),
% with density of coefficients as specified and degree in main variables si
% is degree(1) and degree in dummy variables (si_dum) is degree(2)
%
% Each parameter has the form
% (I_{dim(1)} \otimes Z_{degree(1)}(s)^T) C (I_{dim(2)} \otimes Z_{degree(2)}(s_dum))

% create normalized domains
dom_1 = repmat([0,1],numel(vars_S1),1);
dom_2 = repmat([0,1],numel(vars_S2),1);
dom_3 = repmat([0,1],numel(vars_S3),1);
% find unique variables
% n1= numel(vars_S1);
% n2= numel(vars_S2);
n3= numel(vars_S3);

% cell dimensions for Pout.params. Initialize an empty cell
celldim = repmat(3,1,n3);
%[~,idx] = ismember(vars, intersect(vars_in,vars_out));
%celldim(~idx) = 1;
if n3 == 0
    params = repmat({zeros(dim(1),dim(2))},1,1);
else
    params = repmat({zeros(dim(1),dim(2))},celldim);
end

% number of monomials on each side in quadPoly. Arbitrarily fixed at 10
nMons = [10,10];

var_s = {vars_S2, vars_S3};   % left monomials use S2 and S3 -- common to all terms

% for each parameter, repeat
for i=1:numel(params)
    %linear to multi index
    Aidx = cell(1,numel(celldim));
    [Aidx{:}] = ind2sub(celldim,i); % position of element i
    vars_dum=vars_S3;
    for j=1:n3
        if Aidx{j}~=1
            vars_dum{j} = [vars_dum{j},'_dum'];    %create list of dummy variables, but only if alpha_j neq 1
        end
    end
    var_t = {vars_dum, vars_S1};   % right monomials use S1 and S3_dum

    % crreate a randquadPoly with the monomials chosen; ideally this should
    % create parameters leading to well-defined opvars. For example, no
    % variables corresponding to sj if there is no sj in output, and if sj is
    % in input then only sj_dum appears making sure that sj is integrated out,
    % etc.
    params{i} = quadPoly.randquadPoly(dim,nMons,var_s,var_t,[degree,degree],density);
%    params{Aidx} = reduce(quadPoly.randquadPoly(dim,nMons,var_s,var_t,[degree,degree],density));
end
%out = sopvar(vars_S1,vars_out,dim,dom_in, dom_out,params);
out = sopvar(vars_S3,dom_3,dim,params,vars_S1,dom_1,vars_S2,dom_2);
