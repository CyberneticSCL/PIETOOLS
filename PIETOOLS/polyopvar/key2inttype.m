function alpha = key2inttype(maxdeg,key)
% KEY2INTTYPE determines for a given max distributed monomial degree MAXDEG
% the type of integral associated with the integer value KEY
%
% INPUT
% - maxdeg: Nx2 array specifying for each of the N variables the 
%           distributed monomial degree taken along the input variable
%           var2(ii) (maxdeg(ii,2)) and output variable var1(ii)
%           (maxdeg(ii,1));
% - key:    scalar integer uniquely determining what type of integral 

% OUTPUT
% - alpha:  1xN cell with each element ii a maxdeg(ii,1) x maxdeg(ii,2) 
%           array with elements ranging from 1 to 3, specifying what type
%           of integral must be taken.
%
% NOTES:
% Consider variable s1 with associated dummy variable t1. If
% maxdeg(ii,1)=p and maxdeg(ii,2)=q, we consider a map from the qth
% distributed monomial along variable s1,
%   x(s1_1,...,s1_q) = x(s1_1)...x(s1_q),
% to the pth distributed monomial along variable s1,
%   x(s1_1,...,s1_p) = x(s1_1)...x(s1_p).
% We define such a map by integral operators
%   int_{a1}^{b1}...int_{a1}^{b1} I_{alpha1}(s1_1,...,s1_p,t1_1,...,t1_q)
%       P_{alpha}(s1_1,...,s1_p,t1_1,...,t1_q) x(t1_1,...,t1_q) dt1
% where now alpha1=alpha{1} is a pxq matrix determining the limits of the
% integral. In particular, for each i in 1:p and j in 1:q, we have
%   int_{a1}^{b1} I_{alpha1(i,j)}(s1_i,t1_j) f(t1_j) dt1_j
%        { int_{a1}^{s1} f(t1_j) dt1_j,     alpha1(i,j)=-1;
%       ={ f(s1_i),                         alpha1(i,j)=0;
%        { int_{s1}^{b1} f(t1_j) dt1_j,     alpha1(i,j)=1;


arguments (Input)
    maxdeg      (:,2) {mustBeInteger}
    key     (1,1) {mustBeInteger}
end

% Determine the number of variables
N = size(maxdeg,1);
% For each variable, check how many possible matrices alpha{i} there are
n_alpha = 3.^(maxdeg(:,1).*maxdeg(:,2));
n_alpha_tot = cumprod([1;n_alpha]);
if key>n_alpha_tot(end)
    error("Value of key exceeds total number of possibe integrals.")
end
% Determine which of the possible options of alpha{i} in each i correspond
% to ind
alp_idcs = zeros(N,1);
ind_rem = key;
for ii=1:N-1
    alp_idx = ceil(ind_rem/n_alpha_tot(end-ii))-1;
    alp_idcs(N-ii+1) = alp_idx+1;
    ind_rem = ind_rem - alp_idx*n_alpha_tot(end-ii);
end
alp_idcs(1) = ind_rem;
% Build the corresponding value of alpha{i} for each i
alpha = cell(1,N);
for ii=1:N
    % Generate alpha as sz1 x sz2 array
    alpha_ii = zeros(maxdeg(ii,1),maxdeg(ii,2));
    if isempty(alpha_ii)
        alpha{ii} = alpha_ii;
        continue
    end
    % Eestablish the cumulative number of possible values in elements 1
    % through jj of alpha
    sz_ii = maxdeg(ii,1)*maxdeg(ii,2);
    n_alpha_ii = 3.^(0:sz_ii)';
    % Establish for each element the value associated with index alp_idx_ii
    alp_idx_ii = alp_idcs(ii);
    for jj=1:sz_ii-1
        alp_val_jj = ceil(alp_idx_ii/n_alpha_ii(end-jj))-1;
        alpha_ii(sz_ii-jj+1) = alp_val_jj+1;
        alp_idx_ii = alp_idx_ii - alp_val_jj*n_alpha_ii(end-jj);
    end
    alpha_ii(1) = alp_idx_ii;
    alpha{ii} = alpha_ii-2;
end


end