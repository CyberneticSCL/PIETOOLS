function key = inttype2key(alpha)
% INTTYPE2KEY Takes a cell of arrays specifying what type of integral is
% taken, and converts it tp a unique integer key, identifying this same
% integral in the polyopvar structure
arguments (Input)
    alpha cell
end

% Determine the number of spatial variables
N = numel(alpha);
% Determine the maximal monomial degree
maxdeg = zeros(N,2);
for ii=1:N
    maxdeg(ii,1) = size(alpha{ii},1);
    maxdeg(ii,2) = size(alpha{ii},2);
end

% Determine for each variable how many possible values of alpha{ii} exist
n_alpha = 3.^(maxdeg(:,1).*maxdeg(:,2));
n_alpha_tot = cumprod([1;n_alpha(1:end-1)]);

% Determine for each variable which of the n_alpha(ii) options of alpha is
% specified
idcs  = zeros(1,N);
for ii=1:N
    alpha_ii = alpha{ii}(:)+2;
    if isempty(alpha_ii)
        idcs(ii) = 1;
        continue
    end
    n_alpha_tot_ii = 3.^(0:numel(alpha_ii)-1);
    idcs(ii) = n_alpha_tot_ii*(alpha_ii-1)+1;
end

% Set the key
key = (idcs-1)*n_alpha_tot + 1;


end