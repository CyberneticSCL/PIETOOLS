function Pop = rand_sopvar(matdim,vars,dom,degs,dnsty)
% RAND_SOPVAR Generates a random sopvar object operator with monomials of 
% degree at most degs in each of the variables specified in vars.
%
% INPUTS
% - matdim: 1 x 2 array of integers, specifying the row and column
%           dimensions of the matrix-valued objects;
% - vars:   struct with fields 'in' and 'out'. Each field is a 1 x N (or 1
%           x M) cellstr, specifying the names of the input and output
%           variables of the operator. Can also be set to "vars.in=N" and
%           "vars.in=M" for the function to randomly generate N input and M
%           output variable names.
% - dom:    struct with fields 'in' and 'out', with each element a N x 2
%           (or M x 2) array specifying the intervals on which the spatial
%           variables are defined;
% - degs:   struct with fields 'in' and 'out', with each element a 1 x N
%           (or 1 x M) array specifying the maximal degree of the monomials
%           in the input and output variables. The true monomial vector
%           will be generated randomly based on this maximal degree;
% - dnsty:  scalar between 0 and 1, specifying the density of the sparse
%           coefficient matrices defining the operator parameters;
%
% OUTPUTS
% - Pop:    randomly generated sopvar object in the specified variables.


% Check that the matrix dimensions are properly specified
if isscalar(matdim)
    matdim = [matdim,matdim];
elseif numel(matdim)~=2
    error("Matrix dimensions of the operator should be specified as 1x2 array.")
end

% Check that the variables are properly specified
if ~isa(vars,'struct') || ~isfield(vars,'in') || ~isfield(vars,'out')
    error("Variables should be specified as struct with fields 'in' and 'out'.")
end
if isa(vars.in,'double')
    % Declare M=vars.in input variables
    N = vars.in;
    vars.in = cell(1,N);
    idcs = sort(randperm(2*N,N));
    for i=1:N
        vars.in{i} = ['s',num2str(idcs(i))];
    end
elseif ~iscellstr(vars.in)
    error("Input variables should be specified as 1 x N 'cellstr' object.")
end
vars.in = vars.in(:)';      N = numel(vars.in);
if isnumeric(vars.out) && isscalar(vars.out)
    % Declare M=vars.out output variables
    M = vars.out;
    vars.out = cell(1,M);
    idcs = sort(randperm(2*N+2*M,M));
    for i=1:M
        vars.out{i} = ['s',num2str(idcs(i))];
    end
elseif ~iscellstr(vars.out)
    error("Output variables should be specified as 1 x M 'cellstr' object.")
end
vars.out = vars.out(:)';    M = numel(vars.out);
[~,~,idcs_S3_in] = intersect(vars.out,vars.in);
N3 = numel(idcs_S3_in);

% Check that the domains are appropriately specified
if isnumeric(dom) && all(size(dom)==[1,2])
    dom1 = dom;
    dom = struct();
    dom.in = repmat(dom1,N,1);
    dom.out = repmat(dom1,M,1);
elseif ~isa(dom,'struct') || ~isfield(dom,'in') || ~isfield(dom,'out')
    error("Domains should be specified as struct with fields 'in' and 'out'.")
end
if all(size(dom.in)==[1,2])
    % Assume same domain for all input variables
    dom.in = repmat(dom.in,N,1);
elseif ~all(size(dom.in)==[N,2])
    error("Input domains should be specified as N x 2 array for N input variables.")
end
if all(size(dom.out)==[1,2])
    % Assume same domain for all output variables
    dom.out = repmat(dom.out,M,1);
elseif ~all(size(dom.out)==[M,2])
    error("Output domains should be specified as M x 2 array for M output variables.")
end

% A variable shared by the input and output space is ONE variable and       % MMP, 09/17/2026
% must carry ONE interval: delta(s-s') in a multiplier direction            % MMP, 09/17/2026
% identifies the two copies, and an integral direction takes its bounds     % MMP, 09/17/2026
% from a single interval. Two different intervals for one name is a         % MMP, 09/17/2026
% caller error, so it is REJECTED here rather than reconciled: adopting     % MMP, 09/17/2026
% either side silently would build an operator the caller did not ask       % MMP, 09/17/2026
% for, and the discarded interval would reappear as a wrong integration     % MMP, 09/17/2026
% limit. Agreement costs the caller nothing -- the scalar and 1x2 forms     % MMP, 09/17/2026
% above broadcast one row to both sides, and a caller passing               % MMP, 09/17/2026
% per-variable arrays should index a single domain table by variable        % MMP, 09/17/2026
% name, as both in-repo callers do. NaN matches NaN, as in                  % MMP, 09/17/2026
% 'canonical_var_order', so a placeholder domain is not a clash.            % MMP, 09/17/2026
vS3_dom = intersect(vars.in,vars.out);                                      % MMP, 09/17/2026
if ~isempty(vS3_dom)                                                        % MMP, 09/17/2026
    [~,i_in ] = ismember(vS3_dom,vars.in);                                  % MMP, 09/17/2026
    [~,i_out] = ismember(vS3_dom,vars.out);                                 % MMP, 09/17/2026
    dl = dom.in(i_in,:);        dr = dom.out(i_out,:);                      % MMP, 09/17/2026
    bad = find(any(dl~=dr & ~(isnan(dl)&isnan(dr)),2),1);                   % MMP, 09/17/2026
    if ~isempty(bad)                                                        % MMP, 09/17/2026
        % '+' on a string with NaN yields <missing>, which 'error' rejects, % MMP, 09/17/2026
        % so the endpoints go through %g rather than string concatenation.  % MMP, 09/17/2026
        error("Variable '%s' is shared by the input and output space but "...
              +"was given two different domains, [%g,%g] on the input "...
              +"side and [%g,%g] on the output side.", ...
              vS3_dom{bad},dl(bad,1),dl(bad,2),dr(bad,1),dr(bad,2))         % MMP, 09/17/2026
    end                                                                     % MMP, 09/17/2026
end                                                                         % MMP, 09/17/2026

% Check that the degrees are properly specified
if isnumeric(degs) && isscalar(degs)
    deg1 = degs;
    degs = struct();
    degs.in = repmat(deg1,1,N);
    degs.out = repmat(deg1,1,M);
elseif ~isa(degs,'struct') || ~isfield(degs,'in') || ~isfield(degs,'out')
    error("Monomial degrees should be specified as struct with fields 'in' and 'out'.")
end
if isscalar(degs.in)
    % Assume same degree for all input variables
    degs.in = repmat(degs.in,1,N);
elseif numel(degs.in)~=N
    error("Number of input monomial degrees should match number of input variables.")
end
if isscalar(degs.out)
    % Assume same degree for all output variables
    degs.out = repmat(degs.out,1,M);
elseif numel(degs.out)~=M
    error("Number of output monomial degrees should match number of output variables.")
end

% Declare the monomial vectors
ZR = cell(1,N);     nZR = zeros(1,N);
ZL = cell(1,M);     nZL = zeros(1,M);
for i=1:N
    ZR{i} = unique([0;randi([0,degs.in(i)],[degs.in(i)+1,1])]);
    nZR(i) = numel(ZR{i});
end
for i=1:M
    ZL{i} = unique(randi([0,degs.out(i)],[degs.out(i)+1,1]));
    nZL(i) = numel(ZL{i});
end

% Declare random coefficient matrices
S3 = intersect(vars.in,vars.out);
n_shared = numel(S3);
Ccell = cell([3*ones(1,n_shared),1,1]);
nL = matdim(1)*prod(nZL);
nR = matdim(2)*prod(nZR);
if nargin<=4 || isempty(dnsty)
    dnsty = 2/(prod(nZL)*prod(nZR));
end
for k=1:numel(Ccell)
    Ccell{k} = sprand(nL,nR,dnsty);
    % Determine if the parameter corresponds to a multiplier term
    if N3==0
        continue
    end
    idcs_k = cell(1,N3+1);
    [idcs_k{:}] = ind2sub([3*ones(1,N3),1],k);
    is_mult = cellfun(@(a) a==1,idcs_k(1:N3));
    if ~any(is_mult)
        continue
    end
    % Split column indices per monomials
    cidcs_rtn = 1:nR;
    cidcs_rtn = reshape(cidcs_rtn,fliplr([matdim(2),nZR]));
    % Determine which column indices correspond with constant monomial
    is_mult_full = false(1,N);
    is_mult_full(idcs_S3_in) = is_mult;
    full_idcs = cellfun(@(a) 1:a,num2cell(size(cidcs_rtn)),'UniformOutput',false);
    full_idcs(fliplr([false,is_mult_full])) = {1};
    cidcs_rtn = cidcs_rtn(full_idcs{:});
    cidcs_rtn = cidcs_rtn(:);
    % Maintain only indices associated with constant monomial
    [ridcsC,cidcsC] = find(Ccell{k});
    rtn_idcs = ismember(cidcsC,cidcs_rtn);
    vals_rtn = sprand(nnz(rtn_idcs),1,dnsty);% + vals(rtn_idcs);
    Ccell{k} = sparse(ridcsC(rtn_idcs),cidcsC(rtn_idcs),vals_rtn,nL,nR);
end

% Declare the sopvar object
Pop = sopvar(Ccell,vars,ZL,ZR,dom,matdim);                                  % MMP, 08/29/2026

end

