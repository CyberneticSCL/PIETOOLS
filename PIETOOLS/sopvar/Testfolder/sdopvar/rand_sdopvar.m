function Pout = rand_sdopvar(matdim,vars,dom,degs,ndecvars,dnsty)
% Generates a random sdopvar object given input parameters
% Each kernel coefficient satisfies
%
%   vec(C_alpha(d)) = params{alpha}.A + params{alpha}.B' * d
%
% where d is the column vector of decision variables.

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

% Declare decision variables
Zd = "coef_"+string(num2cell(1:ndecvars));

% Declare random coefficient matrices
% Dimensions of C_alpha
nL = matdim(1)*prod(nZL);
nR = matdim(2)*prod(nZR);
nC = nL*nR;

if nargin <= 5 || isempty(dnsty)
    dnsty = 2/(prod(nZL)*prod(nZR));
end
% Allocate parameter cells
Acell = cell([3*ones(1,N3),1,1]);
Bcell = cell([3*ones(1,N3),1,1]);

for k = 1:numel(Acell)
    % Determine whether alpha-index is a multiplier term
    if N3 == 0
        is_mult = false(1,0);
    else
        idcs_k = cell(1,N3+1);
        [idcs_k{:}] = ind2sub([3*ones(1,N3),1],k);
        is_mult = cellfun(@(a) a == 1,idcs_k(1:N3));
    end

    if ~any(is_mult)
        % no multiplier: all entries of C_alpha(d) are admissible
        A = sprand(nC,1,dnsty);
        B = sprand(ndecvars,nC,dnsty);
    else
        % Multiplier kernel: cannot have theta coefficients is some variables.
        cidcs = reshape(1:nR,fliplr([matdim(2),nZR]));

        is_mult_full = false(1,N);
        is_mult_full(idcs_S3_in) = is_mult;

        full_idcs = cellfun(@(a) 1:a, ...
            num2cell(size(cidcs)), ...
            'UniformOutput',false);

        full_idcs(fliplr([false,is_mult_full])) = {1};

        % Admissible columns of C_alpha
        cidcs = cidcs(full_idcs{:});
        cidcs = cidcs(:);

        nRk = numel(cidcs);
        nCk = nL*nRk;

        % Generate only on the admissible support
        Ak = sprand(nCk,1,dnsty);
        Bk = sprand(ndecvars,nCk,dnsty);

        % Map reduced vectorized coordinates into vec(C_alpha)
        [ia,~,va] = find(Ak);

        rowC = mod(ia-1,nL) + 1;
        colC = cidcs(floor((ia-1)/nL) + 1);
        ia_full = rowC + (colC-1)*nL;

        A = sparse(ia_full,1,va,nC,1);

        [ib,jb,vb] = find(Bk);

        rowC = mod(jb-1,nL) + 1;
        colC = cidcs(floor((jb-1)/nL) + 1);
        jb_full = rowC + (colC-1)*nL;

        B = sparse(ib,jb_full,vb,ndecvars,nC);
    end

    Acell{k} = A;
    Bcell{k} = B;
end

params = struct('A',{Acell},'B',{Bcell});
Pout = sdopvar(params,vars,Zd,ZL,ZR,dom,matdim);
end