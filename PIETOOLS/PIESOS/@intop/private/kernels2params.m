function [params,Z,Zd] = kernels2params(kernels,intmat,vars)
% [PARAMS,Z,ZD] = KERNELS2PARAMS(KERNELS,INTMAT,VARS) expresses the cell
% KERNELS of kernels of type 'double', 'polynomial', or 'dpvar', defining
% an 'intop' object in the variables VARS with intervals of integration
% INTMAT, in terms of the same monomials and decision variables, returning
% the cell PARAMS of sparse coefficient matrices, the monomials Z, and the
% decision variables ZD, see the documentation of the 'intop' class.
%
% INPUTS
% - kernels:    cell of m x n kernels of type 'double', 'polynomial', or
%               'dpvar', with element {i1,...,iN} the kernel associated
%               with the integral specified by row ik of intmat{k}, for
%               each k in {1,...,N}. The kernels may only depend on the
%               dummy variables of the functional;
% - intmat:     1 x N cell, with element i a qi x 2*di array of integers
%               specifying the intervals of integration along variable i,
%               as for the field intmat of an 'intop' object. Only the
%               number di of dummy variables of each variable is used;
% - vars:       struct with field 'in' (as for the field vars of an
%               'intop' object), with vars.in a 1 x N cellstr specifying
%               the names of the spatial variables. The dummy variables of
%               variable vars.in{i} are named
%               vars.in{i}_dum1,...,vars.in{i}_dumdi;
%
% OUTPUTS
% - params:     cell of the same size as kernels, with element {i1,...,iN}
%               the sparse m*(nd+1) x n*nZ coefficient matrix of the kernel
%               kernels{i1,...,iN}, for nd the number of decision variables
%               and nZ the number of monomials in the basis Z;
% - Z:          struct with fields 'in' and 'out'. Z.in is the 1 x N cell
%               with element i the nmons(i) x di array of monomial degrees
%               along variable i, and Z.out is an empty 1 x 0 cell;
% - Zd:         nd x 1 cellstr specifying the names of the decision
%               variables;
%
% NOTES
% - The monomials are the union of the rows of the degree matrices
%   restricted to the dummy variables of each variable, so that the basis
%   is a tensor product of the bases of the different variables, but not of
%   the different dummy variables of one variable.
% - 'double' kernels are expressed in terms of the monomial of degree zero.
% - The decision variables are the union of those of the different kernels,
%   in order of first appearance.
%

% DJ, 10/15/2026: Initial coding;

N = numel(vars.in);
deps = zeros(1,N);
dumnames = cell(1,0);
for i=1:N
    deps(i) = size(intmat{i},2)/2;
    for j=1:deps(i)
        dumnames{end+1} = [vars.in{i},'_dum',num2str(j)]; %#ok<AGROW>
    end
end
dir_of_var = repelem(1:N,deps);                 % direction of each dummy variable
K = numel(kernels);

% Check that the kernels are of appropriate type and dimension
msz = zeros(K,2);
for k=1:K
    if ~isa(kernels{k},'dpvar') && ~isa(kernels{k},'polynomial') && ~isa(kernels{k},'double')
        error("Kernels must be specified as objects of type 'double', 'polynomial', or 'dpvar'.")
    end
    msz(k,:) = size(kernels{k});
end
if K==0
    error("At least one kernel must be specified.")
elseif any(any(msz~=msz(1,:)))
    error("Matrix dimensions of the different kernels must match.")
end
m = msz(1,1);   n = msz(1,2);

% Convert all kernels to triplets, with decision variables and monomials
% expressed per kernel
trp = cell(K,1);
Zin_all = cell(K,N);
Zd_all = cell(K,1);
for k=1:K
    P = kernels{k};
    t = struct();
    if isa(P,'double')
        % Constant kernel: monomial of degree zero
        [r,c,v] = find(P);
        r = r(:);   c = c(:);   v = v(:);
        t.row = r;      t.j = c;        t.z = ones(size(r));    t.d = ones(size(r));
        nZ_k = 1;
        dg = zeros(1,numel(dumnames));
        Zd_all{k} = cell(0,1);
    else
        if isa(P,'dpvar')
            nd_k = numel(P.dvarname);
            [r,c,v] = find(P.C);
            r = r(:);   c = c(:);   v = v(:);   % row vectors if C has one row
            nZ_k = size(P.degmat,1);
            blk = ceil(r/(nd_k+1));
            t.d = r - (blk-1)*(nd_k+1);         % 1: constant term, l+1: decision variable l
            t.row = blk;                        % matrix row
            t.j = ceil(c/nZ_k);                 % matrix column
            t.z = c - (t.j-1)*nZ_k;             % monomial
            Zd_all{k} = P.dvarname(:);
        else
            % polynomial: monomial rows, columns (j-1)*m+i
            [z,c,v] = find(P.C);
            z = z(:);   c = c(:);   v = v(:);   % row vectors if C has one row
            nZ_k = size(P.degmat,1);
            t.row = mod(c-1,m)+1;
            t.j = ceil(c/m);
            t.z = z;
            t.d = ones(size(z));
            Zd_all{k} = cell(0,1);
        end
        % Express the monomials in terms of the full list of dummy variables, ignoring
        % variables that do not appear in any monomial
        dg = zeros(nZ_k,numel(dumnames));
        if nZ_k>0 && ~isempty(P.varname)
            is_used = full(any(P.degmat,1));
            [is_dum,loc] = ismember(P.varname(is_used),dumnames);
            if ~all(is_dum)
                error("Kernels must only depend on the dummy variables of the functional.")
            end
            dg(:,loc) = full(P.degmat(:,is_used));
        end
    end
    t.v = v;    t.nZ = nZ_k;
    for i=1:N
        Zin_all{k,i} = dg(:,dir_of_var==i);
    end
    trp{k} = t;
end

% Declare the unified monomial bases
Zin = cell(1,N);
for i=1:N
    Zi = unique(vertcat(Zin_all{:,i}),'rows');
    if isempty(Zi)
        Zi = zeros(1,deps(i));
    end
    Zin{i} = Zi;
end
[~,strides,nZ] = Zdims(Zin);
Z = struct('in',{Zin},'out',{cell(1,0)});

% Declare the unified list of decision variables
Zd = Zd_all{1};
if ~all(cellfun(@(a) isequal(a,Zd),Zd_all))
    Zd = unique(vertcat(Zd_all{:}),'stable');
end
nD = numel(Zd);

% Build the coefficient matrices
params = cell(size(kernels));
for k=1:K
    t = trp{k};
    % Map monomials to unified basis
    zmap = ones(max(t.nZ,1),1);
    for i=1:N
        [~,loc] = ismember(Zin_all{k,i},Zin{i},'rows');
        zmap(1:numel(loc)) = zmap(1:numel(loc)) + (loc-1)*strides(i);
    end
    % Map decision variables to unified list
    if isequal(Zd_all{k},Zd)
        pos = (1:nD)';
    else
        [~,pos] = ismember(Zd_all{k},Zd);
    end
    d_new = ones(size(t.d));
    is_d = t.d>1;
    d_new(is_d) = 1 + pos(t.d(is_d)-1);
    rows_new = (t.row-1)*(nD+1) + d_new;
    cols_new = (t.j-1)*nZ + zmap(t.z);
    params{k} = sparse(rows_new,cols_new,t.v,m*(nD+1),n*nZ);
end

end
