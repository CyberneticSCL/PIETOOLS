function [P,ZL,ZR] = dpvar2sdvar(D,vars)
% [P,ZL,ZR] = DPVAR2SDOPVAR(D,VARS) takes a dpvar object D, and returns a
% struct P with fields A, B, representing the same decision
% variable matrix as
%   D = (Im o ZL{1}(s1) o ... o ZL{M}(sM))^T unvec(P.A + P.B'*d)
%           (In o ZR{1}(t1) o ... o ZR{N}(tN))
% where d is the vector of decision variables (D.dvarname) and where
%   {s1,...,sM} = vars.out;     {t1,...,tN} = vars.in;
% P also has fields P.m = m*nZL1*...*nZLM and P.n = n*nZLN*...*nZLN;

if ~isa(D,'dpvar')
    error("Input object must be of type 'dpvar'.")
end
if nargin==1
    % Assume all independent variables are output variables
    vars = struct();
    vars.out = D.varname';
    vars.in = cell(1,0);
elseif ~isa(vars,'struct') || ~isfield(vars,'in') || ~isfield(vars,'out')
    error("Input and output variable names should be specified as a struct with fields 'in' and 'out'.")
end

% Decompose the dpvar object
[m,n] = size(D);
dvars = D.dvarname;       ndvars = numel(dvars);
degmat = D.degmat;        nZ = size(D.degmat,1);

% Declare a basis of monomials in the desired variables including all 
% monomials appearing in D
varnames1 = vars.out;           M = numel(varnames1);
varnames2 = vars.in;            N = numel(varnames2);
ZL = repmat({0},1,numel(varnames1));
ZR = repmat({0},1,numel(varnames2));

% Express the parameter in terms of full set of variables
[~,idcs1_new,idcs1_old] = intersect(varnames1,D.varname','stable');
[~,idcs2_new,idcs2_old] = intersect(varnames2,D.varname','stable');
degmat_full = zeros(size(D.degmat,1),M+N);
degmat_full(:,idcs1_new) = degmat(:,idcs1_old);
degmat_full(:,M+idcs2_new) = degmat(:,idcs2_old);

% Extract monomial degrees in common variables
for i=1:numel(idcs1_new)
    degs_i = unique(degmat_full(:,idcs1_old(i)));
    ZL{idcs1_new(i)} = unique([ZL{idcs1_new(i)}; degs_i]);
end
for i=1:numel(idcs2_new)
    degs_i = unique(degmat_full(:,idcs2_old(i)));
    ZR{idcs2_new(i)} = unique([ZR{idcs2_new(i)}; degs_i]);
end
   
% Determine how each of the monomial bases relates to the full bases 0:d
nZL_arr = ones(1,M);   ZL_maps = cell(1,M);
for i=1:M
    nZL_arr(i) = numel(ZL{i});
    is_mon = (0:max(ZL{i}))'==ZL{i}';
    ZL_maps{i} = is_mon*(1:nZL_arr(i))'; 
end
nZR_arr = ones(1,N);   ZR_maps = cell(1,N);
for i=1:N
    nZR_arr(i) = numel(ZR{i});
    is_mon = (0:max(ZR{i}))'==ZR{i}';
    ZR_maps{i} = is_mon*(1:nZR_arr(i))'; 
end
nnZL_arr = cumprod([nZL_arr,1],'reverse');
nnZR_arr = cumprod([nZR_arr,1],'reverse');

% Determine the row index and column index associated with each
% coefficient defining the dpvar object
[rridcs,ccidcs,vals] = find(D.C);
ridcs = ceil(rridcs(:)/(ndvars+1));            % row number
cidcs = ceil(ccidcs(:)/nZ);                    % column number
didcs = rridcs(:) - (ridcs-1)*(ndvars+1);      % decision variable index
Zidcs = ccidcs(:) - (cidcs-1)*nZ;              % monomial index

% Determine the left-monomial index associated with each coefficient
ZLidcs = ones(numel(rridcs),1);
for i=1:M
    % Check which monomial in variable i is considered
    Zi_vals = degmat_full(Zidcs,i);
    ZLidcs_i = ZL_maps{i}(Zi_vals+1);
    % Account for Kronecker product with other monomials
    ZLidcs = ZLidcs + (ZLidcs_i-1)*nnZL_arr(i+1);
end

% Determine the right-monomial index associated with each coefficient
ZRidcs = ones(numel(ccidcs),1);
for i=1:N
    % Check which monomial in variable i is considered
    Zi_vals = degmat_full(Zidcs,i+M);
    ZRidcs_i = ZR_maps{i}(Zi_vals+1);
    % Account for Kronecker product with other monomials
    ZRidcs = ZRidcs + (ZRidcs_i-1)*nnZR_arr(i+1);
end

% Determine which coefficient correspond to fixed term (no dvars)
isA = didcs==1;
isB = ~isA;

% Declare a coefficient matrix associated with the fixed term
rridcsA = (ridcs(isA)-1)*nnZL_arr(1) + ZLidcs(isA);
ccidcsA = (cidcs(isA)-1)*nnZR_arr(1) + ZRidcs(isA);
lidcsA = (ccidcsA-1)*nnZL_arr(1)*m + rridcsA;
valsA = vals(isA);
A = sparse(lidcsA,1,valsA,m*nnZL_arr(1)*n*nnZR_arr(1),1);
% We store A as m*n row vector

% Declare a coefficient matrix acting on the decision variable terms
rridcsB = (ridcs(isB)-1)*nnZL_arr(1) + ZLidcs(isB);
ccidcsB = (cidcs(isB)-1)*nnZR_arr(1) + ZRidcs(isB);
lidcsB = (ccidcsB-1)*nnZL_arr(1)*m + rridcsB;
valsB = vals(isB);
B = sparse(didcs(isB)-1,lidcsB,valsB,ndvars,m*nnZL_arr(1)*n*nnZR_arr(1));
% We store B rather than Bt, which is a nd x m*n matrix

% Return a structure
P = struct();
P.A = A;
P.B = B;
P.m = m*nnZL_arr(1);
P.n = n*nnZR_arr(1);
P.dvarname = dvars;

end