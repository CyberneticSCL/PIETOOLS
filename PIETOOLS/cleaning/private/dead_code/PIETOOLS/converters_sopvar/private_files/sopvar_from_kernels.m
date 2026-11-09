function P = sopvar_from_kernels(K,S2,S3,S1,vars,dom,dims)
% P = SOPVAR_FROM_KERNELS(K,S2,S3,S1,VARS,DOM,DIMS) builds the 'sopvar'
%   (Pu)(S2,S3) = sum_alpha int_{S1} int_{S3dum} K_alpha(S2,S3,S3dum,S1dum) u(S3dum,S1dum),
% for kernels K_alpha specified as 'polynomial' matrices.
%
% INPUTS
% - K:      Cell array of q x p 'polynomial' matrices, one per parameter of
%           the 'sopvar'. If S3 contains n3 variables, K must be of size
%           [3,3,...,3] (n3 times; [1,1] for n3=0, [3,1] for n3=1), where
%           index 1 along the dimension of a variable in S3 is the
%           multiplier, index 2 the integral from the lower boundary to the
%           variable, and index 3 the integral from the variable to the
%           upper boundary, as for the parameters R0,R1,R2 of an 'opvar'.
%           Variables are indexed along the dimensions in order of their
%           NAME (S3 must be sorted by name).
%           Kernels use the primary variables for S2, both primary and
%           dummy variables for S3, and the dummy variables for S1;
% - S2,S3,S1:   Index vectors into the rows of VARS, of the variables that
%           appear only in the output, in both, or only in the input.
%           Each must be sorted by variable name;
% - vars:   nv x 2 array of 'polynomial' variables, specifying the primary
%           and dummy variable of each spatial direction;
% - dom:    nv x 2 array of the domain of each variable;
% - dims:   [q,p] size of the operator;
%
% OUTPUT
% - P:      The 'sopvar' object;
%
% NOTE: local helper of the sopvar converters, building the operator in the
% same way as 'opvar2d2sopvar', for an arbitrary number of variables.
%
% PIETOOLS - sopvar_from_kernels
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding: Claude 10/28/2026

names = pvar2varname(vars(:,1));    names = names(:)';
dnames = pvar2varname(vars(:,2));   dnames = dnames(:)';
for idx = {S2,S3,S1}
    if ~issorted(names(idx{1}))
        error('sopvar_from_kernels:order','Variable groups must be sorted by name.')
    end
end
n3 = numel(S3);
if n3==0,       sz = [1,1];
elseif n3==1,   sz = [3,1];
else,           sz = 3*ones(1,n3);
end
if ~isequal(size(K),sz) && ~(n3==2 && isequal(size(K),[3,3]))
    error('sopvar_from_kernels:size','K has the wrong size for the number of shared variables.')
end

out = [S2(:)',S3(:)'];      in = [S3(:)',S1(:)'];
sv = struct();
sv.out = names(out);        sv.in = names(in);
sdom = struct();
sdom.out = dom(out,:);      sdom.in = dom(in,:);
if isempty(out),    sv.out = {};    sdom.out = zeros(0,2);  end
if isempty(in),     sv.in = {};     sdom.in = zeros(0,2);   end

var1 = names(out);          var2 = dnames(in);
Kc = cell(numel(K),1);
for i=1:numel(K),   Kc{i} = polynomial(K{i});   end
[C,ZL,ZR] = get_params(Kc,var1,var2,dims);
params = reshape(C,sz);
P = sopvar(params,sv,ZL,ZR,sdom,dims);
end


function [C,Zs,Zt] = get_params(Q,var1,var2,dims)
% Coefficients of the polynomial matrices in Q in terms of a common basis.
nvar1 = numel(var1);        nvar2 = numel(var2);
Ccell = cell(size(Q));      Zscell = cell(size(Q));     Ztcell = cell(size(Q));
for i=1:numel(Q)
    Qi = quadPoly.polynomial2quadPoly(Q{i}, var1, var2);
    Ccell{i} = full(Qi.C);
    Zscell{i} = Qi.Zs;      Ztcell{i} = Qi.Zt;
end
if isscalar(Q)
    C = Ccell;      Zs = Zscell{1};     Zt = Ztcell{1};
    return
end
Zsfull = cell(1,nvar1);     Ztfull = cell(1,nvar2);
nZs = zeros(numel(Q),nvar1);    nZt = zeros(numel(Q),nvar2);
for i=1:numel(Q)
    for j=1:nvar1
        nZs(i,j) = size(Zscell{i}{j},1);
        Zsfull{j} = [Zsfull{j}; Zscell{i}{j}];
    end
    for j=1:nvar2
        nZt(i,j) = size(Ztcell{i}{j},1);
        Ztfull{j} = [Ztfull{j}; Ztcell{i}{j}];
    end
end
Zs = cell(1,nvar1);     Ps = 1;
Zt = cell(1,nvar2);     Pt = 1;
for j=1:nvar1
    [Ps_j,Zs_j] = uniquerows_integerTable(Zsfull{j});
    Zs{j} = Zs_j;       Ps = kron(Ps,Ps_j);
end
for j=1:nvar2
    [Pt_j,Zt_j] = uniquerows_integerTable(Ztfull{j});
    Zt{j} = Zt_j;       Pt = kron(Pt,Pt_j);
end
nnZs = [zeros(1,nvar1);cumsum(nZs,1)];
nnZt = [zeros(1,nvar2);cumsum(nZt,1)];
stride1 = [cumprod(nnZs(end,:),2,'reverse'),1];
stride2 = [cumprod(nnZt(end,:),2,'reverse'),1];
idcs1 = cell(numel(Q),1);   idcs2 = cell(numel(Q),1);
for i=1:numel(Q)
    idcs1{i} = 1;           idcs2{i} = 1;
    for j=1:nvar1
        idcs1{i} = idcs1{i}(:)' + stride1(j+1)*(nnZs(i,j):nnZs(i+1,j)-1)';
        idcs1{i} = idcs1{i}(:)';
    end
    for j=1:nvar2
        idcs2{i} = idcs2{i}(:)' + stride2(j+1)*(nnZt(i,j):nnZt(i+1,j)-1)';
        idcs2{i} = idcs2{i}(:)';
    end
end
C = Ccell;
m = dims(1);    n = dims(2);
for i=1:numel(C)
    Ps_i = spIkron(m,Ps(idcs1{i},:));
    Pt_i = spIkron(n,Pt(idcs2{i},:));
    C{i} = full(Ps_i'*Ccell{i}*Pt_i);
end
end
