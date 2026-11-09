% Test_tensopvar_subsref.m -- tests @tensopvar/subsref.m '()' indexing.
% Odd trials (otype=[1,2]): row indexing Aop(ridcs,:) is checked against
% the rows of apply_op(Aop,v), and the column and out-of-range indexing
% against the expected sizes and errors. Even trials (otype=[2,2]): every
% '()' index of an object with row dimension >1 must error with the
% Kronecker-product message.
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/21/2026: Add single-factor trials (otype 2 allowed there);

% NOTES
% - subsref disallows '()' indexing along a Kronecker (otype==2)
%   dimension of size >1, and apply_op disallows otype(2)==1, so the
%   functional test uses every factor with column dimension 1: the
%   overall column dimension (a Kronecker product of 1's) is then 1, and
%   row indexing is unrestricted since otype(1)==1.

ntrials = 20;
for trial=1:ntrials
if mod(trial,2)==1
tol = 1e-8;

d = randi(3)+1;    % at least 2 factors
N = randi(3);
type = [1,2];
dnsty = 0.5;
p = 2;
pdegv = 2;

m = randi(4)+1;    % at least 2 rows, so row-indexing is non-trivial
dims = [repmat(m,d,1), ones(d,1)];  % every factor: m x 1
vars = struct();
vars.in = cell(1,N);
idcs = sort(randperm(2*N,N));
for i=1:N
    vars.in{i} = ['s',num2str(idcs(i))];
end
vars.out = cell(1,0);
dom = struct();
dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];
dom.out = zeros(0,2);
degs = struct('in',randi(2,[1,N]),'out',zeros(1,0));
depmat = struct('in',logical(randi([0,1],d,N)),'out',false(d,0));

Aop = rand_tensopvar(dims,vars,dom,depmat,1,type,degs,dnsty);   % nterms=1
if size(Aop,2)~=1
    error("Trial %d: test setup error: expected overall column dimension 1.",trial)
end

v = cell(1,d);
for j=1:d
    pvars_j = polynomial(Aop.ops{1,j}.vars.in(:));
    n_j = size(Aop.ops{1,j},2);
    v{j} = rand_poly([n_j,p],pvars_j,pdegv,50);
end
Av = apply_op(Aop,v);

%%%% Row indexing (unrestricted since otype(1)==1)
ridcs = [2,1];
Asub_r = Aop(ridcs,:);
if ~isequal(size(Asub_r),[numel(ridcs),1])
    error("Trial %d: row indexing: size of Aop(ridcs,:) is incorrect.",trial)
end
Av_r = apply_op(Asub_r,v);
err_r = cleanpoly(Av_r - Av(ridcs,:),1e-10);
if max(max(abs(err_r.C)))>tol
    error("Trial %d: row indexing: residual exceeds tolerance.",trial)
end

%%%% Column indexing: only valid index is 1 (overall n==1), and ':'
Asub_c1 = Aop(:,1);
Asub_cc = Aop(:,':');
if ~isequal(size(Asub_c1),[m,1]) || ~isequal(size(Asub_cc),[m,1])
    error("Trial %d: column indexing: size of Aop(:,1) or Aop(:,':') is incorrect.",trial)
end

%%%% Out-of-range indices should still error
threw_error = false;
try
    Asub_bad = Aop(m+1,:); %#ok<NASGU>
catch
    threw_error = true;
end
if ~threw_error
    error("Trial %d: row indexing: out-of-range index should have thrown an error.",trial)
end

else
% otype(1)==2 and overall row dimension >1: every '()' index must error (subsref's restriction)
d2 = randi(3)+1;
N2 = randi(3);
type2 = [2,2];
dnsty2 = 0.5;
rdims = randi(2,[d2,1])+1;    % every factor row dim in {2,3}, so Kronecker product (otype(1)==2) exceeds 1
dims2 = [rdims, ones(d2,1)];  % column dim 1 per factor, as in the odd-trial case above
vars2 = struct();
vars2.in = cell(1,N2);
idcs2 = sort(randperm(2*N2,N2));
for i=1:N2
    vars2.in{i} = ['t',num2str(idcs2(i))];
end
vars2.out = cell(1,0);
dom2 = struct();
dom2.in = randi([-2,2],[N2,1]) +[zeros(N2,1),rand([N2,1])];
dom2.out = zeros(0,2);
degs2 = struct('in',randi(2,[1,N2]),'out',zeros(1,0));
depmat2 = struct('in',logical(randi([0,1],d2,N2)),'out',false(d2,0));
depmat2.in(1,~any(depmat2.in,1)) = true;   % every variable used at least once

Bop2 = rand_tensopvar(dims2,vars2,dom2,depmat2,1,type2,degs2,dnsty2);
if size(Bop2,1)<=1
    error("Trial %d: test setup error: expected overall row dimension >1.",trial)
end
threw_error2 = false;
try
    Bsub = Bop2(1,1); %#ok<NASGU>
catch ME2
    threw_error2 = true;
    if ~contains(ME2.message,'Kronecker')
        error("Trial %d: Bop2(1,1) errored, but not with the expected Kronecker-product message (got: %s).",trial,ME2.message)
    end
end
if ~threw_error2
    error("Trial %d: '()'-indexing of an otype(1)==2 object with row dimension >1 should have thrown an error.",trial)
end
end
end

%%%% Single factor: otype is immaterial, so rows and columns can be extracted also for otype=[2,2]
% (u=T*x with a vector-valued state), checked against apply_op: rows of A*v, and A(:,c)*v(c) = A*(v
% with the entries outside c zeroed)
nsingle = 0;
for trial=1:10
    tol = 1e-8;
    type = {[2,2],[1,2]};       type = type{mod(trial,2)+1};
    N = randi(2);
    m = randi(3)+1;             n = randi(3)+1;                 % at least 2 rows and columns
    vars = struct('in',{arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false)},'out',{cell(1,0)});
    dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));
    degs = struct('in',randi(2,[1,N]),'out',zeros(1,0));
    depmat = struct('in',true(1,N),'out',false(1,0));
    Aop = rand_tensopvar([m,n],vars,dom,depmat,1,type,degs,0.5);
    if size(Aop.ops,2)~=1 || ~isequal(size(Aop),[m,n])
        error("Single factor, trial %d: test setup error.",trial)
    end
    pvars = polynomial(vars.in(:));
    v = rand_poly([n,2],pvars,2,50);
    Av = apply_op(Aop,{v});
    ridcs = randperm(m,2);
    cidcs = sort(randperm(n,2));
    % rows
    Asub_r = Aop(ridcs,':');
    if ~isequal(size(Asub_r),[2,n])
        error("Single factor, trial %d (otype [%d,%d]): size of Aop(ridcs,:) is incorrect.",trial,type)
    end
    err = cleanpoly(apply_op(Asub_r,{v}) - Av(ridcs,:),1e-10);
    if max(max(abs(err.C)))>tol
        error("Single factor, trial %d (otype [%d,%d]): row indexing: residual exceeds tolerance.",trial,type)
    end
    % columns
    Asub_c = Aop(':',cidcs);
    if ~isequal(size(Asub_c),[m,2])
        error("Single factor, trial %d (otype [%d,%d]): size of Aop(:,cidcs) is incorrect.",trial,type)
    end
    vz = polynomial(diag(double(ismember(1:n,cidcs))))*v;
    err = cleanpoly(apply_op(Asub_c,{v(cidcs,:)}) - apply_op(Aop,{vz}),1e-10);
    if max(max(abs(err.C)))>tol
        error("Single factor, trial %d (otype [%d,%d]): column indexing: residual exceeds tolerance.",trial,type)
    end
    nsingle = nsingle+1;
end

fprintf('Test_tensopvar_subsref passed (%d trials, otype(1) in {1,2} alternating; %d single-factor trials).\n',ntrials,nsingle)
