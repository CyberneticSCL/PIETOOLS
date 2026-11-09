% Test_tensopvar_reorder_vars.m -- tests @tensopvar/reorder_vars.m. The
% oracle is the permuted fields themselves: reordering the factors by
% 'order' must permute 'ops' (columns), 'dims' and
% 'depmat.in'/'depmat.out' (rows, one per factor) consistently. Cases:
%   (1) otype(2)==1 (pointwise columns): permutation, inverse
%       permutation, identity;
%   (2) single factor: no-op regardless of the column dimension;
%   (3) otype(2)==2 (Kronecker columns), multiple factors: clean error.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - A reordering for Kronecker columns with several factors is not
%   implemented (it needs the kernel arguments within each sopvar factor
%   permuted), so case (3) checks for the explicit error.

d = randi(3)+1;   % at least 2 factors, so factor order is meaningful
N = d;            % one input variable per factor, so depmat.in is diagonal-like
M = randi(2);

%%%% Case 1: otype(2)==1 (pointwise columns) -- reorder must actually permute
% Use otype(1)==2 (Kronecker rows) so per-factor row dims may differ,
% making a row-permutation of 'dims' and 'depmat' visibly distinct from a
% no-op -- otype(2)==1 forces per-factor column dims to match instead.
type = [2,1];
dims = [randi(4,[d,1]), repmat(randi(4),d,1)];
vars = struct();
vars.in = arrayfun(@(i)sprintf('x%d',i),1:N,'UniformOutput',false);
vars.out = arrayfun(@(i)sprintf('y%d',i),1:M,'UniformOutput',false);
dom = struct();
dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];
dom.out = randi([-2,2],[M,1]) +[zeros(M,1),rand([M,1])];
degs = struct('in',randi(2,[1,N]),'out',randi(2,[1,M]));
depmat = struct('in',eye(d,N)>0,'out',logical(randi([0,1],d,M)));
dnsty = 0.5;

Aop = rand_tensopvar(dims,vars,dom,depmat,randi(2),type,degs,dnsty);
order = randperm(d);
Bop = reorder_vars(Aop,order);

if ~isequal(Bop.ops,Aop.ops(:,order))
    error("Case 1: 'ops' not permuted as expected.")
elseif ~isequal(Bop.dims,Aop.dims(order,:))
    error("Case 1: 'dims' not permuted consistently with 'ops'.")
elseif ~isequal(Bop.depmat.in,Aop.depmat.in(order,:)) ||...
        ~isequal(Bop.depmat.out,Aop.depmat.out(order,:))
    error("Case 1: 'depmat' not permuted consistently with 'ops'.")
end

% Reordering by 'order' then by its inverse must recover the original object
[~,inv_order] = sort(order);
Cop = reorder_vars(Bop,inv_order);
if ~isequal(Cop.ops,Aop.ops) || ~isequal(Cop.dims,Aop.dims) || ~isequal(Cop.depmat,Aop.depmat)
    error("Case 1: reorder followed by its inverse did not recover the original object.")
end

% The identity permutation must be a no-op
Dop = reorder_vars(Aop,1:d);
if ~isequal(Dop.ops,Aop.ops) || ~isequal(Dop.dims,Aop.dims) || ~isequal(Dop.depmat,Aop.depmat)
    error("Case 1: identity permutation should leave the object unchanged.")
end

%%%% Case 2: single factor (d==1) -- must be a no-op regardless of matdim(2)
type1 = [1,2];
dims1 = [randi(3),randi(3)+1];  % column dim >1, to stress the old size(Cop,2)-based guard
vars1 = struct('in',{{'x1'}},'out',{{'y1'}});
dom1 = struct('in',[0,1],'out',[0,1]);
degs1 = struct('in',2,'out',2);
depmat1 = struct('in',true(1,1),'out',true(1,1));
Eop = rand_tensopvar(dims1,vars1,dom1,depmat1,1,type1,degs1,dnsty);
Fop = reorder_vars(Eop,1);
if ~isequal(Fop.ops,Eop.ops) || ~isequal(Fop.dims,Eop.dims)
    error("Case 2: reordering a single factor should be a no-op.")
end

%%%% Case 3: otype(2)==2 (Kronecker columns), multiple factors -- must error
% cleanly (not implemented), rather than crash on an undefined property.
type3 = [1,2];
dims3 = [repmat(randi(3),d,1), randi(3,[d,1])];  % uniform row dims (type3(1)==1); per-factor column dims may differ
Gop = rand_tensopvar(dims3,vars,dom,depmat,randi(2),type3,degs,dnsty);
threw_error = false;
try
    reorder_vars(Gop,order); %#ok<NASGU>
catch ME
    threw_error = true;
    if contains(ME.message,'Unrecognized property')
        error("Case 3: reorder_vars crashed on an undefined property instead of erroring cleanly.")
    end
end
if ~threw_error
    error("Case 3: reordering a multi-factor Kronecker-column operator should currently error.")
end

disp('Test_tensopvar_reorder_vars passed.')
