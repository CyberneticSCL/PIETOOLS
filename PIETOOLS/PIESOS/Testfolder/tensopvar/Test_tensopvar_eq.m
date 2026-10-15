% Test_tensopvar_eq.m -- tests @tensopvar/eq.m: reflexivity, inequality
% of independent random operators, and inequality of operators with
% different factor dimensions. otype(1) alternates over {1,2}, otype(2)=2.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - 'Aop==(Aop-Aop+Aop)' and '(Aop-Aop)==0' are not tested: eq only
%   detects zero via an all-zero factor in each term, and plus/minus do
%   not cancel terms, so Aop-Aop (two terms with no zero factor) is
%   reported nonzero. The isequal(A,B) fast path only helps for
%   structurally identical objects.

ntrials = 20;
for trial=1:ntrials
nterms = randi(3);
d = randi(5);
N = randi(3);
M = randi(3);
if mod(trial,2)==1
    type = [1,2];
else
    type = [2,2];
end
dnsty = 0.5;

dims = randi(4,[d,2]);
if type(1)==1
    dims(:,1) = dims(1,1);
end
if type(2)==1
    dims(:,2) = dims(1,2);
end
vars = struct();
% disjoint 'x'/'y' name pools: input and output variables must not share a name (rand_sopvar rejects one name with two domains)
vars.in = arrayfun(@(i)sprintf('x%d',i),1:N,'UniformOutput',false);
vars.out = arrayfun(@(i)sprintf('y%d',i),1:M,'UniformOutput',false);
dom = struct();
dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];
dom.out = randi([-2,2],[M,1]) +[zeros(M,1),rand([M,1])];
degs = struct();
degs.in = randi(3,[1,N]);     degs.out = randi(3,[1,M]);
depmat = struct();
depmat.in = logical(randi([0,1],d,N));
depmat.out = logical(randi([0,1],d,M));

Aop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);
Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);

%%%% A should equal itself
if ~(Aop==Aop)
    error("Trial %d (type=[%d,%d]): Aop==Aop should be true.",trial,type(1),type(2))
end

%%%% A should not equal B (generically, for independently-random A,B)
if Aop==Bop
    % unlikely but possible: warn rather than fail
    warning("Aop==Bop returned true for independently-random operators; this is not itself a bug (could be a rare coincidence), but worth a second look if seen repeatedly.")
end

%%%% Structural mismatch should register as not-equal, not error
% bump the row dim of all factors: pointwise rows (type(1)==1) require uniform dims(:,1)
dims2 = dims; dims2(:,1) = dims2(:,1)+1;
Dop = rand_tensopvar(dims2,vars,dom,depmat,nterms,type,degs,dnsty);
if Aop==Dop
    error("Trial %d (type=[%d,%d]): operators with different factor dimensions should not be equal.",trial,type(1),type(2))
end
end

fprintf('Test_tensopvar_eq passed (%d trials, otype(1) in {1,2} alternating).\n',ntrials)
