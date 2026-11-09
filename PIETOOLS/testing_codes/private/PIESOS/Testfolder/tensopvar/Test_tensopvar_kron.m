% Test_tensopvar_kron.m -- tests @tensopvar/kron.m, repmat.m and
% repelem.m.
%
% STRUCTURAL checks only (matdim of the result, against the expected size
% from the sizes of A and B); a functional check is in
% Test_tensopvar_class (repmat, repelem). otype(1) alternates over {1,2},
% otype(2)=2.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - d>=2 factors: with a single factor otype(1) cannot affect matdim.

ntrials = 20;
for trial=1:ntrials
nterms = randi(2);
d = randi(2)+1;
N = randi(3);
M = randi(2);
if mod(trial,2)==1
    type = [1,2];
else
    type = [2,2];
end
dnsty = 0.5;

dims = randi(3,[d,2]);
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
degs = struct('in',randi(2,[1,N]),'out',randi(2,[1,M]));
depmat = struct('in',logical(randi([0,1],d,N)),'out',logical(randi([0,1],d,M)));

Aop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);
[mA,nA] = size(Aop);

%%%% repmat
r1 = randi(3); r2 = randi(3);
Brm = repmat(Aop,r1,r2);
if ~isequal(size(Brm),[mA*r1,nA*r2])
    error("Trial %d (type=[%d,%d]): repmat: matrix dimensions of output do not match expected [mA*r1,nA*r2].",trial,type(1),type(2))
end

%%%% repelem
Bre = repelem(Aop,r1,r2);
if ~isequal(size(Bre),[mA*r1,nA*r2])
    error("Trial %d (type=[%d,%d]): repelem: matrix dimensions of output do not match expected [mA*r1,nA*r2].",trial,type(1),type(2))
end

%%%% kron
Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);
[mB,nB] = size(Bop);
Cop = kron(Aop,Bop);
if ~isa(Cop,'tensopvar')
    error("Trial %d (type=[%d,%d]): kron: output should be of type 'tensopvar'.",trial,type(1),type(2))
elseif ~isequal(size(Cop),[mA*mB,nA*nB])
    error("Trial %d (type=[%d,%d]): kron: matrix dimensions of output do not match expected [mA*mB,nA*nB].",trial,type(1),type(2))
end
end

fprintf('Test_tensopvar_kron passed (%d trials, otype(1) in {1,2} alternating).\n',ntrials)
