% Test_tensopvar_isvalid.m -- tests @tensopvar/isvalid.m directly: a
% properly initialized object is reported valid, and a handful of
% corrupted copies (non-operator entry in ops, duplicate variable name,
% empty domain, mismatched row dimensions for otype(1)==1) are reported
% invalid without error. otype(1) alternates over {1,2}.
%
% DJ (CC), 10/15/2026: Initial coding;

ntrials = 20;
for trial=1:ntrials
nterms = randi(3);
d = randi(3)+1;
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
vars = struct();
vars.in = arrayfun(@(i)sprintf('x%d',i),1:N,'UniformOutput',false);
vars.out = arrayfun(@(i)sprintf('y%d',i),1:M,'UniformOutput',false);
dom = struct();
dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];
dom.out = randi([-2,2],[M,1]) +[zeros(M,1),rand([M,1])];
degs = struct('in',randi(3,[1,N]),'out',randi(3,[1,M]));
depmat = struct('in',logical(randi([0,1],d,N)),'out',logical(randi([0,1],d,M)));
depmat.in(1,~any(depmat.in,1)) = true;
depmat.out(1,~any(depmat.out,1)) = true;

Aop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);

%%%% A properly-initialized object must be reported valid, without error
[logval,errmsg] = isvalid(Aop);
if ~logval
    error("Trial %d (type=[%d,%d]): a properly-constructed tensopvar object was reported invalid: %s",trial,type(1),type(2),errmsg)
end

%%%% Corrupting 'ops' with a non-cell/non-sopvar entry must be caught
Bop = Aop;
Bop.ops{1,1} = 'not an operator';
[logval,~] = isvalid(Bop);
if logval
    error("Trial %d (type=[%d,%d]): an object with a non-operator entry in 'ops' should be reported invalid.",trial,type(1),type(2))
end

%%%% Corrupting 'vars.in' with a duplicate name must be caught
if N>=2
    Cop = Aop;
    Cop.vars.in{2} = Cop.vars.in{1};
    [logval,~] = isvalid(Cop);
    if logval
        error("Trial %d (type=[%d,%d]): an object with duplicate input variable names should be reported invalid.",trial,type(1),type(2))
    end
end

%%%% Corrupting 'dom' (upper bound below lower bound) must be caught
Dop = Aop;
Dop.dom.in(1,:) = [1,0];
[logval,~] = isvalid(Dop);
if logval
    error("Trial %d (type=[%d,%d]): an object with an invalid (empty) domain should be reported invalid.",trial,type(1),type(2))
end

%%%% Corrupting 'otype' to violate the elementwise-product dimension requirement must be caught
% (type(1)==1 && max(dims(:,1))~=min(dims(:,1)))
if d>=2 && max(Aop.dims(:,1))==min(Aop.dims(:,1))
    Eop = Aop;
    Eop.otype(1) = 1;
    Eop.dims(1,1) = Eop.dims(1,1)+1;
    [logval,~] = isvalid(Eop);
    if logval
        error("Trial %d (type=[%d,%d]): an otype(1)==1 object with mismatched row dimensions should be reported invalid.",trial,type(1),type(2))
    end
end
end

fprintf('Test_tensopvar_isvalid passed (%d trials, otype(1) in {1,2} alternating).\n',ntrials)
