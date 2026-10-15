% Test_tensopvar_cat.m -- tests @tensopvar/blkdiag.m, vertcat.m and
% horzcat.m.
%
% STRUCTURAL checks only (matdim of the resulting 'tensopmat'); the oracle
% is the expected size computed from the random block sizes.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - d=1 factor per operator: with Kronecker columns (type(2)==2) the
%   column dimension of a d-factor operator is the product of the
%   per-factor column dims.

k = 3;
type = [1,2];
dnsty = 0.5;
d = 1;

%%%% blkdiag: no matching requirements between pieces
ms = randi(3,[1,k]);   ns = randi(3,[1,k]);
Aops = cell(1,k);
for b=1:k
    N = randi(2);   M = randi(2);
    dims = randi(3,[d,2]);
    dims(:,1) = ms(b);  dims(:,2) = ns(b);
    vars = struct();
    vars.in = arrayfun(@(v) sprintf('b%d_x%d',b,v),1:N,'UniformOutput',false);
    vars.out = arrayfun(@(v) sprintf('b%d_y%d',b,v),1:M,'UniformOutput',false);
    dom = struct();
    dom.in = randi([-2,2],[N,1])+[zeros(N,1),rand([N,1])];
    dom.out = randi([-2,2],[M,1])+[zeros(M,1),rand([M,1])];
    degs = struct('in',randi(2,[1,N]),'out',randi(2,[1,M]));
    depmat = struct('in',logical(randi([0,1],d,N)),'out',logical(randi([0,1],d,M)));
    Aops{b} = rand_tensopvar(dims,vars,dom,depmat,randi(2),type,degs,dnsty);
end
Pop = blkdiag(Aops{:});
if ~isa(Pop,'tensopmat')
    error("blkdiag: output should be of type 'tensopmat'.")
elseif ~isequal(size(Pop),[sum(ms),sum(ns)])
    error("blkdiag: matrix dimensions of output do not match expected [sum(ms),sum(ns)].")
end

%%%% vertcat: pieces must share the same input space
% input space includes depmat.in: generated once and shared (vertcat requires vars.in, dims.in, depmat.in to match)
N = 2;
vars_shared = struct();
vars_shared.in = arrayfun(@(v) sprintf('vx%d',v),1:N,'UniformOutput',false);
dom_shared_in = randi([-2,2],[N,1])+[zeros(N,1),rand([N,1])];
degs_in_shared = randi(2,[1,N]);
depmat_in_shared = logical(randi([0,1],d,N));
ms = randi(3,[1,k]);
n = randi(3);
Bops = cell(1,k);
for b=1:k
    M = randi(2);
    dims = randi(3,[d,2]);
    dims(:,1) = ms(b);  dims(:,2) = n;
    vars = vars_shared;
    vars.out = arrayfun(@(v) sprintf('vb%d_y%d',b,v),1:M,'UniformOutput',false);
    dom = struct();
    dom.in = dom_shared_in;
    dom.out = randi([-2,2],[M,1])+[zeros(M,1),rand([M,1])];
    degs = struct('in',degs_in_shared,'out',randi(2,[1,M]));
    depmat = struct('in',depmat_in_shared,'out',logical(randi([0,1],d,M)));
    Bops{b} = rand_tensopvar(dims,vars,dom,depmat,randi(2),type,degs,dnsty);
end
Pop = vertcat(Bops{:});
if ~isa(Pop,'tensopmat')
    error("vertcat: output should be of type 'tensopmat'.")
elseif ~isequal(size(Pop),[sum(ms),n])
    error("vertcat: matrix dimensions of output do not match expected [sum(ms),n].")
end

%%%% horzcat: pieces must share the same output row-count (m)
% output space includes dims.out, vars.out and depmat.out: generated once and shared, as for vertcat
m = randi(3);
M = randi(2);
vars_out_shared = arrayfun(@(v) sprintf('vy%d',v),1:M,'UniformOutput',false);
dom_shared_out = randi([-2,2],[M,1])+[zeros(M,1),rand([M,1])];
degs_out_shared = randi(2,[1,M]);
depmat_out_shared = logical(randi([0,1],d,M));
ns = randi(3,[1,k]);
Cops = cell(1,k);
for b=1:k
    N = randi(2);
    dims = randi(3,[d,2]);
    dims(:,1) = m;  dims(:,2) = ns(b);
    vars = struct();
    vars.in = arrayfun(@(v) sprintf('c%d_x%d',b,v),1:N,'UniformOutput',false);
    vars.out = vars_out_shared;
    dom = struct();
    dom.in = randi([-2,2],[N,1])+[zeros(N,1),rand([N,1])];
    dom.out = dom_shared_out;
    degs = struct('in',randi(2,[1,N]),'out',degs_out_shared);
    depmat = struct('in',logical(randi([0,1],d,N)),'out',depmat_out_shared);
    Cops{b} = rand_tensopvar(dims,vars,dom,depmat,randi(2),type,degs,dnsty);
end
Pop = horzcat(Cops{:});
if ~isa(Pop,'tensopmat')
    error("horzcat: output should be of type 'tensopmat'.")
elseif ~isequal(size(Pop),[m,sum(ns)])
    error("horzcat: matrix dimensions of output do not match expected [m,sum(ns)].")
end

disp('Test_tensopvar_cat passed: blkdiag, vertcat and horzcat all produced correctly-sized tensopmat objects.')
