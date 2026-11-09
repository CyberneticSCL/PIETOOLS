% Test_intop_dims.m -- tests the field 'dims' of 'intop' objects.
%
% For random functionals with known matrix dimensions [m,n], with and
% without decision variables, K.dims is a struct with K.dims.in=n and
% K.dims.out=m, as for 'tensopmat' objects. The values are compared with the
% dimensions passed to rand_intop (not with K.matdim), after indexing with
% '()', and for a functional without spatial variables. The field cannot be
% set, and is not changed by assigning the other fields.
%
% DJ (CC), 10/15/2026: Initial coding;

ntrials = 30;
for trial=1:ntrials
    m = randi(3);       n = randi(3);
    N = randi(2);       q = randi([0,3]);
    Kop = rand_intop([m,n],N,[0,1],1,1,q);

    % Field names and values
    D = Kop.dims;
    if ~isstruct(D) || ~isequal(sort(fieldnames(D)),{'in';'out'})
        error("Test_intop_dims failed: dims must be a struct with fields 'in' and 'out' (trial=%d).",trial)
    elseif ~isequal(D.in,n) || ~isequal(D.out,m)
        error("Test_intop_dims failed: dims.in=%d, dims.out=%d, expected %d and %d (trial=%d).",D.in,D.out,n,m,trial)
    end
    % Chained subsref
    if ~isequal(Kop.dims.in,n) || ~isequal(Kop.dims.out,m)
        error("Test_intop_dims failed: chained subsref of dims (trial=%d).",trial)
    end

    % Indexing
    ridcs = unique(randi(m,[1,2]));      cidcs = unique(randi(n,[1,3]));
    Ks = Kop(ridcs,cidcs);
    if ~isequal(Ks.dims.out,numel(ridcs)) || ~isequal(Ks.dims.in,numel(cidcs))
        error("Test_intop_dims failed: dims after indexing (trial=%d).",trial)
    end
end

% Functional without spatial variables
c = randn(2,3);
K0 = intop(c);
if ~isequal(K0.dims.out,2) || ~isequal(K0.dims.in,3)
    error("Test_intop_dims failed: functional without spatial variables.")
end

% Read-only
Kop = rand_intop([2,2],1,[0,1],1,1);
try
    Kop.dims = struct('in',1,'out',1);
    error("Test_intop_dims failed: dims could be set.")
catch ME
    if contains(ME.message,"Test_intop_dims failed")
        rethrow(ME)
    end
end

disp("Test_intop_dims passed.")
