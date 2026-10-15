% Test_tensopmat_initialize.m -- tests the 'tensopmat' constructor for a
% single 'tensopvar' argument (@tensopmat/tensopmat.m, the "elseif
% isa(Cop,'tensopvar')||isa(Cop,'intop')" branch), which must collapse
% depmat via the otype-aware 'Cop.deps' (max for a pointwise row/column,
% sum for a Kronecker one), not 'sum(Cop.depmat,1)'.
% Oracle: the expected depmat is computed BY HAND from the construction
% (not from 'block.deps', which would test the constructor against
% itself); apply_op is used only to check that wrapping in 'tensopmat'
% does not corrupt the computation.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - rand_tensopmat does not exercise this path: it builds from a p x d
%   cell of blocks via @tensopmat/initialize.m;
% - The difference between sum and max is only observable for a
%   MULTI-FACTOR 'tensopvar' (depmat has one row per factor; sum and max
%   coincide for one factor), so the block is built directly via
%   rand_tensopvar with both factors depending on the SAME variable,
%   making sum (=2) and max (=1) differ.

tol = 1e-8;
nreps = 20;
dom = [0,1];

for rep_idx=1:nreps

d_out = randi([1,2]);              % common row dim for both factors (needed for otype1==1 validity)
d_in1 = randi([1,2]);   d_in2 = randi([1,2]);
dims = [d_out,d_in1; d_out,d_in2]; % 2 factors, same output dim
vars = struct('in',{{'s1'}},'out',{{'t1'}});   % BOTH factors share the same sole variable
domst = struct('in',dom,'out',dom);
depmat = struct('in',true(2,1),'out',true(2,1));
K = randi([1,2]);
otype1 = randi(2);
block = rand_tensopvar(dims,vars,domst,depmat,K,[otype1,2],randi([1,3]),[]);

if size(block.depmat.out,1)~=2 || size(block.depmat.in,1)~=2
    error("Test_tensopmat_initialize failed: expected a 2-factor (2-row depmat) block, got %s/%s.", ...
        mat2str(size(block.depmat.out)),mat2str(size(block.depmat.in)))
end

Cop = tensopmat(block);

% Expected: max(1,1)=1 (otype1==1), sum(1,1)=2 (otype1==2); otype(2)==2,
% so depmat.in is always summed
if otype1==1
    expected_out = 1;
else
    expected_out = 2;
end
expected_in = 2;

if ~isequal(Cop.depmat.out,expected_out)
    error("Test_tensopmat_initialize failed: Cop.depmat.out=%s, expected %d for otype(1)=%d.", ...
        mat2str(Cop.depmat.out),expected_out,otype1)
end
if ~isequal(Cop.depmat.in,expected_in)
    error("Test_tensopmat_initialize failed: Cop.depmat.in=%s, expected %d.", mat2str(Cop.depmat.in),expected_in)
end
% Secondary: constructor agrees with the public 'deps' getter
if ~isequal(Cop.depmat.in,block.deps.in) || ~isequal(Cop.depmat.out,block.deps.out)
    error("Test_tensopmat_initialize failed: Cop.depmat does not match block.deps.")
end

[logval,errmsg] = isvalid(Cop);
if ~logval
    error("Test_tensopmat_initialize failed: isvalid(Cop) failed: %s", errmsg)
end

% Wrapping in 'tensopmat' must not change the output of the bare 'tensopvar'
X1 = rand_poly([d_in1,1],{'s1'},1,3);
X2 = rand_poly([d_in2,1],{'s1'},1,3);
Y1 = apply_op(Cop,{{X1,X2}});
Y2 = apply_op(block,{X1,X2});
err = cleanpoly(Y1-Y2,1e-10);
if max(max(abs(full(err.C))))>tol
    error("Test_tensopmat_initialize failed: apply_op(Cop,.) does not match apply_op(block,.) directly.")
end

end

disp('Test_tensopmat_initialize passed.')
