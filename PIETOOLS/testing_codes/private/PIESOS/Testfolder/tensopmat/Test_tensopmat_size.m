% Test_tensopmat_size.m -- tests @tensopmat/size.m.
%
% Checks all calling conventions: single 1x2 output, two separate outputs,
% and a single dimension queried by 'dim', against the ground truth
% matdim = [sum(dims.out),sum(dims.in)] built independently from the
% 'dims' field (not from size.m itself).

% DJ (CC), 10/15/2026: Initial coding;

p = randi([1,4]);  d = randi([1,4]);
N = randi([2,4]);  M = randi([2,4]);
dom = [0,1];

[Cop,~] = rand_tensopmat(p,d,N,M,dom,0.7);
true_dim = [sum(Cop.dims.out),sum(Cop.dims.in)];

sz = size(Cop);
if ~isequal(sz,true_dim)
    error("Test_tensopmat_size failed: size(Cop) does not match [sum(dims.out),sum(dims.in)].")
end

[nr,nc] = size(Cop);
if nr~=true_dim(1) || nc~=true_dim(2)
    error("Test_tensopmat_size failed: [nr,nc]=size(Cop) does not match true dimensions.")
end

if size(Cop,1)~=true_dim(1)
    error("Test_tensopmat_size failed: size(Cop,1) does not match number of rows.")
end
if size(Cop,2)~=true_dim(2)
    error("Test_tensopmat_size failed: size(Cop,2) does not match number of columns.")
end

% Invalid dimension must error
did_error = false;
try
    size(Cop,3); %#ok<NASGU>
catch
    did_error = true;
end
if ~did_error
    error("Test_tensopmat_size failed: size(Cop,3) should have errored.")
end

disp('Test_tensopmat_size passed.')
