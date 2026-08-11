clc; clear;
% TEST_SCRIPT tests addition of sdopvar objects against ndopvar addition.
%
% The test constructs random ndopvar objects, converts them to
% sdopvar, performs sdopvar addition, converts the
% result back to ndopvar, and compares with ndopvar addition.

% Set the spatial dimension and size of the operators.
N = 2;
m = 2;
n = 2;
% Set the maximal monomial degree.
d = 2;
% Declare N variables and their domains.
var1_name = [repmat('s',[N,1]),num2str((1:N)')];
var2_name = [var1_name,repmat('_dum',[N,1])];
var1 = polynomial(mat2cell(var1_name,ones(N,1),size(var1_name,2)));
var2 = polynomial(mat2cell(var2_name,ones(N,1),size(var2_name,2)));
dom = [zeros(N,1),ones(N,1)];
% Declare decision variables.
ndec = 10;
dvarname = cell(ndec,1);
for idx = 1:ndec
    dvarname{idx} = ['dvar_',num2str(idx)];
end
% Generate random ndopvar objects.
A_ndop = rand_ndopvar([m,n],d,dom,var1,var2,dvarname);
B_ndop = rand_ndopvar([m,n],d,dom,var1,var2,dvarname);
% Convert to sdopvar
A_sdop = ndopvar2sdopvar(A_ndop);
B_sdop = ndopvar2sdopvar(B_ndop);
% Addition in sdopvar representation.
C_sdop = A_sdop + B_sdop;
C_from_sdop = sdopvar2ndopvar(C_sdop);

% Direct addition in ndopvar representation.
C_ndop = A_ndop + B_ndop;

% Compare the two ndopvar results.
tol = 1e-10;
if ~eq(C_from_sdop,C_ndop,tol)
    error('sdopvar addition test failed: converted sum differs from ndopvar sum.');
end

disp('sdopvar addition test passed.');
