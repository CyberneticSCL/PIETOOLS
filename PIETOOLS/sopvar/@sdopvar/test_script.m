clc; clear;
% TEST_SCRIPT tests sdopvar addition, transpose, and equality.
%
% The addition test constructs random ndopvar objects, converts them to
% sdopvar, performs sdopvar addition, converts the result back to ndopvar,
% and compares with ndopvar addition. The eq test checks equality after a
% monomial-basis expansion. The transpose test checks that double transpose 
% returns the orignal sdopvar.
rng(1);
tol = 1e-10;

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

% Convert to sdopvar.
A_sdop = ndopvar2sdopvar(A_ndop);
B_sdop = ndopvar2sdopvar(B_ndop);

%% Addition
C_sdop = A_sdop + B_sdop;
C_from_sdop = sdopvar2ndopvar(C_sdop);
C_ndop = A_ndop + B_ndop;
if ~eq(C_from_sdop,C_ndop,tol)
    error('sdopvar addition test failed: converted sum differs from ndopvar sum.');
end
disp('sdopvar addition test passed.');

%% Equality with different monomial bases
vars_eq = struct('in',{{'s'}},'out',{{'s'}});
dom_eq = struct('in',[0 1],'out',[0 1]);
params_eq_1 = struct('A',{{sparse([1;2])}},'B',{{sparse([3 4])}});
params_eq_2 = struct('A',{{sparse([1;2;0])}},'B',{{sparse([3 4 0])}});
P_eq = sdopvar(params_eq_1,vars_eq,{'d1'},{[0;1]},{0},dom_eq,[1 1]);
Q_eq = sdopvar(params_eq_2,vars_eq,{'d1'},{[0;1;2]},{0},dom_eq,[1 1]);
if ~eq(P_eq,Q_eq,tol)
    error('sdopvar eq test failed: equal operators with different monomial bases compare unequal.');
end
if eq(P_eq,2*Q_eq,tol)
    error('sdopvar eq test failed: unequal operators compare equal.');
end
disp('sdopvar eq test passed.');

%% Transpose
At_sdop = A_sdop';
Att_sdop = At_sdop';
Att_ndop = sdopvar2ndopvar(Att_sdop);
if ~eq(Att_ndop,A_ndop,tol)
    error('sdopvar transpose test failed: double transpose differs from original operator.');
end
if ~eq(Att_sdop,A_sdop,tol)
    error('sdopvar transpose test failed: double transpose sdopvar comparison failed.');
end
disp('sdopvar transpose test passed.');
