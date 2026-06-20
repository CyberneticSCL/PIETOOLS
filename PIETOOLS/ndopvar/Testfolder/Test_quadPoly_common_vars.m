
% Generate a random quadPoly object
m1 = 3;      % row dimension
n1 = 2;      % column dimension
m2 = 2;      % row dimension
n2 = 3;      % column dimension
M = 6;      % max number of output variables
N = 6;      % max number of input variables
M1 = randi([0,M]);      % number of output variables
N1 = randi([0,N]);      % number of input variables
M2 = randi([0,M]);      % number of output variables
N2 = randi([0,N]);      % number of input variables
s_idcs1 = sort(randperm(M,M1));
t_idcs1 = sort(randperm(N,N1));
s_idcs2 = sort(randperm(M,M2));
t_idcs2 = sort(randperm(N,N2));
pdeg1 = randi([0,3],N1+M1,1);
pdeg2 = randi([0,3],N2+M2,1);
nZ_max = 50;

% Declare the variables names
var1_full = cell(1,M);
for i=1:M_new
    var1_full{i} = ['s',num2str(i)];
end
svar1 = var1_full(s_idcs1);
svar2 = var1_full(s_idcs2);
var2_full = cell(1,N);
for i=1:N_new
    var2_full{i} = ['t',num2str(i)];
end
tvar1 = var2_full(t_idcs1);
tvar2 = var2_full(t_idcs2);

A_poly = rand_poly([m1,n1],[svar1,tvar1]',pdeg1,nZ_max);
B_poly = rand_poly([m2,n2],[svar2,tvar2]',pdeg2,nZ_max);
A = quadPoly.polynomial2quadPoly(A_poly,svar1,tvar1);
B = quadPoly.polynomial2quadPoly(B_poly,svar2,tvar2);

[A2,B2] = common_vars(A,B);

err = A-A2
err2 = B-B2

