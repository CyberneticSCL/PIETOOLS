
% Generate a random quadPoly object
m = 3;      % row dimension
n = 2;      % column dimension
M_new = 6;  % new number of output variables
N_new = 6;  % new number of input variables
M = randi([0,M_new]);      % number of output variables
N = randi([0,N_new]);      % number of input variables
s_idcs = sort(randperm(M_new,M));
t_idcs = sort(randperm(N_new,N));
pdeg1 = randi([0,3],N,1);
pdeg2 = randi([0,3],M,1);
nZ_max = 50;

% Declare the variables names
var1_new = cell(1,M_new);
for i=1:M_new
    var1_new{i} = ['s',num2str(i)];
end
var1 = var1_new(s_idcs);
var2_new = cell(1,N_new);
for i=1:N_new
    var2_new{i} = ['t',num2str(i)];
end
var2 = var2_new(t_idcs);

A_poly = rand_poly([m,n],[var1,var2]',[pdeg1;pdeg2],nZ_max);
A = quadPoly.polynomial2quadPoly(A_poly,var1,var2);

B = set_vars(A,var1_new,var2_new);

err = A-B

