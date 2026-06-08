
% Declare the independent variables
pvar s1 s1_dum
var1 = s1;
var2 = s1_dum;
% Declare the domain
dom = [0,2];
% Declare the degree of the parameters in the independent variables
deg1 = randi(4)-1;
deg2 = randi(4)-1;
% Declare the input and output dimensions of the operators
m = randi(5);      % row dimension of Pop
n = randi(5);      % column dimension of Pop
q = randi(5);      % column dimension of x
dim_idx1 = randi([1,2]);
dim_idx2 = randi([1,2]);
dim = zeros(2,2);
dim(dim_idx1,1) = m;
dim(dim_idx2,2) = n;

% Generate a random opvar object of the desired dimensions
Pop = rand_opvar(dim,deg1,var1,var2,dom);

% Convert to 'sopvar' object
Psop = opvar2sopvar(Pop);

% Generate a random polynomial
if dim_idx1==1
    x = rand([n,q]);
else
    x = rand_poly([n,q],var1,3,10);
end

% Apply the opvar object
Popx = apply_opvar(Pop,x);
Psopx = apply_sopvar(Psop,x);

err = cleanpoly(polynomial(Popx - Psopx),1e-10);
disp(['nrm = ',num2str(max(max(abs(err.C))))]);