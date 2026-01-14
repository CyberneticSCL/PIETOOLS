
% Set the spatial dimension and size of the operators
N = 3;
m = 1;
k = 1;
n = 1;

% Set the maximal monomial degree
d = 3;

% Declare N variables and their domain
var1_name = [repmat('s',[N,1]),num2str((1:N)')];
var2_name = [var1_name,repmat('_dum',[N,1])];
var1 = polynomial(mat2cell(var1_name,ones(N,1),size(var1_name,2)));
var2 = polynomial(mat2cell(var2_name,ones(N,1),size(var2_name,2)));
dom = [zeros(N,1),ones(N,1)];

% Declare the decision variables
dvarname = cell(0,1);

% Declare random ndopvar objects
Qop = rand_ndopvar([m,k],d,dom,var1,var2,dvarname);
Rop = rand_ndopvar([k,n],d,dom,var1,var2,{});

% Compute the composition
Pop = Qop*Rop;