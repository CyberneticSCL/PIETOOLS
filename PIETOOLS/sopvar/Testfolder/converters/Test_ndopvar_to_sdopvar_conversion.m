clc;clear;
% d = 5;
% n = 3;
% m = 4;
% D = rand(1, d);
% C = rand(d*n, m);
% AAA =  kron(eye(n), D)*C;   % first representation
% BBB = D*reshape(C, d, n*m); % second representation
% max(abs(reshape(AAA, 1, [])- BBB))


% Set the spatial dimension and size of the operators
N = 2;
m = 2;
k = 2;
n = 2;

% Set the maximal monomial degree
d = 2;
% Declare N variables and their domain
var1_name = [repmat('s',[N,1]),num2str((1:N)')];
var2_name = [var1_name,repmat('_dum',[N,1])];
var1 = polynomial(mat2cell(var1_name,ones(N,1),size(var1_name,2)));
var2 = polynomial(mat2cell(var2_name,ones(N,1),size(var2_name,2)));
dom = [zeros(N,1),ones(N,1)];

% Declare the decision variables
dvarname = cell(100,1);
for idx = 1:length(dvarname)
    dvarname{idx} = ['dvar_', num2str(idx)];
end
% Declare random ndopvar objects
Qop = rand_ndopvar([m,k],d,dom,var1,var2,dvarname);
Rop1 = rand_ndopvar([m,k],d,dom,var1,var2,cell(0, 1));
Rop2 = rand_ndopvar([m,k],d,dom,var1,var2,cell(0, 1));

% ndopvar to sdopvar conversion
sdopvar_obj = ndopvar2sdopvar(Qop);
% sdopvar to ndopvar conversion
Qop_2 = sdopvar2ndopvar(sdopvar_obj);

Qop == Qop_2 

% verify product using nopvar vs sopvar
% nopvar to sopvar conversion
Rop_sopvar1 = nopvar2sopvar(Rop1);
Rop_sopvar2 = nopvar2sopvar(Rop2);
prod1 = Rop1*Rop2;
prod2_sopvar = Rop_sopvar1*Rop_sopvar2;
% sopvar to nopvar conversion
prod2 = sopvar2nopvar(prod2_sopvar);
prod1 == prod2