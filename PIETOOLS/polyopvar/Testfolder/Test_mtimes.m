% Test mtimes.m

m = 3; n = 2; deg = 2;

% Instantiate v.
v = polyopvar();
v.varname = 'v';
v.degmat = 1;
v.pvarname = {'s1';'s2'};
v.dom = [0,1; -1,1];
v.varmat = [1, 1];

% Instantiate T.
N = ndims(v.pvarname);
dim = [m n];
dom = v.dom;
var1_name = [repmat('s',[N,1]),num2str((1:N)')];
var2_name = [var1_name,repmat('_dum',[N,1])];
var1 = polynomial(mat2cell(var1_name,ones(N,1),size(var1_name,2)));
var2 = polynomial(mat2cell(var2_name,ones(N,1),size(var2_name,2))); % dummy vars
dvarname = cell(0,1); % no decision vars.
Top = rand_ndopvar(dim,deg,dom,var1,var2,dvarname);


p = mtimes(Top,v);