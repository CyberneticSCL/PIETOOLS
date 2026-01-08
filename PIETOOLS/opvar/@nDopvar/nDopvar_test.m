clc; clear;
n = 6;

dim = [3,2];
vars_in=[];
for i=1:n
    eval("pvar s"+num2str(i));
    vars_in = [vars_in, eval("s"+num2str(i))];
end
vars_out = vars_in;
%%


mons = monomials(vars_in,0:5);

keys = 1:3^numel(vars_in);
params = cell(1,length(keys));
for i=1:length(keys)
    coef = rand(size(mons,1),prod(dim));
    params{i} = polynomial(coef,mons.degmat,mons.varname, dim);
end
paramsDict = dictionary(keys,params);
P = nDopvar(vars_in, vars_out, dim, paramsDict);


% addition 
tic;
Q = P+P;
tadd=toc;

% transpose
tic;
R = P';
tadj = toc;

% multiply
tic;
S = P*P;
tmul = toc;
