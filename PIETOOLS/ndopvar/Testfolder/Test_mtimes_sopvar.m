clc; clear;
% Aop: L2^k[varmid] -> L2^m[varout]
% Bop: L2^n[varin] -> L2^k[varmid]

m = 3; n = 2; k = 1;
maxvars = 6;

nvarsin = randi([1,floor(maxvars/2)]);
nvarsmid = randi([1,floor(maxvars/2)]);
nvarsout = randi([1,floor(maxvars/2)]);

varin = cellstr(sort("s"+string(randperm(maxvars,nvarsin))));
varmid = cellstr(sort("s"+string(randperm(maxvars,nvarsmid))));
varout = cellstr(sort("s"+string(randperm(maxvars,nvarsout))));
varsall = unique([varin,varmid,varout]);
nvarsall = length(varsall);

domall = randn(nvarsall,1)+abs(randn(nvarsall,1)).*[0,1];

Avars = struct('in',{varmid}, 'out',{varout});
Bvars = struct('in',{varin}, 'out',{varmid});
Adom = struct( ...
        'in', domall(ismember(varsall,varmid),:), ...
        'out', domall(ismember(varsall,varout),:));
Bdom = struct( ...
        'in', domall(ismember(varsall,varin),:), ...
        'out', domall(ismember(varsall,varmid),:));

Adegs = struct('in',randi(3,1,length(Avars.in)), ...
               'out',randi(3,1,length(Avars.out)));
Bdegs = struct('in',randi(3,1,length(Bvars.in)), ...
               'out',randi(3,1,length(Bvars.out)));

pdeg_x = 4;             % maximal monomial degree of test function
dnsty = 0.1;           % density of coefficient matrices

% Generate random sopvar operators
Bop = rand_sopvar([k,n],Bvars,Bdom,Bdegs,dnsty);
Aop = rand_sopvar([m,k],Avars,Adom,Adegs,dnsty);

% Generate random polynomial function compatible with Bop
vars_x = polynomial(Bop.vars.in(:));
x = rand_poly([n,1],vars_x,pdeg_x,50);
%%
% Compute composition
Cop = mtimes3(Aop,Bop);

% Compare direct composed application with sequential application
Bx = apply_sopvar(Bop,x);
ABx = apply_sopvar(Aop,Bx);
Cx = apply_sopvar(Cop,x);

err = cleanpoly(Cx - ABx,1e-12);

if isempty(err.C)
    errmax = 0;
else
    errmax = max(max(abs(err.C)));
end

disp("composition err = ")
disp(num2str(errmax));