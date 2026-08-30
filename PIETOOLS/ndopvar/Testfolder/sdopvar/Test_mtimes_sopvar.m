clc; clear;
% Aop: L2^k[varmid] -> L2^m[varout]
% Bop: L2^n[varin] -> L2^k[varmid]

% runtime 
% 3D->3D opvars density(0.01) = 100-300 sec, density(0.9) = 10-60 sec;
% 2D->2D opvars density(0.01) = 0.05-0.16 sec, density(0.9) =  0.07-0.2 sec;
% 1D->1D opvars density(0.01) = 0.0025-0.0040 sec, density(0.9) =  0.0030-0.0090 sec;

% higher m x n scaling, degree scaling

m = 2; n = 2; k = 2;
maxvars = 1;

nvarsin = randi([1,floor(maxvars)]);
nvarsmid = randi([1,floor(maxvars)]);
nvarsout = randi([1,floor(maxvars)]);

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

Adegs = struct('in',randi(2,1,length(Avars.in)), ...
               'out',randi(2,1,length(Avars.out)));
Bdegs = struct('in',randi(2,1,length(Bvars.in)), ...
               'out',randi(2,1,length(Bvars.out)));

pdeg_x = 4;             % maximal monomial degree of test function
dnsty = 0.01;           % density of coefficient matrices

% Generate random sopvar operators
Bop = rand_sopvar([k,n],Bvars,Bdom,Bdegs,dnsty);
Aop = rand_sopvar([m,k],Avars,Adom,Adegs,dnsty);

% Generate random polynomial function compatible with Bop
vars_x = polynomial(Bop.vars.in(:));
x = rand_poly([n,1],vars_x,pdeg_x,50);
%%
% Compute composition

Cop = @() mtimes(Aop,Bop);
t=timeit(Cop);
Cop = Cop();

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