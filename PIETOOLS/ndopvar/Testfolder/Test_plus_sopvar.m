
% % Set parameters for declaring random operator
% Dimensions
m = 3;                  % row dimension of operator
n = 4;                  % column dimension of operator
p = 2;                  % column dimension of function
% Input and output variables
M = randi([0,3]);      % Number of input variables
N = randi([0,3]);      % Number of output variables
vars = struct();
vars.in = cell(1,N);
idcs = sort(randperm(N+M,N));
for i=1:N
    vars.in{i} = ['s',num2str(idcs(i))];
end
vars.out = cell(1,M);
idcs = sort(randperm(N+M,M));
for i=1:M
    vars.out{i} = ['s',num2str(idcs(i))];
end
% Domains of the variables
dom = struct();
dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];
dom.out = randi([-2,2],[M,1]) +[zeros(M,1),rand([M,1])];
% Monomial degrees
degs = struct();
degs.in = randi(3,[1,N]);     degs.out = randi(3,[1,M]);
pdeg_x = 4;             % maximal monomial degree of function 
% Density of the coefficient matrices
dnsty = 0.25;

% Generate random sopvar object
Aop = rand_sopvar([m,n],vars,dom,degs,dnsty);
Bop = rand_sopvar([m,n],vars,dom,degs,dnsty);

% Generate random function
vars = Aop.vars;
pvars_in = polynomial(vars.in(:));
x = rand_poly([n,p],pvars_in,pdeg_x,50);

% Compute the sum
Cop = Aop+Bop;
Ax = apply_sopvar(Aop,x);
Bx = apply_sopvar(Bop,x);
Cx = apply_sopvar(Cop,x);

err = cleanpoly(Cx-(Ax+Bx),1e-12);
disp("err = ")
disp(num2str(max(max(abs(err.C)))));


