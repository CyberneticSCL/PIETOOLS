
% MMP, 09/17/2026: Draw one domain per variable NAME over the union of
% vars.in and vars.out and index dom.in/dom.out out of that table, replacing
% two independent per-side draws. A variable shared by the input and output
% space is one variable carrying one interval, so independent draws gave it
% two and 'rand_sopvar' rejected the operator (see commit e6842c92).

% % Set parameters for declaring random operator
% Dimensions
m = 3;                  % row dimension of operator
n = 4;                  % column dimension of operator
p = 2;                  % column dimension of function
nr = randi([0,5]);      % Number of repetitions along row dimension
nc = randi([0,5]);      % Number of repetitions along row dimension
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
% Domains of the variables. A variable in both vars.in and vars.out is one  % MMP, 09/17/2026
% variable carrying ONE interval, and 'rand_sopvar' rejects a clash, so the % MMP, 09/17/2026
% table is drawn over the union of names and indexed by name per side.      % MMP, 09/17/2026
% Independent per-side draws gave a shared variable two different domains.  % MMP, 09/17/2026
varsall = unique([vars.in,vars.out]);  nvarsall = numel(varsall);           % MMP, 09/17/2026
domall = randi([-2,2],[nvarsall,1]) +[zeros(nvarsall,1),rand([nvarsall,1])];% MMP, 09/17/2026
[~,i_in ] = ismember(vars.in, varsall);                                     % MMP, 09/17/2026
[~,i_out] = ismember(vars.out,varsall);                                     % MMP, 09/17/2026
dom = struct();
%   dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];                 % MMP, 09/17/2026 (was)
%   dom.out = randi([-2,2],[M,1]) +[zeros(M,1),rand([M,1])];                % MMP, 09/17/2026 (was)
dom.in  = domall(i_in,:);                                                   % MMP, 09/17/2026
dom.out = domall(i_out,:);                                                  % MMP, 09/17/2026
% Monomial degrees
degs = struct();
degs.in = randi(3,[1,N]);     degs.out = randi(3,[1,M]);
pdeg_x = 4;             % maximal monomial degree of function 
% Density of the coefficient matrices
dnsty = 0.25;

% Generate random sopvar object
Aop = rand_sopvar([m,n],vars,dom,degs,dnsty);

% Generate random function
vars = Aop.vars;
pvars_in = polynomial(vars.in(:));
x = rand_poly([n,p],pvars_in,pdeg_x,50);

% Perform the repmat
Mop = repmat(Aop,nr,nc);
x_rep = repmat(x,nc,1);

% Compute the sum
Ax = apply_sopvar(Aop,x);
Mx_true = nc*repmat(Ax,nr,1);
Mx_alt = apply_sopvar(Mop,x_rep);

err = cleanpoly(Mx_true-Mx_alt,1e-12);
disp("err = ")
disp(num2str(max(max(abs(err.C)))));


