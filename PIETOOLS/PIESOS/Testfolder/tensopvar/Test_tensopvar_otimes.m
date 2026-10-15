% Test_tensopvar_otimes.m -- tests @tensopvar/otimes.m for single-factor
% operators A, B. The oracle is apply_op(otimes(A,B,type),{v1,v2})
% against Av.*Bv (type(1)==1) or kron(Av,Bv) (type(1)==2), with Av, Bv the
% actions of A and B on random inputs v1, v2. otype(1) alternates over
% {1,2}, otype(2)=2.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - For type(1)==2 the shared output variables of Av, Bv are renamed
%   <name>_1, <name>_2 before the Kronecker product, as apply_op names
%   them;
% - degs2 is initialized as a struct explicitly: a double 'degs2' left in
%   the workspace by another test (Test_intop_otimes,
%   Test_intop_dpvar_ops) would break 'degs2.in = ...' when run in the
%   same session.

ntrials = 20;
for trial=1:ntrials
% Set parameters defining the random operators
% Row and column dimensions
m1 = 2;     n1 = 2;      % row and column dimensions of functional
m2 = 2;     n2 = 3;
% Input and output variables
M1 = 1;     N1 = 1;     % number of input and output variables
M2 = 1;     N2 = 1;     % number of variables
% disjoint 'x'/'y' name pools: input and output variables must not share a name (rand_sopvar rejects one name with two domains)
vars1 = struct();
vars1.in = arrayfun(@(i)sprintf('x%d',i),1:N1,'UniformOutput',false);
vars1.out = arrayfun(@(i)sprintf('y%d',i),1:M1,'UniformOutput',false);
vars2 = struct();
vars2.in = arrayfun(@(i)sprintf('x%d',i),1:N2,'UniformOutput',false);
vars2.out = arrayfun(@(i)sprintf('y%d',i),1:M2,'UniformOutput',false);
% Domains of the variables
dom1 = struct();
dom1.in = randi([-2,2],[N1,1]) +[zeros(N1,1),rand([N1,1])];
dom1.out = randi([-2,2],[M1,1]) +[zeros(M1,1),rand([M1,1])];
dom2 = struct();
dom2.in = randi([-2,2],[N2,1]) +[zeros(N2,1),rand([N2,1])];
dom2.out = randi([-2,2],[M2,1]) +[zeros(M2,1),rand([M2,1])];
% Monomial degrees
degs1 = struct();
degs1.in = randi(3,[1,N1]);     degs1.out = randi(3,[1,M1]);
degs2 = struct();
degs2.in = randi(3,[1,N2]);     degs2.out = randi(3,[1,M2]);
% Density of the coefficient matrices
dnsty1 = 0.25;
dnsty2 = 0.5;

% Set the type of product to take
if mod(trial,2)==1
    type = [1,2];
else
    type = [2,2];    % type(2)==2 required
end
if type(1)==1
    m2 = m1;
end
if type(2)==1
    n2 = n1;
end

% Set parameters defining the function
p1 = 1;              % column dimension of function
p2 = 1;
pdegv = 2;

% Generate random PI operators
Aop = rand_sopvar([m1,n1],vars1,dom1,degs1,dnsty1);
Bop = rand_sopvar([m2,n2],vars2,dom2,degs2,dnsty2);
% Make sure the domain match
[~,ai,bi] = intersect(Aop.vars.in,Bop.vars.in,'stable');
Bop.dom.in(bi,:) = Aop.dom.in(ai,:);
[~,ai,bi] = intersect(Aop.vars.out,Bop.vars.out,'stable');
Bop.dom.out(bi,:) = Aop.dom.out(ai,:);

% Compute the product
Aop = tensopvar(Aop);
Bop = tensopvar(Bop);
Cop = otimes(Aop,Bop,type);

% Generate random functions on which the operators act
vars1 = Aop.vars;
vars2 = Bop.vars;
pvars1_in = polynomial(vars1.in(:));
pvars2_in = polynomial(vars2.in(:));
v1 = rand_poly([n1,p1],pvars1_in,pdegv,50);
v2 = rand_poly([n2,p2],pvars2_in,pdegv,50);
% Determine the variables of the output
if type(1)==2
    vars_out_shared = intersect(vars1.out,vars2.out,'stable');
    pvars_out = polynomial(vars_out_shared(:));
    vars1_out = vars_out_shared;
    vars2_out = vars_out_shared;
    for i=1:numel(vars_out_shared)
        vars1_out{i} = [vars1_out{i},'_',num2str(1)];
        vars2_out{i} = [vars2_out{i},'_',num2str(2)];
    end
    pvars1_out = polynomial(vars1_out(:));
    pvars2_out = polynomial(vars2_out(:));
end

% Apply the operators
Av = apply_op(Aop,v1);
Bv = apply_op(Bop,v2);
Cv = apply_op(Cop,{v1,v2});

% Compare
tol = 1e-8;
if all(type==[1,2])
    err_C = cleanpoly(Cv - Av.*Bv,1e-10);
elseif all(type==[2,2])
    Av1 = subs(Av,pvars_out,pvars1_out);
    Bv1 = subs(Bv,pvars_out,pvars2_out);
    err_C = cleanpoly(Cv - kron(Av1,Bv1),1e-10);
else
    error("Test-case not supported.")
end
disp('err_C = ')
disp(max(max(abs(err_C.C))))
if max(max(abs(err_C.C)))>tol
    error("Trial %d (type=[%d,%d]): outputs do not match: residual exceeds tolerance.",trial,type(1),type(2))
end
end

fprintf('Test_tensopvar_otimes passed (%d trials, otype(1) in {1,2} alternating).\n',ntrials)
