% Test_tensopvar_otimes.m -- tests @tensopvar/otimes.m for single-factor
% operators A, B. The oracle is apply_op(otimes(A,B,type),{v1,v2})
% against Av.*Bv (type(1)==1) or kron(Av,Bv) (type(1)==2), with Av, Bv the
% actions of A and B on random inputs v1, v2. otype(1) alternates over
% {1,2}, otype(2)=2.
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/22/2026: Add products of different types: concatenated for
%   the same type, nested for different types, a row of depmat for each
%   factor of the layer;

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

%% Products of operators of different types
% Structure: a product of operators of the same type concatenates the factors (a row of
% depmat for each); a product of operators of different types, e.g. the pointwise product
% A1.*A2 with B1 in a tensor product, nests them as two factors, with a row of depmat for
% each factor of the layer (the dependencies of the nested operator as a whole), whether or
% not the factors have factors themselves. Values: apply_op of the product on one input per
% leaf factor, against the outputs of each leaf factor renamed to its own copy y_1, y_2, y_3
% of the output variable y, in order, and multiplied
rng(31);
nvars = struct('in',{{'x'}},'out',{{'y'}});
ndom = struct('in',[0,1],'out',[0,1]);
ndegs = struct('in',2,'out',2);
A1 = tensopvar(rand_sopvar([1,1],nvars,ndom,ndegs,0.6));
A2 = tensopvar(rand_sopvar([1,1],nvars,ndom,ndegs,0.6));
B1 = tensopvar(rand_sopvar([1,1],nvars,ndom,ndegs,0.6));
px = polynomial({'x'});
v1 = rand_poly([1,1],px,2,50);      v2 = rand_poly([1,1],px,2,50);      v3 = rand_poly([1,1],px,2,50);
a1 = apply_op(A1,v1);       a2 = apply_op(A2,v2);       b1 = apply_op(B1,v3);
pw = a1.*a2;                % (A1 v1)(y).*(A2 v2)(y)
Apw = otimes(A1,A2,[1,2]);  % pointwise: 2 factors, otype [1,2]
N = otimes(A1,B1,[2,2]);    % tensor product: 2 factors, otype [2,2]
% same type: concatenated, not nested
check_struct(N,[1,2],[0,0],{A1.deps.out;B1.deps.out},'tensor product of two single factors');
N3 = otimes(N,A2,[2,2]);
check_struct(N3,[1,3],[0,0,0],{A1.deps.out;B1.deps.out;A2.deps.out},'tensor product of a tensor product and a single factor');
check_struct(Apw,[1,2],[0,0],{A1.deps.out;A2.deps.out},'pointwise product of two single factors');
% different types: nested, a row of depmat for each factor of the layer
Ca = otimes(Apw,B1,[2,2]);          % pointwise factor, single factor
check_struct(Ca,[1,2],[1,1],{Apw.deps.out;B1.deps.out},'pointwise factor o single factor');
Cb = otimes(B1,Apw,[2,2]);          % single factor, pointwise factor
check_struct(Cb,[1,2],[1,1],{B1.deps.out;Apw.deps.out},'single factor o pointwise factor');
Cc = otimes(N,Apw,[2,2]);           % tensor product of 2 factors, pointwise factor
check_struct(Cc,[1,2],[1,1],{N.deps.out;Apw.deps.out},'tensor product o pointwise factor');
Cd = otimes(Apw,N,[2,2]);           % pointwise factor, tensor product of 2 factors
check_struct(Cd,[1,2],[1,1],{Apw.deps.out;N.deps.out},'pointwise factor o tensor product');
% values (these use the evaluation of nested operators by apply_op)
check_equal(apply_op(Ca,{v1,v2,v3}),rename_copy(pw,'y',1)*rename_copy(b1,'y',2),'pointwise factor o single factor');
check_equal(apply_op(Cb,{v3,v1,v2}),rename_copy(b1,'y',1)*rename_copy(pw,'y',2),'single factor o pointwise factor');
check_equal(apply_op(Cc,{v1,v3,v1,v2}),rename_copy(a1,'y',1)*rename_copy(b1,'y',2)*rename_copy(pw,'y',3),'tensor product o pointwise factor');
check_equal(apply_op(Cd,{v1,v2,v1,v3}),rename_copy(pw,'y',1)*rename_copy(a1,'y',2)*rename_copy(b1,'y',3),'pointwise factor o tensor product');
fprintf('Test_tensopvar_otimes passed (%d trials, otype(1) in {1,2} alternating; products of different types).\n',ntrials)


function p = rename_copy(p,nm,k)
% The output variable nm of the polynomial p, renamed to its k-th copy nm_k
p.varname(strcmp(p.varname,nm)) = {[nm,'_',num2str(k)]};
end

function check_struct(C,sz_ops,is_nested,deps_out,label)
% The factors of C: size of ops, which of them are tensopvar objects (nested), and a row of
% depmat.out for each factor of the layer, with the dependencies deps_out{j} of factor j
if ~isequal(size(C.ops),sz_ops) || ~isequal(cellfun('isclass',C.ops(1,:),'tensopvar'),logical(is_nested))
    error("Test_tensopvar_otimes failed (%s): wrong factors, size(ops)=%s, nested=%s.",label,mat2str(size(C.ops)),mat2str(cellfun('isclass',C.ops(1,:),'tensopvar')))
end
if ~isequal(C.depmat.out,vertcat(deps_out{:}))
    error("Test_tensopvar_otimes failed (%s): depmat.out must have a row for each of the %d factors of the layer.",label,numel(deps_out))
end
end
function check_equal(Y,E,label)
err = cleanpoly(Y-E,1e-10);
if max(max(abs(err.C)))>1e-8
    error("Test_tensopvar_otimes failed (%s): residual %.3e.",label,max(max(abs(err.C))))
end
end
