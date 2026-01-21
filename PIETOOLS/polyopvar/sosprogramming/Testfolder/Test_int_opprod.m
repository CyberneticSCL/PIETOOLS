clear
pvar s1 s1_dum
var1 = s1;
var2 = s1_dum;
dom = [0,1];
deg = 2;        % maximal monomial degree in independent variables
m = 3;      n = 1;
d = 4;          % number of state variables in distributed monomial
m_idx = randi(d);   % index for which operator may include a multiplier

% Declare d random 2-PI operators
Pops_opvar = cell(1,d);
Pops = cell(1,d);
for ii=1:d
    Pop_ii = rand_opvar([0,0;m,n],deg,var1,var2,dom);
    if ii~=m_idx
        % Allow only one operator with a nonzero multiplier
        Pop_ii.R.R0 = zeros([m,n]);
    end
    Pops_opvar{ii} = Pop_ii;
    Pops{ii} = dopvar2ndopvar(Pop_ii);
end

[Pvec,idx_mat,var2] = int_opprod(Pops);

% Generate random functions xi(s)
x_tst = polynomial(zeros(d,1));
y_tst = cell(d,1);
z_tst = 1;
for ii=1:d
    x_tst(ii) = rand_poly([1,1],s1,3);
    y_tst{ii} = apply_opvar(Pops_opvar{ii},x_tst(ii));
    z_tst = z_tst.*y_tst{ii};
end

% Compute the integral
%   int_{a}^{b} (Pop{1}*x1)(s).*...*(Pop{1}*xd)(s) ds
fval = double(int(z_tst,s1,dom(1),dom(2)));

% Compute the integral
% sum_{i=1}^{m} int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
%           Pvec{i}(t1,...,td)*x1(t1)...xd(td) dt_kd ... dt_k1
fval_alt = apply_functional(Pvec,x_tst,idx_mat,var2,dom);

f_err = norm(fval-fval_alt);