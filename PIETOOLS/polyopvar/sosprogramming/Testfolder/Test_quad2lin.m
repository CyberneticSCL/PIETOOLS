clear
pvar s1 s1_dum
var1 = s1;
var2 = s1_dum;
dom = [0,1];
deg = 2;        % maximal monomial degree in independent variables
m = 3;      n = 2;
d1 = 3;          % number of state variables in distributed monomial
d2 = 2;
m_idx = randi(d1+d2);   % index for which operator may include a multiplier

% Declare d1 random 2-PI operators
Lops_opvar = cell(1,d1);
Lops = tensopvar();
Lops.ops{1} = cell(1,d1);
for ii=1:d1
    Lop_ii = rand_opvar([0,0;m,1],deg,var1,var2,dom);
    if ii~=m_idx
        % Allow only one operator with a nonzero multiplier
        Lop_ii.R.R0 = zeros([m,1]);
    end
    Lops_opvar{ii} = Lop_ii;
    Lops.ops{1}{ii} = dopvar2ndopvar(Lop_ii);
end

Rops_opvar = cell(1,d1);
Rops = tensopvar();
Rops.ops = cell(1,d1);
for ii=1:d2
    Rop_ii = rand_opvar([0,0;n,1],deg,var1,var2,dom);
    if ii+d1~=m_idx
        % Allow only one operator with a nonzero multiplier
        Rop_ii.R.R0 = zeros([n,1]);
    end
    Rops_opvar{ii} = Rop_ii;
    Rops.ops{1}{ii} = dopvar2ndopvar(Rop_ii);
end

Pmat = rand([m,n]);

[Pvec,idx_mat,var2] = quad2lin_term(Lops,Pmat,Rops);

% Generate random functions xi(s)
x1_tst = polynomial(zeros(1,d1));
y1_tst = cell(d1,1);
z1_tst = 1;
for ii=1:d1
    x1_tst(ii) = rand_poly([1,1],s1,3);
    y1_tst{ii} = apply_opvar(Lops_opvar{ii},x1_tst(ii));
    z1_tst = z1_tst.*y1_tst{ii};
end
x2_tst = polynomial(zeros(1,d2));
y2_tst = cell(d2,1);
z2_tst = 1;
for ii=1:d2
    x2_tst(ii) = rand_poly([1,1],s1,3);
    y2_tst{ii} = apply_opvar(Rops_opvar{ii},x2_tst(ii));
    z2_tst = z2_tst.*y2_tst{ii};
end

% Compute the integral
%   int_{a}^{b} (Pop{1}*x1)(s).*...*(Pop{1}*xd)(s) ds
fval = double(int(z1_tst'*Pmat*z2_tst,s1,dom(1),dom(2)));

% Compute the integral
% sum_{i=1}^{m} int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
%           Pvec{i}(t1,...,td)*x1(t1)...xd(td) dt_kd ... dt_k1
fval_alt = apply_functional(Pvec,[x1_tst,x2_tst],idx_mat,var2,dom);

f_err = norm(fval-fval_alt)