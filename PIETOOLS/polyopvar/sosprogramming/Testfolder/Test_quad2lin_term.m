clear
pvar s1 s1_dum s2 s2_dum
var1 = s1;
var2 = s1_dum;
dom = [0,1];
deg = 2;        % maximal monomial degree in independent variables
m = 3;      n = 2;
nvars = 3;
deg1 = randi(2,[1,nvars])-1;
if ~any(deg1)
    deg1(randi(nvars))=1;
end
deg2 = randi(2,[1,nvars])-1;
if ~any(deg2)
    deg2(randi(nvars))=1;
end
d1 = sum(deg1);
d2 = sum(deg2);
m_idx = randi(d1+d2);   % index for which operator may include a multiplier
if d1==1 && d2==1
    m_idx2 = randi(d1+d2);
else
    m_idx2 = 0;
end
%m_idx = 2;

% Declare state variable names
xvarname = mat2cell([repmat('x',[nvars,1]),num2str((1:nvars)')],ones(nvars,1));

% Declare d1 random 2-PI operators
Lops_opvar = cell(1,d1);
Lops = tensopvar();
Lops.ops{1} = cell(1,d1);
for ii=1:d1
    Lop_ii = rand_opvar([0,0;m,1],deg,var1,var2,dom);
    if ii~=m_idx && ii~=m_idx2
        % Allow only one operator with a nonzero multiplier
        Lop_ii.R.R0 = zeros([m,1]);
    end
    Lops_opvar{ii} = Lop_ii;
    if d1==1
        Lops.ops{1} = dopvar2ndopvar(Lop_ii);
    else
        Lops.ops{1}{ii} = dopvar2ndopvar(Lop_ii);
    end
end
Lmon = polyopvar();
Lmon.varname = xvarname;
Lmon.pvarname = s1.varname;
Lmon.dom = dom;
Lmon.varmat = ones(nvars,1);
Lmon.degmat = deg1;
Lmon.C = Lops;

Rops_opvar = cell(1,d2);
Rops = tensopvar();
Rops.ops{1} = cell(1,d2);
for ii=1:d2
    Rop_ii = rand_opvar([0,0;n,1],deg,var1,var2,dom);
    if ii+d1~=m_idx && ii+d1~=m_idx2
        % Allow only one operator with a nonzero multiplier
        Rop_ii.R.R0 = zeros([n,1]);
    end
    Rops_opvar{ii} = Rop_ii;
    if d2==1
        Rops.ops{1} = dopvar2ndopvar(Rop_ii);
    else
        Rops.ops{1}{ii} = dopvar2ndopvar(Rop_ii);
    end
end
Rmon = Lmon;
Rmon.degmat = deg2;
Rmon.C = Rops;

Pmat = rand([m,n]);

% Convert to linear form
%[Kop] = quad2lin_term(Pmat,Lmon,Rmon);
[Kop] = quad2lin_term(Pmat,Lmon,Rmon,dom,s1);
Pvec = Kop.params;
idx_mat = Kop.omat;
var2 = Kop.vars;

% Generate random functions xi(s)
x_tst = polynomial(zeros(1,nvars));
y1_tst = cell(d1,1);
y2_tst = cell(d2,1);
z1_tst = 1;
z2_tst = 1;
lnum = 0;       rnum = 0;
for ii=1:nvars
    x_tst(ii) = rand_poly([1,1],s1,3);
    for jj=1:deg1(ii)
        lnum = lnum+1;
        y1_tst{jj} = apply_opvar(Lops_opvar{lnum},x_tst(ii));
        z1_tst = z1_tst.*y1_tst{jj};
    end
    for jj=1:deg2(ii)
        rnum = rnum+1;
        y2_tst{jj} = apply_opvar(Rops_opvar{rnum},x_tst(ii));
        z2_tst = z2_tst.*y2_tst{jj};
    end
    %z1_tst = z1_tst.*(y1_tst{ii}.^deg1(ii));
    %z2_tst = z2_tst.*(y2_tst{ii}.^deg2(ii));
end

% Compute the integral
%   int_{a}^{b} (Pop{1}*x1)(s).*...*(Pop{1}*xd)(s) ds
fval = double(int(z1_tst'*Pmat*z2_tst,s1,dom(1),dom(2)));

% Compute the integral
% sum_{i=1}^{m} int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
%           Pvec{i}(t1,...,td)*x1(t1)...xd(td) dt_kd ... dt_k1
fval_alt = apply_functional(Kop,x_tst,deg1+deg2);

f_err = norm(double(fval-fval_alt))