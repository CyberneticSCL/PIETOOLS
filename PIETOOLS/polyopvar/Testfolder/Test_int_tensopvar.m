clear
pvar s1 s1_dum s2 s2_dum r1 r2
var1 = s1;
var2 = s1_dum;
dom = [0,1];
pdeg = 3;        % maximal monomial degree in independent variables
m = 3;      p = 2;  n = 1;
nvars = 2;
deg1 = randi(3,[1,nvars])-1;
if ~any(deg1)
    deg1(randi(nvars))=1;
end
d1 = sum(deg1);
m_idx = randi(d1);   % index for which operator may include a multiplier
m_idx = 0;

% Declare state variable names
xvarname = mat2cell([repmat('x',[nvars,1]),num2str((1:nvars)')],ones(nvars,1));

% Declare d1 random 2-PI operators
Cops_opvar = cell(1,d1);
Cops = tensopvar();
Cops.ops{1} = cell(1,d1);
for ii=1:d1
    Lop_ii = rand_opvar([0,0;p,n],pdeg,var1,var2,dom);
    if ii~=m_idx
        % Allow only one operator with a nonzero multiplier
        Lop_ii.R.R0 = zeros([p,n]);
    end
    Cops_opvar{ii} = Lop_ii;
    if d1==1
        Cops.ops{1} = dopvar2ndopvar(Lop_ii);
    else
        Cops.ops{1}{ii} = dopvar2ndopvar(Lop_ii);
    end
end
Cmon = polyopvar();
Cmon.varname = xvarname;
Cmon.pvarname = s1.varname;
Cmon.dom = dom;
Cmon.varmat = ones(nvars,1);
Cmon.degmat = deg1;
Cmon.C = Cops;

% Set the domain over which to integrate
use_Ldom = logical(randi(2)-1);
use_Udom = logical(randi(2)-1);
int_dom = polynomial(dom);
if use_Ldom
    int_dom(1) = r1;
end
if use_Udom
    int_dom(2) = r2;
end
vars_r = polynomial(int_dom.varname);

% Generate a random kernel function
vars2_name =  mat2cell([repmat('s1_',[d1,1]),num2str((1:d1)')],ones(d1,1))';
vars2 = polynomial(vars2_name');
Kfun = rand_poly([m,p],[var1;vars_r],pdeg);

% Compute the integral of the tensor-PI operator defined by Cops with the
% kernel Kfun
[KCfun,omat,vars2f_name,dom] = int(Cops,Kfun,var1,int_dom,vars2_name);


% Generate random functions xi(s)
x_tst = polynomial(zeros(n,nvars));
y1_tst = cell(d1,1);
z1_tst = 1;
lnum = 0;
for ii=1:nvars
    x_tst(:,ii) = rand_poly([n,1],s1,3);
    for jj=1:deg1(ii)
        lnum = lnum+1;
        y1_tst{lnum} = apply_opvar(Cops_opvar{lnum},x_tst(:,ii));
        z1_tst = z1_tst.*y1_tst{lnum};
    end
end

% Compute the integral
%   int_{a}^{b} (Pop{1}*x1)(s).*...*(Pop{1}*xd)(s) ds
fval = int(Kfun*z1_tst,s1,int_dom(1),int_dom(2));
if use_Ldom
    xr1_tst = rand_poly([n,1],r1,3);
    fval = int(xr1_tst*fval,r1,dom(1),int_dom(2));
end
if use_Udom
    xr2_tst = rand_poly([n,1],r2,3);
    fval = int(xr2_tst*fval,r2,dom(1),dom(2));
end
fval = double(fval);

% Compute the integral
% sum_{i=1}^{m} int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
%           Pvec{i}(t1,...,td)*x1(t1)...xd(td) dt_kd ... dt_k1
xf_tst = x_tst;
degf = deg1;
if use_Ldom
    xf_tst = [xr1_tst,xf_tst];
    degf = [1,degf];
end
if use_Udom
    xf_tst = [xf_tst,xr2_tst];
    degf = [degf,1];
end
Kop = intvar(KCfun,omat,vars2f_name,dom);
% Get rid of duplicate terms
state_idcs = 1:size(degf,2);
state_idcs = repelem(state_idcs,degf);
Kop_alt = combine_terms(Kop,state_idcs);
fval_alt = apply_functional(Kop,xf_tst,degf);

f_err = norm(double(fval-fval_alt))