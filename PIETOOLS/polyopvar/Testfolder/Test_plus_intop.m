clear
pvar s1 s1_dum s2 s2_dum
var1 = s1;
var2 = s1_dum;
dom = [0,2];

m = 2;      n = 1;
nvars1 = 2;      % Number of state variables in first polynomial
nvars2 = 2;     % Number of state variables in second polynomial
nshare = 2;     % Number of state variables in first polynomial that match those in second


% Set the monomial degrees for both terms
deg1 = randi(4,[1,nvars1])-1;
if ~any(deg1)
    deg1(randi(nvars1))=1;
end
deg2 = randi(3,[1,nvars1])-1;
%deg2 = deg1;
if ~any(deg2)
    deg2(randi(nvars1))=1;
end
d1 = sum(deg1);
d2 = sum(deg2);
m_idx = randi(d1);   % index for which operator may include a multiplier
if d1==1 && d2==1
    m_idx2 = d1+randi(d2);   % index for which operator may include a multiplier
else
    m_idx2 = 0;
end
% m_idx = 0;
% m_idx2 = 0;
% %m_idx = 2;

pdeg = 2;        % maximal monomial degree in independent variables
nshare_pvar = randi(min(d1,d2)+1)-1;

% Declare state variable names
xvarname1 = mat2cell([repmat('x',[nvars1,1]),num2str((1:nvars1)')],ones(nvars1,1));
xvarname2 = mat2cell([repmat('y',[nvars2,1]),num2str((1:nvars2)')],ones(nvars2,1));
share_idx1 = randi(nvars1,[1,nshare]);
share_idx1 = unique(share_idx1);
nshare = numel(share_idx1);
share_idx2 = randi(nvars2,[1,nshare]);
share_idx2 = unique(share_idx2);
nshare = numel(share_idx2);
share_idx1 = share_idx1(1:nshare);
xvarname2(share_idx2) = xvarname1(share_idx1);

% Declare the polynomial/dummy variable names
pvarname1 = mat2cell([repmat('s1_',[d1,1]),num2str((1:d1)')],ones(d1,1));
pvarname2 = mat2cell([repmat('s2_',[d2,1]),num2str((1:d2)')],ones(d2,1));
share_idx1_pvar = randi(d1,[1,nshare_pvar]);
share_idx1_pvar = unique(share_idx1_pvar);
%share_idx1_pvar = 1:d1;
nshare_pvar = numel(share_idx1_pvar);
share_idx2_pvar = randi(d2,[1,nshare_pvar]);
share_idx2_pvar = unique(share_idx2_pvar);
%share_idx2_pvar = d1:-1:1;
nshare_pvar = numel(share_idx2_pvar);
share_idx1_pvar = share_idx1_pvar(1:nshare_pvar);
pvarname2(share_idx2_pvar) = pvarname1(share_idx1_pvar);

% Declare the integral operators
Kop1 = rand_intop([m,n],pvarname1,dom,pdeg);
Kop2 = rand_intop([m,n],pvarname2,dom,pdeg);

% Declare the distributed polynomials
f1 = polyopvar();
f1.varname = xvarname1;
f1.C = tensopmat(Kop1);
f1.degmat = deg1;
f1.pvarname = {'s1'};
f1.varmat = true(nvars1,1);
f1.dom = dom;

f2 = polyopvar();
f2.varname = xvarname2;
f2.C = tensopmat(Kop2);
f2.degmat = deg2;
f2.pvarname = {'s2'};
f2.varmat = true(nvars2,1);
f2.dom = dom;

% Take the sum
f3 = f1+f2;
Kops3 = f3.C.ops;
xvarname3 = f3.varname;
deg3 = f3.degmat;

% Generate random functions xi(s)
x_tst1 = polynomial(zeros(1,nvars1));
for ii=1:nvars1
    x_tst1(ii) = rand_poly([1,1],s1,3);
end
x_tst2 = polynomial(zeros(1,nvars2));
for ii=1:nvars2
    x_tst2(ii) = rand_poly([1,1],s2,2);
end
% Make sure the matching states are defined by the same polynomial
x_tst2(share_idx2) = x_tst1(share_idx1);
% Set the values of the state variables in f3
x_tst3 = polynomial(zeros(1,numel(xvarname3)));
for i=1:numel(f3.varname)
    var_idx = strcmp(xvarname3{i},xvarname1);
    if any(var_idx)
        x_tst3(i) = x_tst1(var_idx);
    else
        var_idx = strcmp(xvarname3{i},xvarname2);
        x_tst3(i) = x_tst2(var_idx);
    end
end

% Evaluate the polynomials for the given choice of states
val1 = apply_functional(Kop1,x_tst1,deg1);
val2 = apply_functional(Kop2,x_tst2,deg2);
val3 = 0;
for i=1:size(deg3,1)
    val3 = val3 + apply_functional(Kops3{i},x_tst3,deg3(i,:));
end

% Compare
err = double(val3-(val1+val2))

