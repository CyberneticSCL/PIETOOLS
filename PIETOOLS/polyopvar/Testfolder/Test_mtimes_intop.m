clear
pvar s1 s1_dum s2 s3 s4 s5 s6
vartab = [s1;s2;s3;s4;s5];
var1 = s1;
var2 = s1_dum;
dom = [0,0.5];
pdeg = 2;        % maximal monomial degree in independent variables
m = 3;      n = 1;
nvars = 2;
nvars2 = 2;
% Declare distributed monomial degrees
d1 = 2;
degmat1 = randi(d1+1,[nchoosek(d1+nvars,d1),nvars])-1;
degmat1 = unique(degmat1,'rows');
d2 = 1;
degmat2 = randi(d1+1,[nchoosek(d1+nvars,d1),nvars2])-1;
degmat2 = unique(degmat2,'rows');

%m_idx = 0;
%m_idx = 2;

% Declare state variable names
xvarname = mat2cell([repmat('x',[nvars,1]),num2str((1:nvars)')],ones(nvars,1));
xvarname2 = mat2cell([repmat('y',[nvars2,1]),num2str((1:nvars2)')],ones(nvars2,1));
nshare = 1;
share_idx1 = randi(nvars,[1,nshare]);
share_idx1 = unique(share_idx1);
nshare = numel(share_idx1);
share_idx2 = randi(nvars2,[1,nshare]);
share_idx2 = unique(share_idx2);
nshare = numel(share_idx2);
share_idx1 = share_idx1(1:nshare);
xvarname2(share_idx2) = xvarname(share_idx1);

% Declare random intop objects acting ont he different monomials
C1 = tensopvar();
for j=1:size(degmat1,1)
    dj = sum(degmat1(j,:));
    if dj==0
        C1.ops{j} = rand([m,1]);
    else
        C1.ops{j} = rand_intop([m,1],vartab(1:dj),dom,pdeg);
    end
end
C2 = tensopvar();
for j=1:size(degmat2,1)
    dj = sum(degmat2(j,:));
    if dj==0
        C2.ops{j} = rand([m,1]);
    else
        C2.ops{j} = rand_intop([m,1],vartab(1:dj),dom,pdeg);
    end
end

% Declare the distributed polynomials
f1 = polyopvar();
f1.pvarname = s1.varname;
f1.varname = xvarname;
f1.varmat = true(numel(xvarname),1);
f1.dom = dom;
f1.degmat = degmat1;
f1.C = C1;

f2 = polyopvar();
f2.pvarname = s1.varname;
f2.varname = xvarname2;
f2.varmat = true(numel(xvarname2),1);
f2.dom = dom;
f2.degmat = degmat2;
f2.C = C2;

% Take the product
f3 = f1*f2;

% Generate random functions xi(s)
x_tst = polynomial(zeros(1,nvars));
for ii=1:nvars
    x_tst(ii) = rand_poly([1,1],s1,3)+1;
end
x_tst2 = polynomial(zeros(1,nvars2));
for ii=1:nvars2
    x_tst2(ii) = rand_poly([1,1],s1,2)+1;
end
x_tst2(share_idx2) = x_tst(share_idx1);
x_tstf = [x_tst,x_tst2(setdiff(1:nvars2,share_idx2))];

% Evaluate f1 at x_tst
f1_val = 0;
for j=1:size(degmat1,1)
    if sum(degmat1(j,:))==0
        f1_val = f1_val+f1.C.ops{j};
    else
        yval = apply_functional(f1.C.ops{j},x_tst,degmat1(j,:));
        f1_val = f1_val + yval;
    end
end

% Evaluate f2 at x_tst2
f2_val = 0;
for j=1:size(degmat2,1)
    if sum(degmat2(j,:))==0
        f2_val = f2_val+f2.C.ops{j};
    else
        yval = apply_functional(f2.C.ops{j},x_tst2,degmat2(j,:));
        f2_val = f2_val + yval;
    end
end

% Evaluate f3 at x_tstf
f3_val = 0;
for j=1:size(f3.degmat,1)
    if sum(f3.degmat(j,:))==0
        f3_val = f3_val+f3.C.ops{j};
    else
        yval = apply_functional(f3.C.ops{j},x_tstf,f3.degmat(j,:));
        f3_val = f3_val + yval;
    end
end

% Compare
f_err = norm(double(f3_val-f1_val.*f2_val))