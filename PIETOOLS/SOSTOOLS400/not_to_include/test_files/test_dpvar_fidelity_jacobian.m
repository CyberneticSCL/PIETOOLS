% Code testing proper functioning of jacobian
% The code builds a set of random dpvar objects and the associated polynomial
% objects, and checks whether the function produces the same results as the
% polynomial operation

% Note that the size of the problems are not exceptionally large, since we
% are not testing for speed here.

clear
clc

test = 0;
ntests = 100;
tol = 1e-15;    % tolerance in difference between dpvar and poly int

if test==0
    test_array = (1:4);
else
    test_array = test;
end

disp('Testing jacobian:')
disp('------------------------------')
for test = test_array
if test==1
    disp('testing standard case...')
    nd_max1 = 10;    % max num of dvars
    np_max1 = 10;    % max num of pvars
    mdim_max1 = 5;   ndim_max1 = 5;   % max dims of object
    deg_max1 = 5;    % max deg of object
    Cdensity_max1 = 0.8;   % maximum density of sparse matrix C, use integer value >1 to enforce full matrix
    Degdensity_max1 = 1; % maximum density of sparse matrix degmat, use integer value >1 to enforce ordinary matrix (not necessarily full)
elseif test==2
    disp('testing with no dvars...')
    nd_max1 = 0;
    np_max1 = 10;
    mdim_max1 = 5;   ndim_max1 = 5;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 0.8;
elseif test==3
    disp('testing with 1x1 object...')
    nd_max1 = 10;
    np_max1 = 10;
    mdim_max1 = 1;   ndim_max1 = 1;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test==4
    disp('testing with 1x1 object and no dvars...')
    nd_max1 = 0;
    np_max1 = 10;
    mdim_max1 = 1;   ndim_max1 = 1;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
end

% % % Test some random examples % % %
for testnumber = 1:ntests

% % % Creating first random object % % %

% % dpvar specifications:

% dvar indices
if nd_max1>0
    dvars1 = unique(randi(3*nd_max1,[1,nd_max1]),'stable');      % indices of decision variables
    nd1 = length(dvars1);     % number of decision vars
else
    dvars1 = [];
    nd1 = 0;
end
% pvar indices
if np_max1>0
    pvars1 = unique(randi(3*np_max1,[1,np_max1]),'stable');      % indices of decision variables
    np1 = length(pvars1);     % number of polynomial vars
else
    pvars1 = [];
    np1 = 0;
end
% dimensions of the dpvar matrix [m n]
mdim1 = randi(mdim_max1); 
ndim1 = randi(ndim_max1);     
matdim1 = [mdim1 ndim1];
% maximal degree of the pvars
if deg_max1>0
    deg1 = randi(deg_max1);       
else
    deg1 = 0;
end
% number of terms appearing in the degmat
%nt1 = randi((deg1+1)^np1);
nt1 = (deg1+1)*np1;
% density of degmat
if Degdensity_max1>1
    Degdensity1 = Degdensity_max1;
else
    Degdensity1 = Degdensity_max1*rand(1);
end
% density of coefficient matrix
if Cdensity_max1>1
    Cdensity1 = Cdensity_max1;
else
    Cdensity1 = Cdensity_max1;%*rand(1);
end

% Construct a random dpvar
[D1] = random_dpvar(dvars1,pvars1,deg1,nt1,Degdensity1,Cdensity1,mdim1,ndim1);


%--------------------------------------------------------------------------


% % % Creating associated polynomial % % %
D1_poly = dpvar2poly(D1);


%--------------------------------------------------------------------------


% % % Test the function % % %
pvarname = D1.varname;
np = length(pvarname);
nsub = randi(np);
x_indx1 = unique(randi(np,nsub,1));         
x_indx2 = unique(randi(np,2*nsub,1)); 
x_indx2 = x_indx2(1:min([length(x_indx2),length(x_indx1)]));
x_indx1 = x_indx1(1:length(x_indx2));
var1_cell = pvarname(x_indx1);
var1_poly = polynomial(var1_cell);


var3 = pvarname(randperm(np)');
var3_poly = polynomial(var3);

pvar dum1 dum2 dum3
var4 = [dum1;var3_poly(1:floor(end/2));dum2;(var3_poly(ceil(end/2):end));dum3];

%%%%

D1_sub_full = [];
D4_sub_full = [];
for col = 1:ndim1

D1_sub1 = jacobian(D1(:,col),pvarname);
D1_sub1p = jacobian(D1(:,col),polynomial(pvarname));
D1_poly_sub = jacobian(D1_poly(:,col),polynomial(pvarname));
error1 = dpvar2poly(D1_sub1)-D1_poly_sub;
if max(error1.coeff,[],'all') >= tol
    error('jacobian using {var1}={var1} does not produce the desired polynomial')
end
error1p = dpvar2poly(D1_sub1p)-D1_poly_sub;
if max(error1p.coeff,[],'all') >= tol
    error('jacobian using var1=var1 does not produce the desired polynomial')
end

D2_sub1 = jacobian(D1(:,col),var1_cell);
D2_sub1p = jacobian(D1(:,col),polynomial(var1_poly));
D2_poly_sub = jacobian(D1_poly(:,col),var1_poly);
error2 = dpvar2poly(D2_sub1)-D2_poly_sub;
if max(error2.coeff,[],'all') >= tol
    error('jacobian using {var1}={var1} does not produce the desired polynomial')
end
error2p = dpvar2poly(D2_sub1p)-D2_poly_sub;
if max(error2p.coeff,[],'all') >= tol
    error('jacobian using var1=var1 does not produce the desired polynomial')
end

D3_sub1 = jacobian(D1(:,col),var3);
D3_sub1p = jacobian(D1(:,col),polynomial(var3_poly));
D3_poly_sub = jacobian(D1_poly(:,col),var3_poly);
error3 = dpvar2poly(D3_sub1)-D3_poly_sub;
if max(error3.coeff,[],'all') >= tol
    error('jacobian using {var1}={var1} does not produce the desired polynomial')
end
error3p = dpvar2poly(D3_sub1p)-D3_poly_sub;
if max(error3p.coeff,[],'all') >= tol
    error('jacobian using var1=var1 does not produce the desired polynomial')
end

D4_sub = jacobian(D1(:,col),var4);
D4_poly_sub = jacobian(D1_poly(:,col),var4);
error4p = dpvar2poly(D4_sub)-D4_poly_sub;
if max(error4p.coeff,[],'all') >= tol
    error('jacobian with non-present variables does not produce the desired polynomial')
end

D1_sub_full = [D1_sub_full,D1_sub1];
D4_sub_full = [D4_sub_full,D4_sub];

end

D1_full_sub = jacobian(D1,pvarname);
errorfull = D1_full_sub - D1_sub_full;
if max(errorfull.C,[],'all') >= tol
    error('jacobian on dpvar with multiple rows does not produce the desired polynomial')
end

D4_full_sub = jacobian(D1,var4);
errorfull4 = D4_full_sub - D4_sub_full;
if max(errorfull4.C,[],'all') >= tol
    error('jacobian on dpvar with multiple rows does not produce the desired polynomial')
end

end

end

disp('No errors were encountered testing jacobian!')