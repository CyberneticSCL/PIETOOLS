% Code testing proper functioning of subs
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

disp('Testing subs:')
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
var1_cell = pvarname(x_indx1);      var2_cell = pvarname(x_indx2);
var1_poly = polynomial(var1_cell);  var2_poly = polynomial(var2_cell);

val_indx = unique(randi(length(var2_poly),1,round(0.5*nsub)));
vals = 2*rand(length(val_indx),1) - 1;
var2_sub = subs(var2_poly,var2_poly(val_indx),vals);

var3 = pvarname(randperm(np)');
var3_poly = polynomial(var3);

D1_sub1 = subs(D1,pvarname,pvarname);
D1_sub1p = subs(D1,polynomial(pvarname),polynomial(pvarname));
error1 = dpvar2poly(D1_sub1-D1);
if max(error1.coeff,[],'all') >= tol
    error('subs using with {var1}={var1} does not produce the desired polynomial')
end
error1p = dpvar2poly(D1_sub1p-D1);
if max(error1p.coeff,[],'all') >= tol
    error('subs using with var1=var1 does not produce the desired polynomial')
end

D1_sub2 = subs(D1,var1_cell,var2_cell);
D1_sub2p = subs(D1,var1_cell,var2_poly);
D1_subp2 = subs(D1,var1_poly,var2_cell);
D1_subpp = subs(D1,var1_poly,var2_poly);
D1_sub2_poly = subs(D1_poly,var1_poly,var2_poly);

error2 = dpvar2poly(D1_sub2-D1_sub2_poly);
if max(error2.coeff,[],'all') >= tol
    error('subs using random pvarnames does not produce the desired polynomial')
end
error2p = dpvar2poly(D1_sub2p-D1_sub2_poly);
if max(error2p.coeff,[],'all') >= tol
    error('subs using with {var1} = pvar2 does not produce the desired polynomial')
end
errorp2 = dpvar2poly(D1_subp2-D1_sub2_poly);
if max(errorp2.coeff,[],'all') >= tol
    error('subs using with pvar1 = {var2} does not produce the desired polynomial')
end
errorpp = dpvar2poly(D1_subpp-D1_sub2_poly);
if max(errorpp.coeff,[],'all') >= tol
    error('subs using with pvar1 = pvar2 does not produce the desired polynomial')
end


D1_sub2val = subs(D1,var1_poly,var2_sub);
D1_sub2_poly_val = subs(D1_poly,var1_poly,var2_sub);
error2_sub = dpvar2poly(D1_sub2val-D1_sub2_poly_val);
if max(error2_sub.coeff,[],'all') >= tol
    error('subs using actual values does not produce the desired polynomial')
end


D1_sub3 = subs(D1,pvarname,var3_poly);
D1_sub3_poly = subs(D1_poly,pvarname,var3_poly);
error3 = dpvar2poly(D1_sub3-D1_sub3_poly);
if max(error3.coeff,[],'all') >= tol
    error('subs using full pvarname does not produce the desired polynomial')
end


end

end

disp('No errors were encountered testing subs!')