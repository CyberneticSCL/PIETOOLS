% Code testing proper functioning of varswap
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

disp('Testing varswap:')
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
x_indx1 = randi(np);        x_indx2 = randi(np);
var1 = pvarname{x_indx1};   var2 = pvarname{x_indx2};

D1_swpd0 = varswap(D1,var1,var1);
D1_swpd1 = varswap(D1,var1,var2);
D1_swpd2 = varswap(D1,pvar(var1),pvar(var2));
D1_poly_swpd = var_swap(D1_poly,pvar(var1),pvar(var2));

error1 = dpvar2poly(D1_swpd1) - D1_poly_swpd;
if max(error1.coeff,[],'all') >= tol
    error('varswap using varnames does not produce the desired polynomial')
end

error2 = dpvar2poly(D1_swpd2) - D1_poly_swpd;
if max(error2.coeff,[],'all') >= tol
    error('varswap using pvars does not produce the desired polynomial')
end

error0 = dpvar2poly(D1_swpd0 - D1);
if max(error0.coeff,[],'all') >= tol
    error('varswap with var1=var2 does not produce the desired polynomial')
end

if length(pvarname)>=4
    var1 = pvarname{1};
    var2 = pvarname{2};
    var3 = pvarname{3};
    var4 = pvarname{4};
end
D1_swpd3 = varswap(D1,{var1;var3},{var2;var4});
D1_swpd4 = varswap(D1,[pvar(var1);pvar(var3)],[pvar(var2);pvar(var4)]);
D1_poly_swpd3 = var_swap(var_swap(D1_poly,pvar(var1),pvar(var2)),pvar(var3),pvar(var4));

error3 = dpvar2poly(D1_swpd3) - D1_poly_swpd3;
if max(error3.coeff,[],'all') >= tol
    error('varswap with multiple varnames does not produce the desired polynomial')
end
error4 = dpvar2poly(D1_swpd4) - D1_poly_swpd3;
if max(error4.coeff,[],'all') >= tol
    error('varswap with multiple pvars does not produce the desired polynomial')
end

end

end

disp('No errors were encountered testing varswap!')