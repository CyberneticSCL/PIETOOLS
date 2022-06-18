% Code testing proper functioning of poly2dpvar
% The code builds a set of random dpvar objects and the associated polynomial
% objects, and checks whether the dpvars produced using poly2dpvar
% are correct.

% Note that the size of the problems are not exceptionally large, since we
% are not testing for speed here.

clear
clc

test = 0;
ntests = 100;
tol = 1e-15;

if test==0
    test_array = (1:8);
else
    test_array = test;
end

disp('Testing poly2dpvar:')
disp('------------------------------')
for test = test_array
if test==1
    disp('testing standard case...')
    nd_max1 = 10;    % max num of dvars
    np_max1 = 10;    % max num of pvars
    mdim_max1 = 5;   ndim_max1 = 5;   % max dims of object
    deg_max1 = 5;    % max deg of object
    Cdensity_max1 = 0.8;   % (maximum) density of sparse matrix C, use integer value >1 to enforce full matrix
    Degdensity_max1 = 1; % maximum density of sparse matrix degmat, use integer value >1 to enforce ordinary matrix (not necessarily full)
elseif test==2
    disp('testing with no dvars...')
    nd_max1 = 0;
    np_max1 = 10;
    mdim_max1 = 5;   ndim_max1 = 5;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test==3
    disp('testing with no pvars...')
    nd_max1 = 10;
    np_max1 = 0;
    mdim_max1 = 5;   ndim_max1 = 5;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test==4
    disp('testing with no dvars or pvars...')
    nd_max1 = 0;
    np_max1 = 0;
    mdim_max1 = 5;   ndim_max1 = 5;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test==5
    disp('testing with 1x1 object...')
    nd_max1 = 10;
    np_max1 = 10;
    mdim_max1 = 1;   ndim_max1 = 1;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test==6
    disp('testing with 1x1 object and no dvars...')
    nd_max1 = 0;
    np_max1 = 10;
    mdim_max1 = 1;   ndim_max1 = 1;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test==7
    disp('testing with 1x1 object and no pvars...')
    nd_max1 = 10;
    np_max1 = 0;
    mdim_max1 = 1;   ndim_max1 = 1;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test==8
    disp('testing with 1x1 object and no dvars or pvars...')
    nd_max1 = 0;
    np_max1 = 0;
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
% Build monomial vectors corresponding to the dpvar
nterms = size(D1.degmat,1);
Z_d = [1;monomials(D1.dvarname,1)];
Z_p = polynomial(eye(nterms),D1.degmat,D1.varname,[nterms,1]);

Z_dd = [];
for i=1:mdim1
    Z_dd = blkdiag(Z_dd,Z_d);
end
Z_pp = [];
for i=1:ndim1
    Z_pp = blkdiag(Z_pp,Z_p);
end

% Build the polynomial
D1_poly = Z_dd'*D1.C*Z_pp;


%--------------------------------------------------------------------------


% % % Comparing dpvar to that produced with poly2dpvar % % %

% test with no dvars specified (noting that dvars are named "coeff_")
D1_poly2dp_1 = poly2dpvar(D1_poly);
dif1 = D1 - D1_poly2dp_1;
if max(dif1.C,[],'all')>tol
    error('poly2dpvar(D1_poly) does not produce the desired dpvar')
end
if ~isempty(setdiff(D1_poly2dp_1.dvarname,D1.dvarname))
    error('poly2dpvar(D1_poly) does not have the desired dvarname')
end
D1_poly2poly_1 = dpvar2poly(D1_poly2dp_1);
dif2 = D1_poly - D1_poly2poly_1;
try
    error2 = double(dif2);
catch
    error('dpvar2poly(poly2dpvar(D1_poly)) does not produce the desired polynomial')
end
if error2>tol
    error('dpvar2poly(poly2dpvar(D1_poly)) does not produce the desired polynomial')
end
D1_dp2dp_1 = poly2dpvar(dpvar2poly(D1));
dif3 = D1 - D1_dp2dp_1;
if max(dif3.C,[],'all')>tol
    error('poly2dpvar(dpvar2poly(D1)) does not produce the desired dpvar')
end


% test with full dvars specified
dvarname = D1.dvarname;
D1_poly2dp_2 = poly2dpvar(D1_poly,dvarname);
dif1 = D1 - D1_poly2dp_2;
if max(dif1.C,[],'all')>tol
    error('poly2dpvar(D1_poly,D1.dvarname) does not produce the desired dpvar')
end
if ~isequal(dvarname,D1_poly2dp_2.dvarname)
    error('poly2dpvar(D1_poly,D1.dvarname) does not have the desired dvarname')
end
D1_dp2dp_2 = poly2dpvar(dpvar2poly(D1));
dif3 = D1 - D1_dp2dp_2;
if max(dif3.C,[],'all')>tol
    error('poly2dpvar(dpvar2poly(D1),D1.dvarname) does not produce the desired dpvar')
end

% test with partial dvar specified
nd = length(dvarname);
if nd~=0
    for i=1:nd
        nd_part = randi(nd);
        dvars = unique(randi(nd,[1,nd_part]),'stable');
        dvarname_part = dvarname(dvars);
        D1_poly2dp_3 = poly2dpvar(D1_poly,dvarname_part);
        dif1 = D1_poly - dpvar2poly(D1_poly2dp_3);
        if max(dif1.C,[],'all')>tol
            error('poly2dpvar(D1_poly,D1.dvarname(i:j)) does not produce the desired dpvar')
        end
        if ~isequal(dvarname_part,D1_poly2dp_3.dvarname)
            error('poly2dpvar(D1_poly,D1.dvarname(i:j)) does not have the desired dvarname')
        end
    end
end    
    
end

end

disp('No errors were encountered testing poly2dpvar!')