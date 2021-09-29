% Code testing proper functioning of subsref
% The code builds a set of random dpvar objects and the associated polynomial
% objects, and checks whether the function produces the same results as the
% polynomial operation

% Note that the size of the problems are not exceptionally large, since we
% are not testing for speed here.

clear
clc

test = 0;
ntests = 100;

if test==0
    test_array = (1:8);
else
    test_array = test;
end

disp('Testing subsref:')
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
D1_poly = dpvar2poly(D1);


%--------------------------------------------------------------------------


% % % Test the function % % %
m = (randi(mdim1,1,2));
m0 = min(m);    m1 = max(m);
n = (randi(ndim1,1,2));
n0 = min(n);    n1 = max(n);

D1_11 = D1(m1,n0);              D1_11_poly = D1_poly(m1,n0);
D1_d1 = D1(:,n0);               D1_d1_poly = D1_poly(:,n0);
D1_1d = D1(m0,:);               D1_1d_poly = D1_poly(m0,:);
D1_dd = D1(:,:);                D1_dd_poly = D1_poly(:,:);

D1_pp = D1(m0:m1,n0:n1);        D1_pp_poly = D1_poly(m0:m1,n0:n1);
D1_pm = D1(m0:m1,n1:-1:n0);     D1_pm_poly = D1_poly(m0:m1,n1:-1:n0);
D1_mp = D1(m1:-1:m0,n0:n1);     D1_mp_poly = D1_poly(m1:-1:m0,n0:n1);
D1_mm = D1(m1:-1:m0,n1:-1:n0);  D1_mm_poly = D1_poly(m1:-1:m0,n1:-1:n0);

rm = randperm(m1,m0);           rn = randperm(n1,n0);
rm = [rm,randperm(m1,m0)];      rn = [rn,randperm(n1,n0)];
D1_1r = D1(m1,rn);              D1_1r_poly = D1_poly(m1,rn);
D1_r1 = D1(rm,n1);              D1_r1_poly = D1_poly(rm,n1);
D1_rr = D1(rm,rn);              D1_rr_poly = D1_poly(rm,rn);

% % Indexing using "end" is not supported:
%D1_dp = D1(:,n0:end);           D1_dp_poly = D1_poly(:,n0:end);
%D1_pd = D1(m0:end,:);           D1_pd_poly = D1_poly(m0:end,:);
%D1_ee = D1(end,end);            D1_ee_poly = D1_poly(end,end);
%D1_dm = D1(:,end:-1:n0);         D1_dm_poly = D1_poly(:,end:-1:n0);
%D1_md = D1(end:-1:m0,:);         D1_md_poly = D1_poly(end:-1:m0,:);

D1_dp2poly = dpvar2poly(D1);

try 
    error11 = double(dpvar2poly(D1_11) - D1_11_poly);
    errord1 = double(dpvar2poly(D1_d1) - D1_d1_poly);
    error1d = double(dpvar2poly(D1_1d) - D1_1d_poly);
    errordd = double(dpvar2poly(D1_dd) - D1_dd_poly);
    errorpp = double(dpvar2poly(D1_pp) - D1_pp_poly);
    errorpm = double(dpvar2poly(D1_pm) - D1_pm_poly);
    errormp = double(dpvar2poly(D1_mp) - D1_mp_poly);
    errormm = double(dpvar2poly(D1_mm) - D1_mm_poly);
    error1r = double(dpvar2poly(D1_1r) - D1_1r_poly);
    errorr1 = double(dpvar2poly(D1_r1) - D1_r1_poly);
    errorrr = double(dpvar2poly(D1_rr) - D1_rr_poly);
    %errordp = double(dpvar2poly(D1_dp) - D1_dp_poly);
    %errorpd = double(dpvar2poly(D1_pd) - D1_pd_poly);
    %erroree = double(dpvar2poly(D1_ee) - D1_ee_poly);
    %errordm = double(dpvar2poly(D1_dm) - D1_dm_poly);    
    %errormd = double(dpvar2poly(D1_md) - D1_md_poly);
catch
    error('subsref does not produce the desired polynomial')
end

if norm(errordd)~=0 || norm(errord1)~=0 || norm(error1d)~=0 || norm(error11)~=0
    error('subsref does not produce the desired polynomial, check single and ":" indexing')
end
if norm(errorpp)~=0 || norm(errorpm)~=0 || norm(errormp)~=0 || norm(errormm)~=0
    error('subsref does not produce the desired polynomial, check increasing and decreasing indexing')
end
if norm(error1r)~=0 || norm(errorr1)~=0 || norm(errorrr)~=0
    error('subsref does not produce the desired polynomial, check indexing using arrays')
end


end

end

disp('No errors were encountered testing subsref!')
disp('However, at this time, indexing using "end" is not supported');