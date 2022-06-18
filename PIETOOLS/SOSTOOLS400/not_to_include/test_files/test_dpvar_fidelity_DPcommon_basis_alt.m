% Code testing proper functioning of DPcommon_basis
% The code builds a set of random dpvar objects, applied DPcommon_basis,
% and checks that the results produce the same polynomial, have the same
% degmats, and the same varnames and dvarnames

% Note that the size of the problems are not exceptionally large, since we
% are not testing for speed here.

clear
clc

test = 0;

if test==0  % run through all the tests
    test1_array = (1:8)';
    test2_array = (1:4)';
    testp2_array = (1:9)';
    test_array = [repmat(test1_array,(4*9),1),...
                  repmat(kron(test2_array,ones(8,1)),9,1),...
                  kron(kron(testp2_array,ones(8,1)),ones(4,1))];
    ntests = 10;
else % run through some test case
    test_array = [1,1,1];
    ntests = 100;
end

disp('Testing DPcommon_basis:')
disp('==============================')
for testnums = 1:size(test_array,1)
    
test1 = test_array(testnums,1);   % determines dimension and dvar/pvar size of D1
test2 = test_array(testnums,2);   % determines dimension and dvar/pvar size of D2
testp2 = test_array(testnums,3);  % determines whether D2 has same of different dvars and pvars

if test1==1
    disp('------------------------------')
    disp('testing standard case...')
    nd_max1 = 10;    % max num of dvars in object 1
    np_max1 = 10;    % max num of pvars in object 1
    mdim_max1 = 5;   ndim_max1 = 5;   % max dims of object 1
    deg_max1 = 5;    % max deg of object 1
    Cdensity_max1 = 0.8;   % (maximum) density of sparse matrix C, use integer value >1 to enforce full matrix
    Degdensity_max1 = 1; % maximum density of sparse matrix degmat, use integer value >1 to enforce ordinary matrix (not necessarily full)
elseif test1==2
    disp('testing with no dvars in D1...')
    nd_max1 = 0;
    np_max1 = 10;
    mdim_max1 = 5;   ndim_max1 = 5;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test1==3
    disp('testing with no pvars in D1...')
    nd_max1 = 10;
    np_max1 = 0;
    mdim_max1 = 5;   ndim_max1 = 5;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test1==4
    disp('testing with no dvars or pvars in D1...')
    nd_max1 = 0;
    np_max1 = 0;
    mdim_max1 = 5;   ndim_max1 = 5;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test1==5
    disp('testing with 1x1 object D1...')
    nd_max1 = 10;
    np_max1 = 10;
    mdim_max1 = 1;   ndim_max1 = 1;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;  
    Degdensity_max1 = 1;
elseif test1==6
    disp('testing with 1x1 object and no dvars in D1...')
    nd_max1 = 0;
    np_max1 = 10;
    mdim_max1 = 1;   ndim_max1 = 1;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test1==7
    disp('testing with 1x1 object and no pvars in D1...')
    nd_max1 = 10;
    np_max1 = 0;
    mdim_max1 = 1;   ndim_max1 = 1;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test1==8
    disp('testing with 1x1 object and no dvars or pvars in D1...')
    nd_max1 = 0;
    np_max1 = 0;
    mdim_max1 = 1;   ndim_max1 = 1;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;
    Degdensity_max1 = 1;
end

if test2==1
    disp('      ...with standard D2...')
    nd_max2 = 10;    % max num of dvars in object 2
    np_max2 = 10;    % max num of pvars in object 2
    mdim_max2 = 5;   ndim_max2 = 5;   % max dims of object 2
    deg_max2 = 5;    % max deg of object 2
    Cdensity_max2 = 0.8;   % (maximum) density of sparse matrix C, use integer value >1 to enforce full matrix
    Degdensity_max2 = 1; % maximum density of sparse matrix degmat, use integer value >1 to enforce ordinary matrix (not necessarily full)
elseif test2==2
    disp('      ...with no dvars in D2...')
    nd_max2 = 0;
    np_max2 = 10;
    mdim_max2 = 5;   ndim_max2 = 5;
    deg_max2 = 5;
    Cdensity_max2 = 0.8;   
    Degdensity_max2 = 1;
elseif test2==3
    disp('      ...with no pvars in D2...')
    nd_max2 = 10;
    np_max2 = 0;
    mdim_max2 = 5;   ndim_max2 = 5;
    deg_max2 = 5;
    Cdensity_max2 = 0.8;   
    Degdensity_max2 = 1;
elseif test2==4
    disp('      ...with no dvars or pvars in D2...')
    nd_max2 = 0;
    np_max2 = 0;
    mdim_max2 = 5;   ndim_max2 = 5;
    deg_max2 = 5;
    Cdensity_max2 = 0.8;
    Degdensity_max2 = 1;
end

if testp2==1
    disp('         ...with mixed dvars and pvars...')
    iseq_dvar = 0;   iseq_pvar = 0;   iseq_degmat = 0;
    isneq_dvar = 0;  isneq_pvar = 0;
elseif testp2==2
    disp('         ...with all dvars equal...')
    iseq_dvar = 1;   iseq_pvar = 0;   iseq_degmat = 0;
    isneq_dvar = 0;  isneq_pvar = 0;
elseif testp2==3
    disp('         ...with all pvars equal...')
    iseq_dvar = 0;   iseq_pvar = 1;   iseq_degmat = 0;
    isneq_dvar = 0;  isneq_pvar = 0;
elseif testp2==4
    disp('         ...with all pvars and degmat equal...')
    iseq_dvar = 0;   iseq_pvar = 1;   iseq_degmat = 1;
    isneq_dvar = 0;  isneq_pvar = 0;
elseif testp2==5
    disp('         ...with all dvars and pvars equal...')
    iseq_dvar = 1;   iseq_pvar = 1;   iseq_degmat = 0;
    isneq_dvar = 0;  isneq_pvar = 0;
elseif testp2==6
    disp('         ...with all dvars and pvars and degmat equal...')
    iseq_dvar = 1;   iseq_pvar = 1;   iseq_degmat = 1;
    isneq_dvar = 0;  isneq_pvar = 0;
elseif testp2==7
    disp('         ...with all dvars different...')
    iseq_dvar = 0;   iseq_pvar = 0;   iseq_degmat = 0;
    isneq_dvar = 1;  isneq_pvar = 0;
elseif testp2==8
    disp('         ...with all pvars different...')
    iseq_dvar = 0;   iseq_pvar = 0;   iseq_degmat = 0;
    isneq_dvar = 0;  isneq_pvar = 1;
elseif testp2==9
    disp('         ...with all dvars and pvars different...')
    iseq_dvar = 0;   iseq_pvar = 0;   iseq_degmat = 0;
    isneq_dvar = 1;  isneq_pvar = 1;
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
    Cdensity1 = Cdensity_max;
else
    Cdensity1 = Cdensity_max1;%*rand(1);
end

% Construct a random dpvar
[D1] = random_dpvar(dvars1,pvars1,deg1,nt1,Degdensity1,Cdensity1,mdim1,ndim1);


%--------------------------------------------------------------------------


% % % Creating second random object % % %

% % % D2 is dpvar % % %

% % dpvar specifications:

% dvar indices
if iseq_dvar
    dvars2 = circshift(dvars1,randi(max([1,nd1])));   % rearrange the dvars a bit
    nd2 = length(dvars2);
elseif nd_max2>0
    dvars2 = unique(randi(3*nd_max2,[1,nd_max2]),'stable');      % indices of decision variables
    if isneq_dvar
        dvars2 = setdiff(dvars2,dvars1);
    end
    nd2 = length(dvars2);     % number of decision vars
else
    dvars2 = [];
    nd2 = 0;
end
% pvar indices
if iseq_pvar
    pvars2 = circshift(pvars1,max([np1,1]));   % rearrange the dvars a bit
    np2 = length(pvars2);
elseif np_max2>0
    pvars2 = unique(randi(3*np_max2,[1,np_max2]),'stable');      % indices of decision variables
    if isneq_pvar
        pvars2 = setdiff(pvars2,pvars1);
    end
    np2 = length(pvars2);     % number of polynomial vars
else
    pvars2 = [];
    np2 = 0;
end
% dimensions of the dpvar matrix [m n]
mdim2 = randi(mdim_max1); 
ndim2 = ndim1;              %ndim2 = randi(ndim_max2);     
matdim2 = [mdim2 ndim2];
% maximal degree of the pvars
if deg_max2>0
    deg2 = randi(deg_max2);       
else
    deg2 = 0;
end
% number of terms appearing in the degmat
%nt2 = randi((deg2+1)^np2);
nt2 = (deg2+1)*np2;
% density of degmat
if Degdensity_max2>1
    Degdensity2 = Degdensity_max2;
else
    Degdensity2 = Degdensity_max2*rand(1);
end
% density of coefficient matrix
if Cdensity_max2>1
    Cdensity2 = Cdensity_max2;
else
    Cdensity2 = Cdensity_max2;%*rand(1);
end

% Construct a random dpvar
[D2] = random_dpvar(dvars2,pvars2,deg2,nt2,Degdensity2,Cdensity2,mdim2,ndim2);
if iseq_degmat
    degmat2 = D1.degmat;
    varname2 = D2.varname;
    dvarname2 = D2.dvarname;
    nterms = size(degmat2,1);
    if Cdensity2>=1
        C2 = rand((nd2+1)*mdim2,nterms*ndim2);
    else
        C2 = sprand((nd2+1)*mdim2,nterms*ndim2,Cdensity2);
    end
    D2 = dpvar(C2,degmat2(:,end-np1+1:end),varname2,dvarname2,matdim2);
end


%--------------------------------------------------------------------------


% % % Testing DPcommon_basis for defined dpvars % % %
% Create corresponding objects with the same bases
[D1new,D2new] = DPcommon_basis_sparse(D1,D2);

% Check difference between old dvars and new dvars, and compare the new
% degmats
try 
    error1 = double(dpvar2poly(D1new) - dpvar2poly(D1));
catch
    error('D1new does not match D1')
end
if norm(error1)~=0
    error('D1new does not match D1')
end

try 
    error2 = double(dpvar2poly(D2new) - dpvar2poly(D2));
catch
    error('D2new does not match D2')
end
if norm(error2)~=0
    error('D2new does not match D2')
end

try 
    error3 = double(D1new.degmat - D2new.degmat);
catch
    error('The produced dpvars have different degmats')
end
if nnz(error2)~=0
    error('The produced dpvars have different degmats')
end

if ~isequal(D1new.varname,D2new.varname)
    error('The produced dpvars have different polynomial varnames')
end
if ~isequal(D1new.dvarname,D2new.dvarname)
    error('The produced dpvars have different decision varnames')
end



end

end

disp('No errors were encountered testing DPcommon_basis!')