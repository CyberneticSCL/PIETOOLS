% Code testing proper functioning of plus and minus
% The code builds a set of random dpvar objects and the associated polynomial
% objects, and checks whether the function performs correctly in each case

% Note that the size of the problems are not exceptionally large, since we
% are not testing for speed here.

clear
clc

test = 0;
if test==0  % run through (almost) all the tests
    test1_array = (1:8)';
    test2_array = (1:8)';
    testm1_array = 1;
    testm2_array = (1:4)';
    testp2_array = (1:9)';
    test_array = [repmat(test1_array,(8*4*9),1),...
                  repmat(kron(test2_array,ones(8,1)),4*9,1),...
                  ones(8*8*4*9,1),...
                  repmat(kron(kron(testm2_array,ones(8,1)),ones(8,1)),9,1),...
                  kron(kron(kron(testp2_array,ones(8,1)),ones(8,1)),ones(4,1))];
    %test_array_m1 = [ones(8,1),test2_array,4*ones(8,1),ones(8,1),ones(8,1)];
    %test_array_m2 = [test1_array,ones(8,1),ones(8,1),4*ones(8,1),ones(8,1)];
    %test_array = [test_array; test_array_m1; test_array_m2];
    ntests = 3;
else % run through some test case
    test_array = [5,5,1,2,1];
    ntests = 100;
end

disp('Testing plus and minus:')
disp('==============================')
for testnums = 1:size(test_array,1)
    
test1 = test_array(testnums,1);   % determines dimension and dvar/pvar size of D1
test2 = test_array(testnums,2);   % determines dimension and dvar/pvar size of D2
testm1 = test_array(testnums,3);  % determines whether D1 is dpvar, polynomial, or matrix
testm2 = test_array(testnums,4);  % determines whether D2 is dpvar, polynomial, or matrix
testp2 = test_array(testnums,5);  % determines whether D2 has same of different dvars and pvars

if test1==1
    disp('------------------------------')
    disp('testing standard case...')
    nd_max1 = 10;    % max num of dvars in object 1
    np_max1 = 10;    % max num of pvars in object 1
    mdim_max1 = 5;   ndim_max1 = 5;   % max dims of object 1
    deg_max1 = 5;    % max deg of object 1
    Cdensity_max1 = 0.8;   % maximum density of sparse matrix C, use integer value >1 to enforce full matrix
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

if testm1==1
    disp('   ...where D1 is a dpvar...')
    ismat1 = 0;   ispoly1 = 0;   ispoly1_alt = 0;
elseif testm1==2
    disp('   ...where D1 is polynomial...')
    ismat1 = 0;   ispoly1 = 1;   ispoly1_alt = 0;
elseif testm1==3
    disp('   ...where D1 is dpvar converted to polynomial...')
    ismat1 = 0;   ispoly1 = 0;   ispoly1_alt = 1;
elseif testm1==4
    disp('   ...where D1 is a matrix...')
    ismat1 = 1;   ispoly1 = 0;   ispoly1_alt = 0;
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
elseif test2==5
    disp('      ...with 1x1 object D2...')
    nd_max2 = 10;
    np_max2 = 10;
    mdim_max2 = 1;   ndim_max2 = 1;
    deg_max2 = 5;
    Cdensity_max2 = 0.8;   
    Degdensity_max2 = 1;
elseif test2==6
    disp('      ...with 1x1 object and no dvars in D2...')
    nd_max2 = 0;
    np_max2 = 10;
    mdim_max2 = 1;   ndim_max2 = 1;
    deg_max2 = 5;
    Cdensity_max2 = 0.8;   
    Degdensity_max2 = 1;
elseif test2==7
    disp('      ...with 1x1 object and no pvars in D2...')
    nd_max2 = 10;
    np_max2 = 0;
    mdim_max2 = 1;   ndim_max2 = 1;
    deg_max2 = 5;
    Cdensity_max2 = 0.8;   
    Degdensity_max2 = 1;
elseif test2==8
    disp('      ...with 1x1 object and no dvars or pvars in D2...')
    nd_max2 = 0;
    np_max2 = 0;
    mdim_max2 = 1;   ndim_max2 = 1;
    deg_max2 = 5;
    Cdensity_max2 = 0.8;   
    Degdensity_max2 = 1;
end

if testm1==4 || testm2==4   % if either object is a matrix
    testp2 = 0;
    iseq_dvar = 0;   iseq_pvar = 0;   iseq_degmat = 0;
    isneq_dvar = 0;  isneq_pvar = 0;
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

if testm2==1 || testm1~=1
    disp('            ...and where D2 is a dpvar.')
    ismat2 = 0;   ispoly2 = 0;   ispoly2_alt = 0;
elseif testm2==2
    disp('            ...and where D2 is polynomial.')
    ismat2 = 0;   ispoly2 = 1;   ispoly2_alt = 0;
elseif testm2==3
    disp('            ...and where D2 is dpvar converted to polynomial.')
    ismat2 = 0;   ispoly2 = 0;   ispoly2_alt = 1;
elseif testm2==4
    disp('            ...and where D2 is a matrix.')
    ismat2 = 1;   ispoly2 = 0;   ispoly2_alt = 0;
end





% % % Test some random examples % % %
for testnumber = 1:ntests

% % % Creating first random object % % %
if ~ispoly1 && ~ismat1
% % % D1 is dpvar % % %

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

% Convert to polynomial if applicable
if ispoly1_alt
    D1 = dpvar2poly(D1);
end

elseif ispoly1
% % D1 is polynomial % %    
% pvar indices
if np_max1>0
    pvars1 = unique(randi(3*np_max1,[1,np_max1]),'stable');      % indices of decision variables
    np1 = length(pvars1);     % number of polynomial vars
else
    pvars1 = [];
    np1 = 0;
end
% dimensions of the polynomial matrix [m n]
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
% Build the polynomial
D1 = random_polynomial(pvars1,deg1,nt1,Degdensity1,Cdensity1,mdim1,ndim1);

elseif ismat1
% % D1 is a matrix % %
% Create random coefficient matrix
mdim1 = randi(mdim_max1); 
ndim1 = randi(ndim_max1);     
if Cdensity_max1>=1
    D1 = rand(mdim1,ndim1);
else
    Cdensity1 = Cdensity_max1;%*rand(1);
    D1 = sprand(mdim1,ndim1,Cdensity1);
end 

end


%--------------------------------------------------------------------------


% % % Creating second random object % % %
if ~ispoly2 && ~ismat2
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
mdim2 = mdim1;              %mdim2 = randi(mdim_max1); 
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

% Convert to polynomial if applicable
if ispoly2_alt
    D2 = dpvar2poly(D2);
end

elseif ispoly2
% % D2 is polynomial % %    
% pvar indices
if iseq_pvar
    pvars2 = circshift(pvars1,randi(max([np1,1])));   % rearrange the dvars a bit
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
% dimensions of the polynomial matrix [m n]
mdim2 = mdim1;              %mdim2 = randi(mdim_max1); 
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
% Build the polynomial
D2 = random_polynomial(pvars2,deg2,nt2,Degdensity2,Cdensity2,mdim2,ndim2);
if iseq_degmat
    degmat2 = D1.degmat;
    varname2 = D2.varname;
    nterms = size(degmat2,1);
    if Cdensity2>=1
        coeff2 = rand(nterms,mdim2*ndim2);
    else
        coeff2 = sprand(nterms,mdim2*ndim2,Cdensity2);
    end
    D2 = polynomial(coeff2,degmat2(:,end-np1+1:end),varname2,matdim2);
end

elseif ismat2
% % D2 is a matrix % %
% Create random coefficient matrix
mdim2 = mdim1;              %mdim2 = randi(mdim_max1); 
ndim2 = ndim1;              %ndim2 = randi(ndim_max2);     
if Cdensity_max2>=1
    D2 = rand(mdim2,ndim2);
else
    Cdensity2 = Cdensity_max2;%*rand(1);
    D2 = sprand(mdim2,ndim2,Cdensity2);
end 

end


%--------------------------------------------------------------------------


% % % Take the sum and compare to polynomial sum % % %
sum1 = D1 + D2;
if ~ispoly1 && ~ispoly1_alt && ~ismat1
    D1_poly = dpvar2poly(D1);
else
    D1_poly = D1;
end
if ~ispoly2 && ~ispoly2_alt && ~ismat2
    D2_poly = dpvar2poly(D2);
else
    D2_poly = D2;
end
sum1_poly = D1_poly + D2_poly;
try 
    error1 = double(dpvar2poly(sum1) - sum1_poly);
catch
    error('D1 + D2 does not produce the desired polynomial')
end
if norm(error1)~=0
    error('D1 + D2 does not produce the desired polynomial')
end

sum2 = D2 + D1;
try 
    error2 = double(dpvar2poly(sum2) - sum1_poly);
catch
    error('D2 + D1 does not produce the desired polynomial')
end
if norm(error2)~=0
    error('D2 + D1 does not produce the desired polynomial')
end

% % % Take the difference and compare to polynomial difference % % %
dif1 = D1 - D2;
dif1_poly = D1_poly - D2_poly;
try 
    error3 = double(dpvar2poly(dif1) - dif1_poly);
catch
    error('D1 - D2 does not produce the desired polynomial')
end
if norm(error3)~=0
    error('D1 - D2 does not produce the desired polynomial')
end

dif2 = D2 - D1;
try 
    error4 = double(dpvar2poly(dif2) + dif1_poly);
catch
    error('D2 - D1 does not produce the desired polynomial')
end
if norm(error4)~=0
    error('D2 - D1 does not produce the desired polynomial')
end

end

end

disp('No errors were encountered testing plus and minus!')