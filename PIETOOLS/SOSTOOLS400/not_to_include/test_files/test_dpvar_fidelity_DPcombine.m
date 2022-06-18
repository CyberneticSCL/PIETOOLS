% Code testing proper functioning of DPcompress
% The code builds a set of random dpvar objects and the associated polynomial
% objects, and checks whether the function produces the same results as the
% polynomial operation

% Note that the size of the problems are not exceptionally large, since we
% are not testing for speed here.

clear
clc

tol = 1e-13;
test = 0;

if test==0
    test1_array = (1:8)';
    testd1_array = (1:3)';
    testdeg1_array = (1:3)';
    testp1_array = (1:3)';
    test_array = [repmat(test1_array,3*3*3,1),...
                  repmat(kron(testd1_array,ones(8,1)),3*3,1),...
                  repmat(kron(kron(testdeg1_array,ones(8,1)),ones(3,1)),3,1),...
                  kron(kron(kron(testp1_array,ones(8,1)),ones(3,1)),ones(3,1))];
    ntests = 10;
else
    test_array = [1,1,2,2];
    ntests = 100;
end

disp('Testing DPcompress:')
disp('------------------------------')
for testnums = 1:size(test_array,1)
    
test = test_array(testnums,1);      % determines dimension and dvar/pvar size of D1
testd = test_array(testnums,2);     % determines what dvars should be removed
testdeg = test_array(testnums,3);   % determines what monomial terms should be removed
testp = test_array(testnums,4);     % determines what pvars should be removed
    
if test==1
    disp('testing standard case...')
    nd_max1 = 10;    % max num of dvars
    np_max1 = 10;    % max num of pvars
    mdim_max1 = 9;   ndim_max1 = 9;   % max dims of object
    deg_max1 = 5;    % max deg of object
    Cdensity_max1 = 0.8;   % (maximum) density of sparse matrix C, use integer value >1 to enforce full matrix
    Degdensity_max1 = 1; % maximum density of sparse matrix degmat, use integer value >1 to enforce ordinary matrix (not necessarily full)
elseif test==2
    disp('testing with no dvars...')
    nd_max1 = 0;
    np_max1 = 10;
    mdim_max1 = 9;   ndim_max1 = 9;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test==3
    disp('testing with no pvars...')
    nd_max1 = 10;
    np_max1 = 0;
    mdim_max1 = 9;   ndim_max1 = 9;
    deg_max1 = 5;
    Cdensity_max1 = 0.8;   
    Degdensity_max1 = 1;
elseif test==4
    disp('testing with no dvars or pvars...')
    nd_max1 = 0;
    np_max1 = 0;
    mdim_max1 = 9;   ndim_max1 = 9;
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


drem = 0;
if testd==1
    disp('   ...removing no dvars...')
    drem = 0;
    drem_all = 0;
elseif testd==2
    disp('   ...removing random dvars...')
    drem = 1;
    drem_all = 0;
elseif testd==3 && test~=2 && test~=4 && test~=6 && test~=8
    disp('   ...removing all dvars...')
    drem_all = 1;
end

degrem = 0;
if testdeg==1
    disp('      ...removing no monomial terms...')
    degrem = 0;
    degrem_all = 0;
elseif testdeg==2
    disp('      ...removing random monomial terms...')
    degrem = 1;
    degrem_all = 0;
elseif testdeg==3 && test~=3 && test~=4 && test~=7 && test~=8
    disp('      ...removing all monomial terms...')
    degrem_all = 1;
end

prem = 0;
if testp==1
    disp('         ...and removing no pvars.')
    prem = 0;
    prem_all = 0;
elseif testp==2
    disp('         ...and removing random pvars.')
    prem = 1;
    prem_all = 0;
elseif testp==3 && test~=3 && test~=4 && test~=7 && test~=8
    disp('         ...and removing all pvars.')
    prem_all = 1;
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


% Artificially create non-unique terms, dvars, and pvars
C = D1.C;
degmat = D1.degmat;
dvarname = D1.dvarname;
pvarname = D1.varname;

if nd1>=2 && drem
    if drem_all==1
        dvar_rem = (1:nd1);
        dvar_ret = 1;
    else
        nd_rem = max(randi(nd1)-1,1);                        % maximum number of dvars to remove
        dvar_rem = unique(randi(nd1,[1,nd_rem]));   % indices of dvars to remove
        drest = setdiff(1:nd1,dvar_rem);
        nd_rem = length(dvar_rem);  ndrest = length(drest);
        dvar_ret = kron(drest,ones(1,ceil(nd_rem/ndrest)));
        dvar_ret = dvar_ret(1:nd_rem);
    end
    dvarname(dvar_rem) = dvarname(dvar_ret);
end

ntt = size(degmat,1);
if ntt>=2 && degrem
    if degrem_all==1
        degrow_ret = 1;
        degrow_rem = (1:ntt);
        degmat = repmat(degmat(1,:),ntt,1);
    else
        nmon_rem = max(randi(ntt-1)-1,1);                        % maximum number of monomials to remove
        degrow_rem = unique(randi(ntt,[1,nmon_rem]));   % indices of rows in degmat to remove
        rest = setdiff(1:ntt,degrow_rem);
        ndeg_rem = length(degrow_rem);  nrest = length(rest);
        degrow_ret = kron(rest,ones(1,ceil(ndeg_rem/nrest)));
        degrow_ret = degrow_ret(1:ndeg_rem);
        degmat(degrow_rem,:) = degmat(degrow_ret,:);
    end
end

if np1>=2 && prem
    if prem_all==1
        pvar_rem = (1:np1);
        pvar_ret = 1;
    else
        np_rem = max(randi(np1)-1,1);                        % maximum number of dvars to remove
        pvar_rem = unique(randi(np1,[1,np_rem]));   % indices of dvars to remove
        prest = setdiff(1:np1,pvar_rem);
        np_rem = length(pvar_rem);  nprest = length(prest);
        pvar_ret = kron(prest,ones(1,ceil(np_rem/nprest)));
        pvar_ret = pvar_ret(1:np_rem);
    end
    pvarname(pvar_rem) = pvarname(pvar_ret);
end

D1 = dpvar(C,degmat,pvarname,dvarname,D1.matdim);

% % % Creating associated polynomial % % %
D1_poly = dpvar2poly(D1);


%--------------------------------------------------------------------------



% % % Test the function % % %

%[newC,newdvar] = DPVuniquevar(C,dvarname);
%[newnewC,newdeg] = DPVuniqueterm(newC,degmat);
%[newnewdeg,newpvar] = PVuniquepvar(newdeg,pvarname);
%D1new = dpvar(newnewC,newnewdeg,newpvar,newdvar,D1.matdim);

D1new = DPcombine_extended(D1);

newdvar = D1new.dvarname;
newdeg = D1new.degmat;
newpvar = D1new.varname;

if nd1>1 && drem && ~isequal(unique(D1.dvarname,'sorted'),newdvar)
    error('DPcompress did not remove all double dvars!')
end

if ntt>1 && degrem && ~isequal(size(unique(D1.degmat,'rows'),1),size(newdeg,1)) %&& size(newdeg,1)~=nrest
    error('DPcompress did not remove all double monomial terms!')
end

if np1>1 && prem && ~isequal(unique(D1.varname,'sorted'),newpvar)   % && length(newpvar)~=nprest
    error('DPcompress did not remove all double pvars!')
end

error1 = dpvar2poly(D1new) - D1_poly;
if max(error1.coeff,[],'all') >= tol
    error('DPcompress returns a different polynomial')
end


end

end

disp('No errors were encountered testing DPcompress!')