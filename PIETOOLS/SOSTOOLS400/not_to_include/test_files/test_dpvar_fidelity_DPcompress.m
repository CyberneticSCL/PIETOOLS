% Code testing proper functioning of DPcompress
% The code builds a set of random dpvar objects and the associated polynomial
% objects, and checks whether the function produces the same results as the
% polynomial operation

% Note that the size of the problems are not exceptionally large, since we
% are not testing for speed here.

clear
clc

tol = 1e-14;
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
    ntests = 100;
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


if testd==1
    disp('   ...removing random dvars...')
    drem = 1;
    drem_all = 0;
    Crem_all = 0;
elseif testd==2 && test~=2 && test~=4 && test~=6 && test~=8
    disp('   ...removing all dvars...')
    drem_all = 1;
    Crem_all = 0;
elseif testd==3
    disp('   ...with all zero C...')
    drem_all = 0;
    degrem = 0; degrem_all = 0;
    Crem_all = 1;
    testdeg = 0;
end

if testdeg==1
    disp('      ...removing no monomial terms...')
    degrem = 0;
    degrem_all = 0;
elseif testdeg==2
    disp('      ...removing random monomial terms...')
    degrem = 1;
    degrem_all = 0;
elseif testdeg==2 && test~=3 && test~=4 && test~=7 && test~=8
    disp('      ...removing all monomial terms...')
    degrem_all = 1;
end

if testp==1
    disp('         ...and removing no pvars.')
    prem = 0;
    prem_all = 0;
elseif testp==2
    disp('         ...and removing random pvars.')
    prem = 1;
    prem_all = 0;
elseif testp==2 && test~=3 && test~=4 && test~=7 && test~=8
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


% Artificially remove dvars and monomial terms from the object
C = D1.C;
degmat = D1.degmat;
if Crem_all
    C = spalloc(size(C,1),size(C,2),0);
elseif nd1>=1 && drem
    if drem_all==1
        dvar_rem = (1:nd1);
    else
        nd_rem = randi(nd1);                                % maximum number of dvars to remove
        dvar_rem = unique(randi(nd1,[1,nd_rem]));           % indices of dvars to remove
    end
    Crow_rem = dvar_rem + 1 + (0:nd1+1:size(C,1)-1)';       % rows in C to remove
    C(vec(Crow_rem),:) = 0;
end

ntt = size(degmat,1);
if ntt>=2 && degrem && ~Crem_all
    if degrem_all==1
        degrow_rem = (1:ntt);
    else
        nmon_rem = randi(ntt-1);                        % maximum number of monomials to remove
        degrow_rem = unique(randi(ntt,[1,nmon_rem]));   % indices of rows in degmat to remove
    end
    Ccol_rem = degrow_rem + (0:ntt:size(C,2)-1)';   % columns in C to remove
    C(:,vec(Ccol_rem)) = 0;
end

if np1>=1 && prem
    if prem_all==1
        pvar_rem = (1:np1);
    else
        np_rem = randi(np1);                        % maximum number of pvars to remove
        pvar_rem = unique(randi(np1,[1,np_rem]));   % indices of pvars to remove
    end
    degmat(:,pvar_rem) = 0;
end
%D1.C = C;

D1 = dpvar(C,degmat,D1.varname,D1.dvarname,D1.matdim);

% % % Creating associated polynomial % % %
D1_poly = dpvar2poly(D1);


%--------------------------------------------------------------------------



% % % Test the function % % %
D1new = compress(D1);

if nd1>=1 && drem && ~isempty(D1new.dvarname)
dvarname = D1.dvarname;
dvar_lost = setdiff(dvarname,D1new.dvarname);
% every dvar_rem must be in dvar_lost, but the opposite need not hold
if ~isempty(setdiff(dvarname(dvar_rem),dvar_lost))
    error('DPcompress did not remove the appropriate dvars!')
end
end

ntt = size(D1new.degmat,1);
ndn = length(D1new.dvarname)+1;
if ~isempty(D1new.C)
[Ci,Cj,~] = find(D1new.C);
if ~isempty(Ci)
dnums = setdiff(unique(mod(Ci-1,ndn)+1),1); 
MTi = setdiff(2:ndn,dnums);
if ~isempty(MTi)
    error('DPcompress did not remove all zero rows of C!')
end
degrow = unique(mod(Cj-1,ntt)+1);
MTj = setdiff(1:ntt,degrow);
if ~isempty(MTj)
    error('DPcompress did not remove all zero columns of C!')
end
end
end

if ~isempty(D1new.degmat)
[~,degj,~] = find(D1new.degmat);
MTj = setdiff(1:size(D1new.degmat,2),degj);
if ~isempty(MTj)
    error('DPcompress did not remove all zero columns of the degmat!')
end
end

if np1>=1 && ~isempty(D1new.varname) && prem
varname = D1.varname;
pvar_lost = setdiff(varname,D1new.varname);
% every dvar_rem must be in dvar_lost, but the opposite need not hold
if ~isempty(setdiff(varname(pvar_rem),pvar_lost))
    error('DPcompress did not remove the appropriate pvars!')
end
end

error1 = dpvar2poly(D1new) - D1_poly;
if max(error1.coeff,[],'all') >= tol
    error('DPcompress returns a different polynomial')
end


end

end

disp('No errors were encountered testing DPcompress!')