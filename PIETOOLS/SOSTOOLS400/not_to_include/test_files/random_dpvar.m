function [D] = random_dpvar(dvars,pvars,deg,nt,Degdensity,Cdensity,mdim,ndim)
% This function builds a random dpvar based on the following inputs:
%
% dvars: Indices for the decision variables. If dindcs(1) = i, then the
% first dvar will be 'coeff_i'. Indices must be unique!
% pvars: Indices for the polynomial variables. If pindcs(1) = j, then the
% first pvar will be 'x_j'. Indices must be unique!
% deg: maximal degree of each pvar appearing in the final object.
% nt: number of terms appearing in the degmat.
% Degdensity: Value in [0,1], indicating density of the sparse degmat. Set
% Degdensity = 1 to include all possible monomials up to degree deg. Set
% Degdensity > 1 to enforce full degmat with individual degrees at most deg
% (possibly expensive!)
% Cdensity: Value in [0,1], indicating density of the sparse matrix C. Set
% Cdensity > 1 to enforce a full rather than sparse matrix.
% mdim: Number of rows in the final object
% ndim: Number of columns in the final object
%
% The randomness comes from the coefficients used in Cdensity and the
% degrees used in degmat.   

% Object dimensions
matdim = [mdim,ndim];

% Construct decision variables
nd = length(dvars);
if nd>0
    cd = cell(nd,1);
    for i=1:length(dvars)
        cd{i} = ['coeff_', num2str(dvars(i))];
    end
    dvarname = cd;
else
    dvarname = {};
end

% Construct polynomial variables and degmat
np = length(pvars);
if np>0
    cp = cell(np,1);
    for i=1:length(pvars)
        cp{i} = ['x_', num2str(pvars(i))];
    end
    pvarname = cp;
    
    if Degdensity==1
        Z = monomials(pvarname,0:deg);  % note that here, deg indicates maximal degree of monomials, rather than individual vars
        degmat = Z.degmat;
    elseif Degdensity>1
        degmat = unique(randi(deg+1,[nt,np])-1,'rows','stable');
        zcols = (sum(degmat,1)==0);
        zrows = (sum(degmat,2)==0);
        if sum(zrows>1)
            degmat = [degmat(~zrows,:);zeros(1,np)];
        end
        degmat(1,zcols) = 1;
    else
        degmat = unique(ceil(deg*sprand(nt,np,Degdensity)),'rows','stable');
        [~,degj,~] = find(degmat);
        zcols = setdiff(vec(1:np),vec(degj));
        degmat(1,zcols) = 1;
    end
    nterms = size(degmat,1);
else
    pvarname = {};
    degmat = zeros(1,0);
    nterms = 1;
end

% Create random coefficient matrix
if Cdensity>=1
    C = rand((nd+1)*mdim,nterms*ndim);
else
    C = sprand((nd+1)*mdim,nterms*ndim,Cdensity);
end

% Construct the dpvar
D = dpvar(C,degmat,pvarname,dvarname,matdim);
