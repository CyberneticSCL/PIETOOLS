function Vx = quad2lin(Pmat,ZopL,ZxL,ZopR,ZxR)
% VX = QUAD2LIN(PMAT,ZOP,ZX) converts a distributed polynomial in quadratic
% form,
%   V(x) = <ZopL*ZL(x), Pmat*ZopR*ZR(x)>_{L2}
% to the linear form, computing a functional operator Cop such that
%   V(x) = Cop*unique(ZL(x)oZR(x))
%
% INPUTS
% - Pmat:   m x n 'double' or 'dpvar' object representing the coefficients
%           parameterizing V(x) in the quadratic format;
% - ZopL:   m x d1 'tensopvar' object representing the coefficient
%           operators acting on the distributed monomials ZxL;
% - ZxL:    d1 x 1 'polyopovar' object representing a basis of distributed
%           monomials;
% - ZopR:   n x d1 'tensopvar' object representing the coefficient
%           operators acting on the distributed monomials ZxR. Defaults to
%           ZopL if not specified;
%
% OUTPUTS
% - Vx:     1 x 1 'polyopvar' object representing the function
%               V(x) = <ZopL*ZL(x), Pmat*ZopR*ZR(x)>_{L2}
%           in the linear format;


% Assume a symmetric operator if only one half is specified.
if nargin<=3 && size(Pmat,1)==size(Pmat,2)
    ZopR = ZopL;
    ZxR = ZxL;
elseif nargin<=3
    error("Insufficient input arguments.")
end

% Extract the degrees of the distributed monomials
degmatL = ZxL.degmat;
nZL = size(degmatL,1);
degmatR = ZxR.degmat;
nZR = size(degmatR,1);



% Compute the list of monomials that appear in the distributed polynomial.
degmat_full = repmat(degmatL,[nZR,1]) + kron(degmatR,ones(nZR,1));
[M,degmat_lin] = uniquerows_integerTable(degmat_full);  % degmat_full = M*degmat_lin
old2new_idcs = M*(1:size(degmat_lin,1))';               % degmat_full = degmat_lin(old2new_idcs,:);
nZ_lin = size(degmat_lin,1);

% Get the dimensions of the different operators
blkdim = zeros(nZL,1);
for ii=1:nZL
    blkdim(ii) = ZopL(ii,ii).dim(1);
end
blkdim_cum = [0;cumsum(blkdim)];

Kop = tensopvar();
Kop.ops = cell(1,nZ_lin);
for ii=1:nZL
    % Extract the monomial operators acting on the ith distributed
    % monomials
    Zop_ii = ZopL(ii,ii);
    for jj=1:nZR
        % Extract the monomial operators acting on the jth distributed
        % monomial
        Zop_jj = ZopR(jj,jj);
        % Extract the block of decision variables acting on monomials Zi
        % and Zj
        Pij = Pmat(blkdim_cum(ii)+1:blkdim_cum(ii+1), blkdim_cum(jj)+1:blkdim_cum(jj+1));
        % Convert to the linear format
        Kij = quad2lin_term(Zop_ii,Pij,Zop_jj);
        % We need to re-order variables to account for the orderof the 
        % state variables in the distributed monomials
        %   i.e. degmat=[2,2] indicates x(s1)*x(s2)*v(s3)*v(s4)
        % but quad2lin_term returns
        %  int x(t1)*v(t2)*x(t3)*v(t4)
        % Match the obtained kernels with the correct distributed monomial
        lidx = old2new_idcs((ii-1)*nZL+jj);
        if isempty(Kop.ops{lidx})
            Kop.ops{lidx} = Kij;
        else
            for kk=1:numel(Kop.ops{lidx})
                Kop.ops{lidx}{kk} = Kij{kk};
            end
        end
    end
end

Vx = ZxL;
Vx.degmat = degmat_lin;
Vx.C = Kop;


end