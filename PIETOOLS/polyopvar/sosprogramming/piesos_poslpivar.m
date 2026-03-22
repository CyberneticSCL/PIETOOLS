function [prog,Pmat,Zop] = piesos_poslpivar(prog,Z,pdegs,opts)
% [PROG,CMAT,ZOP] = PIESOS_POSLPIVAR(PROG,Z,PDEGS,OPTS) generates a
% posisitve semidefinite matrix PMAT of decision variables and a basis of
% PI operators ZOP representing a positive semidefinite PI operator
%   Pop = ZOP' PMAT ZOP
%
% INPUTS
% - prog:   lpiprogram structure representing an LPI optimization program;
% - Z:      m x 1 'polyopvar' object representing the basis of monomials in
%           the state variables on which the operator Pop should act;
% - pdegs:  m x 1 cell, with each element pdegs{i} a N x 2 array specifying
%           the maximal monomial degree of the primary and dummy
%           independent variables to use in the operator Zop, associated
%           with the distributed monomial Z{i}. If specified as a N x 2
%           array, the same degree will be used for all distributed
%           monomials. If specified as a scalar, pdegs*[1,1/d] will be used
%           as degree for each distributed monomials, where d is the degree
%           of the distributed monomial;
% - opts:   struct with fields
%           - 'exclude': 3x1 boolean array indicating whether to exclude 
%               the multiplier term (exclude(1)=1), lower integral 
%               (exclude(2)=1), or upper integral (exclude(3)=1) from the 
%               definition of basis operator Zop;
%           - 'psatz': 1xp array of integers specifying whether to use a
%               psatz multiplier to enforce positivity only on the interval
%               over which the operator integrates. If psatz=0, no such
%               multiplier is used. If psatz=1, the matrix Pmat will be
%               multiplied by (s-a)*(b-s), where [a,b]=Z.dom. If
%               psatz=[0,1], we have Pmat = Pmat1 + (s-a)*(b-s)*Pmat2;
%           - 'sep': scalar boolean specifying whether to enforce
%               separability of the operator;
%
% OUTPUTS
% - prog:   lpiprogram structure representing the same LPI optimization
%           program as the input, but now with the positive semidefinite
%           decision variable matrix Pmat added to the structure;
% - Pmat:   n x n 'dpvar' object representing a symmetric positive
%           semidefinite decision variable operator;
% - Zop:    n x m 'tensopvar' object representing the monomial basis
%           operator acting on the distributed monomial basis defined by Z,
%           parameterizing a positive semidefinite operator as
%               Pop = Zop' * Pmat * Zop,
%           and thus parameterizing a distributed SOS functional as
%               V(x) = <Zop*Z(x) , Pmat*Zop*Z(x)>_{L2},
%           where Z(x) is the basis of distributed monomials in the state
%           variables x;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_poslpivar
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 01/19/2026: Initial coding


% Extract the degrees of the distributed monomials
xdegs = Z.degmat;
nZ = size(xdegs,1);
pvarname = Z.pvarname;
pvarname_full = [pvarname,cellfun(@(a)[a,'_dum'],pvarname,'UniformOutput',false)];
pvars = polynomial(pvarname_full);
dom = Z.dom;
N = size(pvars,1);
if N>=2
    error("Multivariate case is currently not supported.")
end


% Check that the independent variables match
vartab = polynomial(prog.vartable);
if any(~ismember(pvarname,vartab.varname))
    error("Independent variables in the distributed polynomial must match those in the LPI program.")
end
%pvars = [pvar1(1:N,1),pvars(:,2)];
var1 = pvars(:,1);


% % Check that the degrees are properly specified
if nargin<=2
    % Set default degree of 0
    pdegs = 0;
end
if isa(pdegs,'double')
    % Assume same degree of spatial variables for all distributed monomials
    pdegs = repmat({pdegs},size(xdegs,1),1);
elseif ~isa(pdegs,'cell')
    error("Degrees of spatial variables must be specified as Nx2 array of m x 1 cell of Nx2 arrays.")
elseif numel(pdegs)~=size(xdegs,1)
    error("Number of specified monomial degrees in spatial variables should match number of distributed monomials.")
end
% Make sure degrees for each distributed monomial are appropriate
for i=1:numel(pdegs)
    if ~isa(pdegs{i},'double') || any(round(abs(pdegs{i}))~=pdegs{i})
        error("Degrees of spatial variables must be specified as nonnegative integers.")
    end
    if size(pdegs{i},2)==1
        % Assume pdegs{i} is cummulative degree in all dummy variables
        dtot = sum(xdegs(i,:));
        pdegs{i} = [pdegs{i},ceil(pdegs{i}/dtot)];
    end
    if size(pdegs{i},1)==1
        % Assume same degree in all spatial variables
        pdegs{i} = repmat(pdegs{i},N,1);
    end
    if ~all(size(pdegs{i})==[N,2])
        error("Degrees of spatial variables for each distributed monomial must be specified as Nx2 array.")
    end
end

% % Set options for excluding parameters or enforcing separability
excludeL = zeros([3*ones(1,N),1]);
psatz = 0;
use_sep = 0;
if nargin>=3 && isfield(opts,'exclude')
    excludeL = opts.exclude;
end
if nargin>=3 && isfield(opts,'psatz')
    psatz = reshape(opts.psatz,1,[]);
end
if nargin>=3 && isfield(opts,'sep')
    use_sep = opts.sep;
end
if use_sep
    excludeL(3) = 1;
end

% % For each distributed monomial, generate a basis of operators acting on
% % that distributed monomial
Pdim = 0;
Pdim_arr = zeros(nZ,1);
%d_ii = prod(pdegs+1);
Zop = tensopvar();
for ii=1:nZ
    xdeg_ii = xdegs(ii,:);
    dtot = sum(xdeg_ii);
    d_ii = prod(pdegs{ii}(:,2)+1);
    if dtot==0
        % For constant monomial, just multiply with vector in s
        Zop.ops{ii,ii} = 1;
        Pdim = Pdim + 1;
        Pdim_arr(ii) = 1;
    elseif dtot==1
        % Include the basis of 3-PI operators acting on the single state
        % Zd(s) o Zd(t)
        Z_ii = nopvar();
        Z_ii.vars = pvars;
        Z_ii.dom = dom;
        Z_ii.deg = pdegs{ii}(:,2);
        Ccell = cell(3,1);
        m_arr = [d_ii; d_ii^2; d_ii^2].*(~excludeL(:));
        m_tot = sum(m_arr);
        r_shft = 0;
        if ~excludeL(1)
            % Set coefficients defining multiplier operator 1
            c_idcs = 1;
            r_idcs = 1;     % Only first term is nonzero
            Ccell{1} = sparse(r_idcs+r_shft,c_idcs,1,m_tot,1);
            r_shft = r_shft + d_ii;
        end
        if ~excludeL(2)
            % Set coefficients defining integral int_{0}^{s} dt Zd(t)
            c_idcs = (1:d_ii)';
            r_idcs = (1:d_ii:d_ii^2)';
            Ccell{2} = sparse(r_idcs+r_shft,c_idcs,1,m_tot,d_ii);
            r_shft = r_shft + d_ii^2;
            if use_sep
                Ccell{3} = Ccell{2};
            end
        end
        if ~excludeL(3) && ~use_sep
            % Set coefficients defining integral int_{s}^{1} dt Zd(t)
            c_idcs = (1:d_ii)';
            r_idcs = (1:d_ii:d_ii^2)';
            Ccell{3} = sparse(r_idcs+r_shft,c_idcs,1,m_tot,d_ii);
        end
        Z_ii.C = Ccell;
        Zop.ops{ii,ii} = Z_ii;
        Pdim = Pdim + Z_ii.dim(1);
        Pdim_arr(ii) = Z_ii.dim(1);
    else
        % The distributed monomial involves multiple state variables
        % --> include a basis of 2-PI operators acting on each state
        % variable,
        %   (Z{1}*x1)(s)*(Z{2}*x2)(s)
        Z_ii = cell(1,dtot);
        degs_ii = pdegs{ii}(:,2);
        d_ii = prod(degs_ii+1);
        m_arr = d_ii*(~excludeL(2)+~excludeL(3))*ones(1,dtot);
        m_tot = prod(m_arr)*d_ii;
        for jj=1:dtot
            % Declare a nopvar object acting on the jth factor in the
            % monomial
            Zjj = nopvar();
            Zjj.vars = pvars;
            Zjj.dom = dom;
            Zjj.deg = degs_ii;
            Ccell = cell(3,1);
            if ~excludeL(2)
                % Set the parameters of the integral int_{a}^{s}
                c_idcs = kron((1:d_ii)',ones(prod(m_arr(jj+1:end)),1));       % index of monomial in dummy var
                c_idcs = repmat(c_idcs,[prod(m_arr(1:jj-1)),1]);         % account for Kronecker product with other bases
                r_idcs = kron(1,ones(d_ii*prod(m_arr(jj+1:end)),1));    % index of monomial in primary var
                r_idcs = r_idcs + d_ii*(0:d_ii*prod(m_arr(jj+1:end))-1)';     % account for location in vector-valued object, (Im o Zd(s))
                r_idcs = r_idcs + d_ii*prod(m_arr(jj:end))*(0:prod(m_arr(1:jj-1))-1); % account for Kronecker product with other bases
                Ccell{2} = sparse(r_idcs(:),c_idcs,1,m_tot,d_ii);
                if use_sep
                    Ccell{3} = Ccell{2};
                end
            else
                Ccell{2} = sparse(m_tot,d_ii);
            end
            if ~excludeL(3) && ~use_sep
                % Set the parameters of the integral int_{s}^{b}
                c_idcs = kron((1:d_ii)',ones(prod(m_arr(jj+1:end)),1));
                c_idcs = repmat(c_idcs,[prod(m_arr(1:jj-1)),1]);
                r_idcs = kron(1,ones(d_ii*prod(m_arr(jj+1:end)),1));
                r_idcs = r_idcs + d_ii*(0:d_ii*prod(m_arr(jj+1:end))-1)'+~excludeL(2)*prod(m_arr(jj+1:end))*d_ii^2;
                r_idcs = r_idcs + d_ii*prod(m_arr(jj:end))*(0:prod(m_arr(1:jj-1))-1);
                Ccell{3} = sparse(r_idcs(:),c_idcs,1,m_tot,d_ii);
            else
                Ccell{3} = sparse(m_tot,d_ii);
            end
            Zjj.C = Ccell;
            Z_ii{jj} = Zjj;
        end
        % Store the operators as a diagonal
        Zop.ops{ii,ii} = Z_ii;
        Pdim = Pdim + prod(m_arr);
        Pdim_arr(ii) = prod(m_arr);
    end
end

% Declare the positive semidefinite matrix variable
Pmat = num2cell(zeros(nZ,nZ));
for p=psatz
    % Declare the monomial basis
    Zs = cell(nZ,1);
    for i=1:nZ
        deg2_i = pdegs{i}(:,1);
        Zs{i} = var1(1).^(0:max(deg2_i(1)-p,0));
        for k=2:N
            Zs{i} = kron(Zs{i},var1(k).^(0:max(deg2_i(k)-p,0))');
        end
    end
    % Declare an SOS matrix in terms of this basis
    [prog,Pmat_j] = sosquadvar(prog,Zs,Zs,Pdim_arr,Pdim_arr,'pos');
    for k=1:numel(Pmat)
        Pmat{k} = Pmat{k} + ((var1(1)-dom(1))*(dom(2)-var1(1))).^p*Pmat_j{k};
    end
end

end