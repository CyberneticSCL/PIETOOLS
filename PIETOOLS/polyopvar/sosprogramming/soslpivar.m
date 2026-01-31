function [prog,Pmat,Zop] = soslpivar(prog,Z,pdegs,opts)
% [PROG,CMAT,ZOP] = SOSLPIVAR(PROG,Z,PDEGS,OPTS) generates a positive
% semidefinite matrix CMAT of decision variables and a basis of PI
% operators ZOP representing a positive semidefinite PI operator
%   Pop = ZOP' PMAT ZOP
%
% INPUTS
% - prog:   lpiprogram structure representing an LPI optimization program;
% - Z:      m x 1 'polyopvar' object representing the basis of monomials in
%           the state variables on which the operator Pop should act;
% - pdegs:  scalar integer specifying the maximal monomial degree of the
%           independent variables in the definition of the monomial basis
%           operator Zop;
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
% PIETOOLS - soslpivar
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
pvars = polynomial(Z.pvarname);
dom = Z.dom;
N = numel(pvars);
if N>=2
    error("Multivariate case is currently not supported.")
end

% Check that the independent variables match
var1 = prog.vartable;
if any(~ismember(pvars.varname,var1.varname))
    error("Independent variables in the distributed polynomial must match those in the LPI program.")
end
vars = [var1(1:N),var1(N+1:2*N)];
dom = prog.dom;

excludeL = zeros([3*ones(1,N),1]);
psatz = 0;
if nargin>=3 && isfield(opts,'exclude')
    excludeL = opts.exclude;
end
if nargin>=3 && isfield(opts,'psatz')
    psatz = reshape(opts.psatz,1,[]);
end

% For each distributed monomial, generate a basis of operators acting on
% that distributed monomial
Pdim = 0;
d1 = prod(pdegs+1);
%Zop = cell(nZ,1);
Zop = tensopvar();
for ii=1:nZ
    xdeg_ii = xdegs(ii,:);
    dtot = sum(xdeg_ii);
    if dtot==0
        error("Degree-0 monomials are currently not supported.")
    elseif dtot==1
        % Include the basis of 3-PI operators acting on the single state
        % Zd(s) o Zd(t)
        Zii = nopvar();
        Zii.vars = vars;
        Zii.dom = dom;
        Zii.deg = pdegs;
        Ccell = cell(3,1);
        m_arr = [d1^2; d1^3; d1^3].*(~excludeL(:));
        m_tot = sum(m_arr);
        r_shft = 0;
        if ~excludeL(1)
            c_idcs = ones(d1,1);
            r_idcs = (1:d1)' + d1*(0:d1-1)';
            Ccell{1} = sparse(r_idcs+r_shft,c_idcs,1,m_tot,1);
            r_shft = r_shft + d1^2;
        end
        if ~excludeL(2)
            c_idcs = repmat((1:d1)',[d1,1]);
            r_idcs = kron((1:d1)',ones(d1,1)) + d1*(0:d1^2-1)';
            Ccell{2} = sparse(r_idcs+r_shft,c_idcs,1,m_tot,d1);
            r_shft = r_shft + d1^3;
        end
        if ~excludeL(3)
            c_idcs = repmat((1:d1)',[d1,1]);
            r_idcs = kron((1:d1)',ones(d1,1)) + d1*(0:d1^2-1)';
            Ccell{3} = sparse(r_idcs+r_shft,c_idcs,1,m_tot,d1);
        end
        Zii.C = Ccell;
        Zop.ops{ii,ii} = Zii;
        Pdim = Pdim + Zii.dim(1);
    else
        % The distributed monomial involves multiple state variables
        % --> include a basis of 2-PI operators acting on each state
        % variable,
        %   (Z{1}*x1)(s)*(Z{2}*x2)(s)
        Zii = cell(1,dtot);
        m_arr = d1^2*(~excludeL(2)+~excludeL(3))*ones(1,dtot);
        m_tot = prod(m_arr)*d1;
        for jj=1:dtot
            Zjj = nopvar();
            Zjj.vars = vars;
            Zjj.dom = dom;
            Zjj.deg = pdegs;
            Ccell = cell(3,1);
            if ~excludeL(2)
                c_idcs = kron((1:d1)',ones(prod(m_arr(jj+1:end)),1));       % index of monomial in dummy var
                c_idcs = repmat(c_idcs,[d1*prod(m_arr(1:jj-1)),1]);         % account for Kronecker product with other bases
                r_idcs = kron((1:d1)',ones(d1*prod(m_arr(jj+1:end)),1));    % index of monomial in primary var
                r_idcs = r_idcs + d1*(0:d1^2*prod(m_arr(jj+1:end))-1)';     % account for location in vector-valued object, (Im o Zd(s))
                r_idcs = r_idcs + d1*prod(m_arr(jj:end))*(0:prod(m_arr(1:jj-1))-1); % account for Kronecker product with other bases
                Ccell{2} = sparse(r_idcs(:),c_idcs,1,m_tot,d1);
            end
            if ~excludeL(3)
                c_idcs = kron((1:d1)',ones(prod(m_arr(jj+1:end)),1));
                c_idcs = repmat(c_idcs,[d1*prod(m_arr(1:jj-1)),1]);
                r_idcs = kron((1:d1)',ones(d1*prod(m_arr(jj+1:end)),1));
                r_idcs = r_idcs + d1*(0:d1^2*prod(m_arr(jj+1:end))-1)'+~excludeL(2)*prod(m_arr(jj+1:end))*d1^3;
                r_idcs = r_idcs + d1*prod(m_arr(jj:end))*(0:prod(m_arr(1:jj-1))-1);
                Ccell{3} = sparse(r_idcs(:),c_idcs,1,m_tot,d1);
            end
            Zjj.C = Ccell;
            Zii{jj} = Zjj;
        end
        Zop.ops{ii,ii} = Zii;
        Pdim = Pdim + prod(m_arr);
    end
end

% Declare the positive semidefinite matrix variable
Pmat = 0;
for j=psatz
    [prog,Pmat_j] = sosposmatr(prog,Pdim);
    Pmat = Pmat + ((var1(1)-dom(1))*(dom(2)-var1(1))).^j*Pmat_j;
end



end