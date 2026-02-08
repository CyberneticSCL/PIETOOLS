classdef (InferiorClasses={?polynomial,?dpvar,?nopvar,?ndopvar})intvar
% INTVAR is a class used to represent functionals acting on distributed
% monomials, taking the form
%
% Kop*(x1^{d1} o ... o xp^{dp})
%   = sum_{i=1}^{q} int_{a}^{b} int_{t_k1(i)}^{b} ... int_{t_kd(i)}^{b} 
%   K_{i}(t1,...,td)*x1(t1)*...*x1(t_{d1})*...*xp(t_{d-dp+1})*...*xp(t_{d})
%                                                 dt_{kd(i)} ... dt_{k1(i)} 
% where kj = omat(i,j) for j=1,...,d and i=1,...,q, where
% d = d1 + ... + dp.
%
% CLASS properties
% - params: m x q*n 'polynomial' or 'dpvar' object representing the kernels
%           of the different integrals appearing in the functional.
%           Specifically, for each term i in the functional, elements 
%           params(1:m,(i-1)*n+1:i*n) specify the kernel K_{i};
% - omat:   q x d array of integers, where d is the cumulative degree of
%           the distributed monomial on which the functional acts,
%           specifying the order of the dummy variables in the integral
%           defining each of the q terms of the functional. Specifically, 
%           for each i in {1,...,q}, if omat(i,:) = [k1,k2,...,kd], then
%           this implies a <= t_{k1} <= t_{k2} <= ... <= t_{kd} <= b,
%           and thus term i is defined by the integral
%                  int_{a}^{b} int_{t_{k1}}^{b} ... int_{t_{kd}}^{b};
% - pvarname:   1 x p cellstr array, specifying the names of the dummy
%               the dummy variables used in definition of the kernels in 
%               Kop.params;
% - dom:    1 x 2 array specifying the interval [a,b] over which the
%           functional integrates;
% - matdim: 1 x 2 array specifying the dimensions [m,n] of the kernels;
%
% NOTES
% - The class does not yet support functionals on monomials defined by
%   multivariate states
% - The order of the variables in pvarname is important, as the elements of
%   omat indicate the order in which these variables are integrated, e.g.
%   if pvarname = {'s1','s2'} and omat(i,:) = [2,1] then this means s2<=s1;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - intvar
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
% DJ, 02/03/2026: Initial coding


properties
    params = polynomial(zeros(0,0));
    omat = zeros(0,0);
    pvarname = {''};
    dom = zeros(0,2);
end
properties (Dependent)
    matdim;
    vars;
end

methods
    function Kop = intvar(params,omat,pvarname,dom)
        % KOP = INTVAR(PARAMS,OMAT,PVARNAME,DOM) constructs an 
        % instance of the intvar class with specified values;
        if nargin==0
            return
        end
        Kop.params = params;
        Kop.omat = omat;
        Kop.pvarname = pvarname;
        Kop.dom = dom;
    end

    function matdim = get.matdim(Kop)
        % Extract the dimensions of the matrix-valued functional
        m = size(Kop.params,1);
        q = size(Kop.omat,1);
        n = size(Kop.params,2)/q;
        if n~=round(n)
            n = nan;
        end
        matdim = [m,n];
    end

    function vars = get.vars(Kop)
        % Extract the variables defining the kernel as 'pvar' objects
        vars = polynomial(Kop.pvarname);
    end

end


end