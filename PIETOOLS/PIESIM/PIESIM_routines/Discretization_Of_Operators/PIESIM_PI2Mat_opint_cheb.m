%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_opint_cheb.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation for the integrative Q1 operator  
% Discretizes a FULL OPERATOR, but does not provide a PIE-2-PDE component
% (since it maps to ODE states)
%
%   (Q1 v) = int_{-1}^{1} Rop(s) v(s) ds,
%
% A(i,k+1) = int_{-1}^{1} Rop(i,j)(s) T_k(s) ds, k = 0,...,N-p(j), in the
% column block of state j, T_k the Chebyshev polynomials; the integrals
% follow from the Chebyshev coefficients of Rop(i,j), computed on the grid
% of deg+2 nodes (deg = Rop.maxdeg). Rop may depend on one variable, of any
% name.
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop - Q1 operator of dimension no x ns
% p - vector of dimension 1xns -
% a "degree of smoothness" structure for PDE: p(j) is the degree of
% differentiability of state j, whose Chebyshev expansion has degree N-p(j)
%
% Outputs:
% A - discretization matrix of the Q1 operator: 
% dimension no x n0(N+1) n1N n2(N-1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_PI2Mat_opint_cheb
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
% Initial coding YP  - 6_30_2022% 

function A=PIESIM_PI2Mat_opint_cheb(N, Rop, p);
% Default variable s, for a polynomial without variables
pvar var s

% Numbers of states (columns) and outputs (rows)
ns=length(p);
no=size(Rop,1);

% Degree of Rop
if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

% Chebyshev coefficients of each entry of Rop, from its values on the grid
% of deg+2 nodes
Reval=zeros(no,ns,deg+2);
acheb=zeros(no,ns,deg+2);

chebgrid=cos(pi*(0:deg+1)/(deg+1));
% Polynomial Rop: values at the nodes in its variable (s if it has none)
if isa(Rop,'polynomial')
for j=1:ns
    for i=1:no
        if (isempty(Rop.varname))
            var=s;
        else
            var=Rop.varname;
        end
    Reval(i,j,:)=(subs(Rop(i,j),var,chebgrid))';
    acheb(i,j,:)=fcht(double(Reval(i,j,:)));
    end 
end 
else
% Constant Rop
for j=1:ns
for i=1:no
    Reval(i,j,:)=Rop(i,j);
    acheb(i,j,:)=fcht(double(Reval(i,j,:)));
end
end
end


  
  % Column block j: integrals of Rop(i,j) against T_k, k = 0,...,N-p(j)
  for j=1:ns
        Norder=N-p(j);
       
  A=zeros(no,Norder);
  
        for k=0:Norder
            for i=1:no
 vecint=0;
 for m=0:deg+1
% Integrate Rop (T_(k+m) and T_|k-m|) polynomials on [-1,1]

index=[k+m,abs(k-m)];
for kk=index
 if (kk==0)
     vecint=vecint+0.5*acheb(i,j,m+1)*(1)^1-0.5*acheb(i,j,m+1)*(-1)^1;
 elseif (kk==1)
     vecint=vecint+0.5*0.25*acheb(i,j,m+1)*((1)^0+(1)^2)-0.5*0.25*acheb(i,j,m+1)*((-1)^0+(-1)^2);
 else
    vecint=vecint+0.5*0.5*acheb(i,j,m+1)/(kk+1)*(1)^(kk+1)-0.5*0.5*acheb(i,j,m+1)/(kk-1)*(1)^(kk-1)...
    -0.5*0.5*acheb(i,j,m+1)/(kk+1)*(-1)^(kk+1)+0.5*0.5*acheb(i,j,m+1)/(kk-1)*(-1)^(kk-1);
end % if
end % for kk=index
end % for m
       
        A(i,k+1)=vecint;
            end
        end
        
        A_cell{j}=double(A);
  end
        
        % Assemble the column blocks
        A=cell2mat(A_cell);

  



