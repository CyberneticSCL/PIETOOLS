%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_integral_projection1D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, Abig - a dicrete matrix representation corresponding to a
% 1D integration (over [-1 sx], [sx 1], [-1 sy] or [sy 1]) of a polynomial operator already discretized in a
% Chebyshev space
%
% Called by 'PIESIM_9PI2Mat_cheb_opint_discretize_area.m' and 'PIESIM_9PI2Mat_cheb_opint_discretize_line.m'  
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rcheb -  Chebyshev discretization of an operator to be integrated
% csize - the number of columns in the resulting A matrix block (number of the PIE state Chebyshev coefficients per dimension in the right-hand side of the ODE matrix system corresponding to the discrete block in question)
% lim  - limits of 1D integraton ([-1 sx], [sx 1], [-1 sy] or [sy 1]) 
%
% Output:
% A - a non-square matrix of size (N+1)xcsize that represents a Chebyshev
% discretizaiton of the integration of the operator representated by Rcheb over lim
% Abig - a square matrix of size (N+1)x(N+1) that represents a Chebyshev
% discretizaiton of the integration of the operator representated by Rcheb over lim
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Correct limits of for loop j=0:deg and m=0:deg;

function [A,Abig]=PIESIM_integral_projection1D(N, Rcheb, csize, lim)
pvar s theta x y nu theta;

if(lim(1)==-1)
    sign=-1;
else
    sign=1;
end

Abig(1:N+1,1:N+1)=0;
deg=size(Rcheb);

for k=0:N
    % Finding a column vector of Chebyshev coefficient for the column k+1
    % of matrix A
vecmult=zeros(N+1,1);
for j=0:deg(1)-1                                                            % DJ, 12/16/2024
    for m=0:deg(2)-1
% Integrate R (T_(k+m) and T_|k-m|) polynomials on [lim1,lim2]

vecint=zeros(N+1,1);

i=[k+m,abs(k-m)];
for kk=i
 if (kk==0)
     vecint(2)=vecint(2)-sign*0.5*Rcheb(j+1,m+1);
     vecint(1)=vecint(1)+sign*0.5*Rcheb(j+1,m+1)*(sign)^1;

  elseif (kk==1)
     vecint(1)=vecint(1)-sign*0.5*0.25*Rcheb(j+1,m+1);
     vecint(3)=vecint(3)-sign*0.5*0.25*Rcheb(j+1,m+1);
     vecint(1)=vecint(1)+sign*0.5*0.25*Rcheb(j+1,m+1)*((sign)^0+(sign)^2);
else
    if (kk+1+1<=N+1)
    vecint(kk+1+1)=vecint(kk+1+1)-sign*0.5*0.5*Rcheb(j+1,m+1)/(kk+1);
    end
    if (kk-1+1<=N+1)
    vecint(kk-1+1)=vecint(kk-1+1)+sign*0.5*0.5*Rcheb(j+1,m+1)/(kk-1);
    end
    vecint(1)=vecint(1)+sign*0.5*0.5*Rcheb(j+1,m+1)/(kk+1)*(sign)^(kk+1)-sign*0.5*0.5*Rcheb(j+1,m+1)/(kk-1)*(sign)^(kk-1);
end % if
end % for kk=i

% Multiply by T_j(s)
   if (j==0) 
       vecmult=vecmult+vecint;
   else 
       for i=1:N+1
           iorder=i-1;
       if (iorder+j+1<=N+1)
       vecmult(iorder+j+1)=vecmult(iorder+j+1)+0.5*vecint(i);
       end
       if (abs(iorder-j)+1<=N+1)
       vecmult(abs(iorder-j)+1)=vecmult(abs(iorder-j)+1)+0.5*vecint(i);
       end
       end % i
    end % j
    end % m
end % j
Abig(:,k+1)=vecmult(:);
end % k

A=Abig(1:N+1,1:csize);

