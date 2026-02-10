%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI_opint_discretize.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - a dicrete matrix representation of an integrative portion
% of the 3PI operator
% Discretizes an operator component (NOT FULL OPERATOR) 
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R1 - integrative operator R1 of size 1x1
% R2 - integrative operator R2 of size 1x1
% p - scalar - a "degree of smoothness" for the function on which
% R1, R2 operators act
%
% Output:
% A - a  matrix of size (N+1)x(N+1) that represents a Chebyshev
% discretizaiton of the integrative portion of the PI operator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_28_2022
% DJ, 12/28/2024: Remove redundant variables s, theta;
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI2Mat_cheb_opint_discretize.m
function A=PIESIM_3PI2Mat_cheb_opint_discretize(N, R1, R2, p)

Norder=N-p;

A(1:N+1,1:Norder+1)=0;
========
% YP, 2/1/2026: changed to square operator, truncation done on assembled
% operator
function A=PIESIM_3PI_opint_discretize(N, R1, R2)

A(1:N+1,1:N+1)=0;
>>>>>>>> Stashed changes:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
========
% YP, 2/1/2026: changed to square operator, truncation done on assembled
% operator
function A=PIESIM_3PI_opint_discretize(N, R1, R2)

A(1:N+1,1:N+1)=0;
>>>>>>>> Stashed changes:PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
========
% YP, 2/1/2026: changed to square operator, truncation done on assembled
% operator
function A=PIESIM_3PI_opint_discretize(N, R1, R2)

A(1:N+1,1:N+1)=0;
>>>>>>>> Stashed changes:PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m

if isa(R1,'polynomial')
    deg1=R1.maxdeg;
else
    deg1=1;
end

if isa(R2,'polynomial')
    deg2=R2.maxdeg;
else
    deg2=1;
end

maxdeg=max(deg1,deg2);

<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI2Mat_cheb_opint_discretize.m
Rcheb1=PIESIM_PIop_discretize(R1,maxdeg);
Rcheb2=PIESIM_PIop_discretize(R2,maxdeg);

for k=0:Norder
========
========
>>>>>>>> Stashed changes:PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
========
>>>>>>>> Stashed changes:PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
Rcheb1=PIESIM_PI_discretize(R1,maxdeg);
Rcheb2=PIESIM_PI_discretize(R2,maxdeg);

for k=0:N
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
>>>>>>>> Stashed changes:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
========
>>>>>>>> Stashed changes:PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
========
>>>>>>>> Stashed changes:PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI_opint_discretize.m
    % Finding a column vector of Chebyshev coefficient for the column k+1
    % of matrix A
vecmult1=zeros(N+1,1);
vecmult2=zeros(N+1,1);
for j=0:maxdeg+1
    for m=0:maxdeg+1
% Integrate R1 (T_(k+m) and T_|k-m|) polynomials on [-1,s]
% Integrate R2 (T_(k+m) and T_|k-m|) polynomials on [s,1]

vecint1=zeros(N+1,1);
vecint2=zeros(N+1,1);

i=[k+m,abs(k-m)];
for kk=i
 if (kk==0)
     vecint1(2)=vecint1(2)+0.5*Rcheb1(j+1,m+1);
     vecint2(2)=vecint2(2)-0.5*Rcheb2(j+1,m+1);
     vecint1(1)=vecint1(1)-0.5*Rcheb1(j+1,m+1)*(-1)^1;
     vecint2(1)=vecint2(1)+0.5*Rcheb2(j+1,m+1)*(1)^1;

  elseif (kk==1)
     vecint1(1)=vecint1(1)+0.5*0.25*Rcheb1(j+1,m+1);
     vecint1(3)=vecint1(3)+0.5*0.25*Rcheb1(j+1,m+1);
     vecint2(1)=vecint2(1)-0.5*0.25*Rcheb2(j+1,m+1);
     vecint2(3)=vecint2(3)-0.5*0.25*Rcheb2(j+1,m+1);
     vecint1(1)=vecint1(1)-0.5*0.25*Rcheb1(j+1,m+1)*((-1)^0+(-1)^2);
     vecint2(1)=vecint2(1)+0.5*0.25*Rcheb2(j+1,m+1)*((1)^0+(1)^2);
else
    if (kk+1+1<=N+1)
    vecint1(kk+1+1)=vecint1(kk+1+1)+0.5*0.5*Rcheb1(j+1,m+1)/(kk+1);
    vecint2(kk+1+1)=vecint2(kk+1+1)-0.5*0.5*Rcheb2(j+1,m+1)/(kk+1);
    end
    if (kk-1+1<=N+1)
    vecint1(kk-1+1)=vecint1(kk-1+1)-0.5*0.5*Rcheb1(j+1,m+1)/(kk-1);
    vecint2(kk-1+1)=vecint2(kk-1+1)+0.5*0.5*Rcheb2(j+1,m+1)/(kk-1);
    end
    vecint1(1)=vecint1(1)-0.5*0.5*Rcheb1(j+1,m+1)/(kk+1)*(-1)^(kk+1)+0.5*0.5*Rcheb1(j+1,m+1)/(kk-1)*(-1)^(kk-1);
    vecint2(1)=vecint2(1)+0.5*0.5*Rcheb2(j+1,m+1)/(kk+1)*(1)^(kk+1)-0.5*0.5*Rcheb2(j+1,m+1)/(kk-1)*(1)^(kk-1);
end % if
end % for kk=i

% Multiply by T_j(s)
   if (j==0) 
       vecmult1=vecmult1+vecint1;
       vecmult2=vecmult2+vecint2;
   else 
       for i=1:N+1
           iorder=i-1;
       if (iorder+j+1<=N+1)
       vecmult1(iorder+j+1)=vecmult1(iorder+j+1)+0.5*vecint1(i);
       vecmult2(iorder+j+1)=vecmult2(iorder+j+1)+0.5*vecint2(i);
       end
       if (abs(iorder-j)+1<=N+1)
       vecmult1(abs(iorder-j)+1)=vecmult1(abs(iorder-j)+1)+0.5*vecint1(i);
       vecmult2(abs(iorder-j)+1)=vecmult2(abs(iorder-j)+1)+0.5*vecint2(i);
       end
       end % i
    end % j
    end % m
end % j
A(:,k+1)=vecmult1(:)+vecmult2(:);
end % k
