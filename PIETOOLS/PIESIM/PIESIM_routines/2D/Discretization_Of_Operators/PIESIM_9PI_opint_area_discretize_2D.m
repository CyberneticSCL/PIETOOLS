%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_9PI_opint_area_discretize_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of
% double-integrated blocks (R11, R12, R21, R22) of 9PI operator in 2D
% Discretizes operator component of 9PI operator (NOT FULL OPERATOR)
% Called by 'PIESIM_9PI2Mat_cheb_2D.m' 
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial block (R11, R12, R21 or R22) of 9PI operator in 2D
% corresponding to a single solution state
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI2Mat_cheb_opint_discretize_area.m
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m
=======
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI2Mat_cheb_opint_discretize_area.m
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m
% rsize - the number of rows in the resulting A matrix block (number of the PIE state Chebyshev coefficients per dimension in the solution state on the left-hand side of the ODE matrix system corresponding to the discrete block in question)
% csize - the number of columns in the resulting A matrix block (number of the PIE state Chebyshev coefficients per dimension in the right-hand side of the ODE matrix system corresponding to the discrete block in question)
========
>>>>>>>> Stashed changes:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m
=======
========
>>>>>>>> Stashed changes:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% var2: 2x1 pvar array specifying the dummy variables for integration
%           (theta,nu)
% lim  - limits of 2D integraton (for example, [-1 sx;-1 sy]) 
%
% Outputs:
% A - Chebyshev
% discretizaiton of R11, R12, R21 or R22 block that represents a block of a
% total matrix operator - size (N+1)x(N+1) - truncation done on assembled
% operators
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Remove hard-coded variables. Instead, pass variables
%                   defining R as additional inputs.
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI2Mat_cheb_opint_discretize_area.m
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m
=======
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI2Mat_cheb_opint_discretize_area.m
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m

function [Afull, Afull_nonsquare]=PIESIM_9PI2Mat_cheb_opint_discretize_area(N, R, rsize, csize, var1, var2, lim)
========
% YP, 1/31/2026 - changed to square consturction, truncation is done on
% assembled operators

function A=PIESIM_9PI_opint_area_discretize_2D(N, R, var1, var2, lim)
>>>>>>>> Stashed changes:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m
=======
========
% YP, 1/31/2026 - changed to square consturction, truncation is done on
% assembled operators

function A=PIESIM_9PI_opint_area_discretize_2D(N, R, var1, var2, lim)
>>>>>>>> Stashed changes:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI_opint_area_discretize_2D.m

% Extract spatial and dummy variables defining the operator R               % DJ, 12/16/2024
s1 = var1(1);   s1_dum = var2(1);
s2 = var1(2);   s2_dum = var2(2);

if isa(R,'polynomial')
    deg=R.maxdeg;
else
    deg=1;
end

Reval=zeros(deg+2,deg+2,deg+2,deg+2);

chebgrid=cos(pi*(0:deg+1)/(deg+1));

if isa(R,'polynomial')
if ismember(s2_dum.varname{1},R.varname)
    Ry=subs(R,s2_dum,chebgrid);
else
    Ry=R*ones(1,deg+2);
end
if ismember(s2.varname{1},R.varname)
    Rynu=subs(Ry,s2,chebgrid);
else
    Rynu=repmat(Ry,deg+2,1);
end
if ismember(s1_dum.varname{1},R.varname) && ismember(s1.varname{1},R.varname)
    for j=1:deg+2
        for i=1:deg+2
        Rmat=subs(subs(Rynu(i,j),s1_dum,chebgrid),s1,chebgrid);
        Reval(:,:,i,j)=Rmat(:,:);
        end
    end
elseif ismember(s1.varname{1},R.varname)
   for j=1:deg+2
        for i=1:deg+2
        Re=subs(Rynu(i,j),s1,chebgrid);
        Reval(:,:,i,j)=repmat(Re',1,deg+2);
        end
   end
   else
       for j=1:deg+2
        for i=1:deg+2
       Re=subs(Rynu(i,j),s1_dum,chebgrid);
       Reval(:,:,i,j)=repmat(Re,deg+2,1);
        end
       end
   end
else
    Reval=R*ones(deg+2,deg+2,deg+2,deg+2);
end

Reval=double(Reval);
Rcheb=fcgltran4d(Reval);

for q=1:deg+2
    for p=1:deg+2
        A(1:N(1)+1,1:N(1)+1,p,q)=PIESIM_integral_projection1D(N(1), Rcheb(:,:,p,q), lim(1,:));
    end
end

 for k=1:N+1
     for i=1:N+1
         Rycheb(:,:)=A(i,k,:,:);
         A2D(i,k,1:N(2)+1,1:N(2)+1)=PIESIM_integral_projection1D(N(2), Rycheb, lim(2,:));
     end
 end

acheb=permute(A2D,[1 3 2 4]);
A=reshape(acheb,prod(N+1),prod(N+1));




