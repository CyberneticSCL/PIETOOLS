%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_9PI_opint_line_discretize_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representations of
% singly integrated blocks (R01, R02, R10, R20) of 9PI operator in 2D
% Discretizes an operator component of 9PI operator (NOT FULL OPERATOR)
% Called by 'PIESIM_9PI2Mat_cheb_2D.m' 
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial multiplicative block (R01, R02, R10 or R20) of 9PI operator in 2D
% corresponding to a single solution state
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% var2: 2x1 pvar array specifying the dummy variables for integration
%           (theta,nu)
% dir - direction of integration ('x' or 'y')
% sign  - (-1 for [-1,x] or [-1,y] integration and 1 for [x, 1] or [y, 1] integation) 
%
% Outputs:
% A - Chebyshev
% discretizaiton of R01, R02, R10 or R20 block that represents a block of a
% total matrix operator for time-advancement of the
% spatially-discretized PIE solution (spatially-discrertized propagator
% matrix) - consruct square prod(N+1)^2 x prod(N+1)^2 matrices only - truncaiton is done
% on assembled operators
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Remove hard-coded variables. Instead, pass variables
%                   defining R as additional inputs.
% YP, 1/31/2026 - restructured matrix conversion, made it for square
% operators

function A=PIESIM_9PI_opint_line_discretize_2D(N, R, var1, var2, dir, sign)

% Extract spatial and dummy variables defining the operator R               % DJ, 12/16/2024
s1 = var1(1);   s1_dum = var2(1);
s2 = var1(2);   s2_dum = var2(2);

% snative, snstr correspond to the variable not affected by integration (not in the differential or limits) 
% s is variable in limits, theta is the dummy variable of integration (differential)
R_native=false;

if (dir=='x')
    s=s1;
    theta=s1_dum;
    snative=s2;
    snstr=s2.varname{1};
    idir=1;
else
    s=s2;
    theta=s2_dum;
    snative=s1;
    snstr=s1.varname{1};
    idir=2;
end

if isa(R,'polynomial')
    deg=R.maxdeg;
else
    deg=1;
end

chebgrid=cos(pi*(0:deg+1)/(deg+1));
if isa(R,'polynomial')
        check_native=ismember(R.varname,snstr);
        if (any(check_native))
            R_native=true;
            for k=1:R.nterms
            Rstrip=polynomial(R.coeff(k),R.degmat(k,:),R.varname,R.matdim);
            Rstrip=subs(Rstrip,snative,1);
            Reval{k}=subs(subs(Rstrip,theta,chebgrid),s,chebgrid); 
            end
        else
        Reval{1}=subs(subs(R,theta,chebgrid),s,chebgrid); 
        end
    else % not a polymoial 
    Reval{1}=R*ones(deg+2);
    end


for k=1:numel(Reval)
Revald{k}=double(Reval{k});
Rchebcell{k}=fcgltran2d(Revald{k},1);
end

for kcell=1:numel(Reval)
    Rcheb=Rchebcell{kcell};

    Aint=PIESIM_integral_projection1D(N(idir), Rcheb, sign);

% Construct multiplicative operator in the other direction

    Multop=eye(prod(N+1));
    Idx=eye(N(2)+1);
    Idy=eye(N(1)+1);
   
if (R_native)
    index = find(strcmp(R.varname, snstr));
    Rop=polynomial(1,R.degmat(kcell,index),R.varname(index),R.matdim);
    Multop=PIESIM_9PI_opmult_discretize_2D(N, Rop, var1);
end

if (dir=='x')
    if (exist('A','var'))
    A=A+kron(Idx,Aint)*Multop;
   else
    A=kron(Idx,Aint)*Multop;
    end
else
      if (exist('A','var'))
    A=A+kron(Aint,Idy)*Multop;
     else
    A=kron(Aint,Idy)*Multop;
      end
end

end % kcell




