%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI_op_discretize_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of 1D to 2D
% operators (R2x, R2y)
%
% Called by 'PIESIM_1Dto2D2Mat_cheb_2D.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial block (one of the components of the 3PI operator -
% multiplicative or integraitve) of R2x or R2y
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% dir - direction ('x' for R2x and 'y' for R2y)
% flag - determines if R0, R1, or R2 component of 3PI operator is given
% flag = 0 - multiplicative (R0) operator
% flag = [1 0] - integrative R1 operator
% flag = [0 1] - integrative R2 operator
%
%
% Outputs:
% A - Chebyshev
% discretizaiton of R0, R1 or R2 component of the R2x/R2y operator
% (depending on the flag)
% A_nonsquare - Chebyshev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Remove hard-coded variables. Instead, pass variables
%                   defining R as additional inputs
% YP, 2/2/2026: Changed to a component representation and a square matrix construction

function A=PIESIM_3PI_op_discretize_2D(N, R, var1, dir, flag)

% Extract spatial variables (x,y) defining the operator R                   % DJ, 12/16/2024
s1 = var1(1);   s2 = var1(2);

if dir=='x'
    var=s2;
    snstr=s2.varname{1};
    Ncol=N(1);
    Nrow=N(2);
else
    var=s1;
    snstr=s1.varname{1};
    Ncol=N(2);
    Nrow=N(1);
end

A=zeros(prod(N+1),Ncol+1);

    for k=1:R.nterms
        Rstrip=polynomial(R.coeff(k),R.degmat(k,:),R.varname,R.matdim);
            Rcomp=subs(Rstrip,var,1);
            if (flag==0)
            % Multiplicative operator in "dir" direction 
            op_comp=PIESIM_3PI_opmult_discretize(Ncol, Rcomp);
            else 
            % Integrative 3PI operator in "dir" direction 
            Rcomp=Rcomp*flag;
            op_comp=PIESIM_3PI_opint_discretize(Ncol, Rcomp(1), Rcomp(2));
            end
            % Poly2mat in the other direction
            index = find(strcmp(Rstrip.varname, snstr));
            Rpoly=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            op_poly=PIESIM_Poly2Mat_cheb(Nrow, Rpoly, 0);
            if (dir=='x')
            for jcol=1:Ncol+1
            A(:,jcol)=A(:,jcol)+kron(op_poly,op_comp(:,jcol));
            end
            else
            for jcol=1:Ncol+1
            A(:,jcol)=A(:,jcol)+kron(op_comp(:,jcol),op_poly);
            end
            
            end
    end % k