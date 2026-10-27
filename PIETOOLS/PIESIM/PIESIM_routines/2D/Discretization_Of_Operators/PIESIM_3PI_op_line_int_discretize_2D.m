 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI_opfull_line_int_discretize_2D.m    PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representations of single-component operators
% (line integral of 3PI operators) entering 2D to 1D operators (Rx2, Ry2)
% Dicretizes R0, R1 or R2 components (NOT FULL OPERATOR) depending on the
% `flag' variable
%
% Called by 'PIESIM_2Dto1D2Mat_cheb_2D.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial block (one of the components of the 3PI operator -
% multiplicative or integraitve) of Rx2 or Ry2
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% dir - direction ('x' for Rx2 and 'y' for Ry2)
% flag - determines if R0, R1, or R2 component of 3PI operator is given
% flag = 0 - multiplicative (R0) operator
% flag = [1 0] - integrative R1 operator
% flag = [0 1] - integrative R2 operator
%
%
% Outputs:
% A - Chebyshev discretizaiton of the R0, R1 or R2 component of the 3PI
% operator integrated in the `dir' direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Remove hard-coded variables. Instead, pass variables
%                   defining R as additional inputs.
% YP, 2/2/2026: Changed to a component representation and a square matrix construction

function A=PIESIM_3PI_op_line_int_discretize_2D(N, R, var1,dir, flag)

% Extract spatial variables (x,y) defining the operator R                   % DJ, 12/16/2024
s1 = var1(1);   s2 = var1(2);

if dir=='x'
    var=s2;
    snstr=s2.varname{1};
    Nrow=N(1);
    Ncol=N(2);
else
    var=s1;
    snstr=s1.varname{1};
    Nrow=N(2);
    Ncol=N(1);
end
        A=zeros(Nrow+1,prod(N+1));

    for k=1:R.nterms
        Rstrip=polynomial(R.coeff(k),R.degmat(k,:),R.varname,R.matdim);
            Rcomp=subs(Rstrip,var,1);
            if (flag==0)
            % Multiplicative operator in "dir" direction 
            op_comp=PIESIM_3PI_opmult_discretize(Nrow, Rcomp);
            else 
            % Integrative 3PI operator in "dir" direction 
            Rcomp=Rcomp*flag;
            op_comp=PIESIM_3PI_opint_discretize(Nrow, Rcomp(1), Rcomp(2));
            end
            % Integrate in the other direction over the interval [-1,1]
            index = find(strcmp(Rstrip.varname, snstr));
            Rint=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            op_int=PIESIM_PI2Mat_opint_cheb(Ncol, Rint, 0);
            if (dir=='x')
            for jrow=1:Nrow+1
            A(jrow,:)=A(jrow,:)+kron(op_int,op_comp(jrow,:));
            end
            else
            for jrow=1:Nrow+1
            A(jrow,:)=A(jrow,:)+kron(op_comp(jrow,:),op_int);
            end
            end
    end % k



