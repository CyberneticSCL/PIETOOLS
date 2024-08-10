%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_cheb_opmixed_discretize_1to2.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of 1D to 2D
% operators (R2x, R2y)
%
% Called by 'PIESIM_1Dto2D2Mat_cheb.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial block (one of the components of the 3PI operator -
% multiplicative or integraitve) of R2x or R2y
% p - vector with the degrees of differentiability for the 2D states
% p1 - vector with the degrees of differentiability for the 1D states
% (dependent on x for R2x and on y for R2y)
% dir - direction ('x' for R2x and 'y' for R2y)
% flag - determines if R0, R1, or R2 component of 3PI operator is given
% flag = 0 - multiplicative (R0) operator
% flag = [1 0] - integrative R1 operator
% flag = [0 1] - integrative R2 operator
%
%
% Outputs:
% A - Chebyshev
% discretizaiton of R0, R1 or R2 component of the R2x/R2y block that represents a block of a
% square total matrix operator for time-advancement of the
% spatially-discretized PIE solution (square ODE system matrix)
% A_nonsquare - Chebyshev
% discretizaiton of R0, R1 or R2 component of the R2x/R2y block that represents a block of a
% nonsquare total matrix operator for reconstruction of the primary (PDE)
% solution (nonsquare transformation matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024

function [A, A_nonsquare]=PIESIM_PI2Mat_cheb_opmixed_discretize_1to2(N, R, p, p1, dir, flag)
pvar s1 s2 var


if dir=='x'
    var=s2;
    snstr='s2';
else
    var=s1;
    snstr='s1';
end
    ns_row=size(p,2);
    ns_col=size(p1,2);

for i=1:ns_row
    rsize=N-p(i)+1;
    for j=1:ns_col
        Rloc=R(i,j);
        csize=N-p1(j)+1;
        sumop=zeros(rsize^2,csize);
        sumop_nonsquare=zeros((N+1)^2,csize);
    for k=1:Rloc.nterms
        Rstrip=polynomial(Rloc.coeff(k),Rloc.degmat(k,:),Rloc.varname,Rloc.matdim);
            R_mult=subs(Rstrip,var,1);
            if (flag==0)
            % Multiplicative operator in "dir" direction 
            [op_mult,op_mult_nonsquare]=PIESIM_PI2Mat_cheb_opmult_discretize(N, rsize, R_mult, p1(j));
            else 
            % Integrative 3PI operator in "dir" direction 
            R_mult=R_mult*flag;
            op_mult_nonsquare=PIESIM_3PI2Mat_cheb_opint_discretize(N, R_mult(1), R_mult(2), p1(j));
            op_mult=op_mult_nonsquare(1:rsize,:);
            end
            % Poly2mat in the other direction
            index = find(strcmp(Rstrip.varname, snstr));
            R_poly=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            [op_poly,op_poly_nonsquare]=PIESIM_Poly2Mat_cheb(N, 1, R_poly, p(i));
            if (dir=='x')
            for jcol=1:csize
            sumop(:,jcol)=sumop(:,jcol)+kron(op_poly,op_mult(:,jcol));
            sumop_nonsquare(:,jcol)=sumop_nonsquare(:,jcol)+kron(op_poly_nonsquare,op_mult_nonsquare(:,jcol));
            end
            else
            for jcol=1:csize
            sumop(:,jcol)=sumop(:,jcol)+kron(op_mult(:,jcol),op_poly);
            sumop_nonsquare(:,jcol)=sumop_nonsquare(:,jcol)+kron(op_mult_nonsquare(:,jcol),op_poly_nonsquare);
            end
            
            end
    end % k
    A_cell{i,j}=sumop;
    A_cell_nonsquare{i,j}=sumop_nonsquare;
    end % j
end %i

 A=cell2mat(A_cell);
 A_nonsquare=cell2mat(A_cell_nonsquare);







