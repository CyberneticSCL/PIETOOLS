 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_cheb_opmixed_discretize_2to1.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of 2D to 1D
% operators (Rx2, Ry2)
%
% Called by 'PIESIM_2Dto1D2Mat_cheb.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial block (one of the components of the 3PI operator -
% multiplicative or integraitve) of Rx2 or Ry2
% p - vector with the degrees of differentiability for the 2D states
% p1 - vector with the degrees of differentiability for the 1D states
% (dependent on x for Rx2 and on y for Ry2)
% dir - direction ('x' for Rx2 and 'y' for Ry2)
% flag - determines if R0, R1, or R2 component of 3PI operator is given
% flag = 0 - multiplicative (R0) operator
% flag = [1 0] - integrative R1 operator
% flag = [0 1] - integrative R2 operator
%
%
% Outputs:
% A - Chebyshev discretizaiton of R0, R1 or R2 component of the Rx2/Ry2 block that represents a block of a
% square total matrix operator for time-advancement of the
% spatially-discretized PIE solution (square ODE system matrix)
% A_nonsquare - Chebyshev
% discretizaiton of R0, R1 or R2 component of the Rx2/Ry2 block that represents a block of a
% nonsquare total matrix operator for reconstruction of the primary (PDE)
% solution (nonsquare transformation matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
function [A, A_nonsquare]=PIESIM_PI2Mat_cheb_opmixed_discretize_2to1(N, R, p, p1, dir, flag)
pvar s1 s2 var


if dir=='x'
    var=s2;
    snstr='s2';
else
    var=s1;
    snstr='s1';
end
    ns_row=size(p1,2);
    ns_col=size(p,2);

for i=1:ns_row
    rsize=N-p1(i)+1;
    for j=1:ns_col
        Rloc=R(i,j);
        csize=N-p(j)+1;
        int=zeros(rsize,csize^2);
        int_nonsquare=zeros(N+1,csize^2);
    for k=1:Rloc.nterms
        Rstrip=polynomial(Rloc.coeff(k),Rloc.degmat(k,:),Rloc.varname,Rloc.matdim);
            R_mult=subs(Rstrip,var,1);
            if (flag==0)
            % Multiplicative operator in "dir" direction 
            [op_mult,op_mult_nonsquare]=PIESIM_PI2Mat_cheb_opmult_discretize(N, rsize, R_mult, p(i));
            else 
            % Integrative 3PI operator in "dir" direction 
            R_mult=R_mult*flag;
            op_mult_nonsquare=PIESIM_3PI2Mat_cheb_opint_discretize(N, R_mult(1), R_mult(2), p(i));
            op_mult=op_mult_nonsquare(1:rsize,:);
            end
            % Integrate in the other direction over the interval [-1,1]
            index = find(strcmp(Rstrip.varname, snstr));
            R_int=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            op_int=PIESIM_PI2Mat_cheb_opint_discretize(N, 1, R_int, p(j));
            if (dir=='x')
            for jrow=1:rsize
            int(jrow,:)=int(jrow,:)+kron(op_int,op_mult(jrow,:));
            end
            for jrow=1:N+1
            int_nonsquare(jrow,:)=int_nonsquare(jrow,:)+kron(op_int,op_mult_nonsquare(jrow,:));
            end
            else
            for jrow=1:rsize
            int(jrow,:)=int(jrow,:)+kron(op_mult(jrow,:),op_int);
            end
            for jrow=1:N+1
            int_nonsquare(jrow,:)=int_nonsquare(jrow,:)+kron(op_mult_nonsquare(jrow,:),op_int);
            end
            end % if
    end % k
    A_cell{i,j}=int;
    A_cell_nonsquare{i,j}=int_nonsquare;
    end % j
end %i

 A=cell2mat(A_cell);
 A_nonsquare=cell2mat(A_cell_nonsquare);







