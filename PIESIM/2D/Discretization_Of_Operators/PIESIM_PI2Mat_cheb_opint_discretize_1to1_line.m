%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_cheb_opint_discretize_1to1_line.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of 1D to 1D
% cross-state operators (Rxy, Ryx)
%
% Called by 'PIESIM_fullPI2Mat_cheb_2D.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial block corresponding to Rxy or Ryx operator
% px - vector with the degrees of differentiability for the 1D states dependent only on x direction 
% py - vector with the degrees of differentiability for the 1D states dependent only on y direction 
% dir - direction of integration ('y' for Rxy and 'x' for Ryx)
%
%
% Outputs:
% A - Chebyshev
% discretizaiton of Rxy or Ryx block that represents a block of a
% square total matrix operator for time-advancement of the
% spatially-discretized PIE solution (square ODE system matrix)
% A_nonsquare - Chebyshev
% discretizaiton of Rxy or Ryx block that represents a block of a
% nonsquare total matrix operator for reconstruction of the primary (PDE)
% solution (nonsquare transformation matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024

function [A, A_nonsquare]=PIESIM_PI2Mat_cheb_opint_discretize_1to1_line(N, R, px, py, dir)
pvar s1 s2 var

if dir=='x'
    var=s2;
    pint=px;
    prow=py;
    snstr='s2';
else
    var=s1;
    pint=py;
    prow=px;
    snstr='s1';
end

    ns_row=size(prow,2);
    ns_col=size(pint,2);

for i=1:ns_row
    rsize=N-prow(i)+1;
    for j=1:ns_col
        Rloc=R(i,j);
        csize=N-pint(j)+1;
        int=zeros(rsize,csize);
        int_nonsquare=zeros(N+1,csize);
    for k=1:Rloc.nterms
        Rstrip=polynomial(Rloc.coeff(k),Rloc.degmat(k,:),Rloc.varname,Rloc.matdim);
            Rint=subs(Rstrip,var,1);
            % Integrate in pint direction over the interval [-1,1]
            op_int=PIESIM_PI2Mat_cheb_opint_discretize(N, 1, Rint, pint(j));
            index = find(strcmp(Rstrip.varname, snstr));
            R_row=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            % Convert to a matrix operator in prow direction 
            [op_row,op_row_nonsquare]=PIESIM_Poly2Mat_cheb(N, 1, R_row, prow(i));
            int=int+kron(op_row,op_int);
            int_nonsquare=int_nonsquare+kron(op_row_nonsquare,op_int);
    end % k
    A_cell{i,j}=int;
    A_cell_nonsquare{i,j}=int_nonsquare;
    end % j
end %i

 A=cell2mat(A_cell);
 A_nonsquare=cell2mat(A_cell_nonsquare);







