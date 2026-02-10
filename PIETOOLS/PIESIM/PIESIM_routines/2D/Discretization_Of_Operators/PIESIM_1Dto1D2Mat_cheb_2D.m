%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_1Dto1D2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_2PDEstate - discrete matrix representations of 1D to 1D
% cross-state operators (Rxy, Ryx)
% Discretized FULL OPERATOR of 1D-to-1D states
% Called by 'PIESIM_fullPI2Mat_cheb_2D.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial block corresponding to Rxy or Ryx operator
% var1: 2x1 pvar array specifying the primary spatial variables (x,y)
% px - vector with the degrees of differentiability for the 1D states dependent only on x direction 
% py - vector with the degrees of differentiability for the 1D states dependent only on y direction 
% dir - direction of integration ('y' for Rxy and 'x' for Ryx)
%
%
% Outputs:
% A - Chebyshev
% discretizaiton of Rxy or Ryx block that represents a block of a
% total matrix operator for time-advancement of the
% spatially-discretized PIE solution 
% A_2PDEstate - Chebyshev
% discretizaiton of Rxy or Ryx block that represents a block of a
%  total matrix operator for reconstruction of the primary (PDE)
% solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Remove hard-coded variables. Instead, pass variables
%                   defining R as additional inputs.
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto1D2Mat_cheb_2D.m
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_PI2Mat_cheb_opint_discretize_1to1_line.m

function [A, A_nonsquare]=PIESIM_PI2Mat_cheb_opint_discretize_1to1_line(N, R, var1, px, py, dir)
========
=======
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto1D2Mat_cheb_2D.m
% YP, 2/1/2026: rewrote for square operators, truncation made on assembled
% operators

function [A, A_2PDEstate]=PIESIM_1Dto1D2Mat_cheb_2D(N, R, var1, prow, pcol, dir)
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto1D2Mat_cheb_2D.m
>>>>>>>> Stashed changes:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto1D2Mat_cheb_2D.m
=======
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto1D2Mat_cheb_2D.m

% Extract spatial variables (x,y) defining the operator R                   % DJ, 12/16/2024
s1 = var1(1); s2 = var1(2);

if dir=='x'
    var=s2;
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto1D2Mat_cheb_2D.m
<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_PI2Mat_cheb_opint_discretize_1to1_line.m
    pint=px;
    prow=py;
    snstr=s2.varname{1};
else
    var=s1;
    pint=py;
    prow=px;
    snstr=s1.varname{1};
end

    ns_row=length(prow);
    ns_col=length(pint);
========
=======
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto1D2Mat_cheb_2D.m
    snstr=s2.varname{1};
    Nrow=N(2);
    Ncol=N(1);
else
    var=s1;
    snstr=s1.varname{1};
    Nrow=N(1);
    Ncol=N(2);
end

    ns_row=size(prow,2);
    ns_col=size(pcol,2);
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto1D2Mat_cheb_2D.m
>>>>>>>> Stashed changes:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto1D2Mat_cheb_2D.m
=======
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto1D2Mat_cheb_2D.m

for i=1:ns_row
    rsize=Nrow-prow(i)+1;
    for j=1:ns_col
        Rloc=R(i,j);
        csize=Ncol-pcol(j)+1;
        int=zeros(Nrow+1,Ncol+1);
    for k=1:Rloc.nterms
        Rstrip=polynomial(Rloc.coeff(k),Rloc.degmat(k,:),Rloc.varname,Rloc.matdim);
            Rint=subs(Rstrip,var,1);
            % Integrate in pint direction over the interval [-1,1]
            op_int=PIESIM_PI2Mat_opint_cheb(Ncol, Rint, 0);
            index = find(strcmp(Rstrip.varname, snstr));
            Rrow=polynomial(1,Rstrip.degmat(1,index),Rstrip.varname(index),Rstrip.matdim);
            % Convert to a matrix operator in prow direction 
            op_row=PIESIM_Poly2Mat_cheb(Nrow, Rrow, 0);
            int=int+kron(op_row,op_int);
    end % k

    % Truncate the total operator to the correct dimensions

% First, PIE-to-PDE state matrix
     int_2PDEstate=int(:,1:csize);
% Next, PIE time propagator matrix
     int_propagator=int_2PDEstate(1:rsize,:);

    A_cell{i,j}=int_propagator;
    A_cell_2PDEstate{i,j}=int_2PDEstate;
    end % j
end %i

 A=cell2mat(A_cell);
 A_2PDEstate=cell2mat(A_cell_2PDEstate);







