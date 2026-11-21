%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_fullPI2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs a discrete version of the full 2D PI operator  
% R maps R^n0 x L2^nx[x] x L2^ny[y] x L2^n2[x,y] to the same kind of
% space; A is the block matrix
%   [R00 R0x R0y R02; Rx0 Rxx Rxy Rx2; Ry0 Ryx Ryy Ry2; R20 R2x R2y R22]
% of the discretized blocks, acting on Chebyshev coefficients. Rows of
% states (flag >= 2) and columns of states (flag odd) keep the
% coefficients resolved for the degree of smoothness of the state.
%
% Inputs:
% Rop - 2D PI operator 
%    ('struct' with the kernel fields R00, R0x, R0y, R02, Rx0, Rxx, Rxy,
%    Rx2, Ry0, Ryx, Ryy, Ry2, R20, R2x, R2y, R22 – cells, along each
%    variable on which both spaces depend, of the multiplier and the
%    integrals from the lower bound and to the upper bound – and var1,
%    var2, the 2x1 variables and dummy variables)
% psize - size of the PIE problem: all variables defining the size of the PIE problem
%    ('struct' with fields N = [N_x, N_y], nx, ny, n)
% flag = 0 if a structure maps disturbances or control inputs to regulated or observed outputs (for D11, D12,
% D21, D22 operators)
% flag = 1 if a structure maps solution states (ODE+PDE) to regulated or observed outputs (for C1, C2 operators)
% flag = 2 if a structure maps disturbances or control inputs to solution states (ODE+PDE) (for Tw, Tu, B1 and B2
% operators)
% flag = 3 if a structure maps solution states (ODE+PDE) to solution states (ODE+PDE) (for A and T)
%
% Outputs:
% A - discretization matrix of the time-dependent PIE propagator
% A_2PDEstate -discretization matrix for transform between
% fundamental and primary solution
%    (as A, with all coefficients in the rows of states; returned only
%    if requested)
%
% NOTE: this routine also uses 1D discretization routines 
% (for operators that involve 0D-1D, 1D-0D maps or 1D-1D maps involving a single
% direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_fullPI2Mat_cheb_2D
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
% Initial coding YP  -  4/16/2024
% DJ, 12/16/2024: Explicitly pass the variables Rop.var1 and Rop.var2 on
%                   which Rop is defined to the discretization subroutines,
%                   in case these variables are not equal to (s1,s2) and
%                   (s1_dum,s2_dum);
% YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output A_nonsquare
%   renamed A_2PDEstate; flag 1 now marks state-to-output and flag 2
%   input-to-state operators (swapped); outputs may be infinite-dimensional
%   – every row block is discretized when its kernel is nonempty, output
%   rows with degree of smoothness zero; N is [N_x, N_y] and each 2-D state
%   has its own [p_x; p_y].

function [A, A_2PDEstate]=PIESIM_fullPI2Mat_cheb_2D(Rop, psize, flag)

% Chebyshev orders [N_x, N_y]
N=psize.N;

% Define degree of smoothness p of 2D-1var states
px = repelem(0:length(psize.nx)-1, psize.nx);
py = repelem(0:length(psize.ny)-1, psize.ny);
% Define degree of smoothness p of 2D-2var states
[cols, rows] = meshgrid(0:size(psize.n,2)-1, 0:size(psize.n,1)-1);
p = [repelem(cols(:), psize.n(:))'; repelem(rows(:), psize.n(:))'];

 % Iniatilize the matrix blocks of the full PI operator as empty
 R00block=[];
 Rx0block{2}=[];
 Ry0block{2}=[];
 R20block{2}=[];

 R0xblock=[];
 Rxxblock{2}=[];
 Ryxblock{2}=[];
 R2xblock{2}=[];

 R0yblock=[];
 Rxyblock{2}=[];
 Ryyblock{2}=[];
 R2yblock{2}=[];

 R02block=[];
 Rx2block{2}=[];
 Ry2block{2}=[];
 R22block{2}=[];

% Discretize the first column
% In each block, prow (pcol) is the degree of smoothness of the row
% (column) states, zero for output rows, flag <= 1 (input columns, flag
% even)

if size(Rop.R00,2)>0

% R00block 
 R00block=double(Rop.R00);

% Rx0block 

    if min(size(Rop.Rx0))>0  

    if (flag<=1)
        prow=zeros(1,size(Rop.Rx0,1),1);
    else
        prow=px;
    end

    if (nargout==1)
    Rx0block{1}= PIESIM_Poly2Mat_cheb(N(1), Rop.Rx0, prow);
    else
    [Rx0block{1}, Rx0block{2}]=PIESIM_Poly2Mat_cheb(N(1), Rop.Rx0, prow);
    end
    end   % Rx0 block


% Ry0block 

    if min(size(Rop.Ry0))>0  

    if (flag<=1)
        prow=zeros(1,size(Rop.Ry0,1));
    else
        prow=py;
    end

    if (nargout==1)
    Ry0block{1}= PIESIM_Poly2Mat_cheb(N(2),Rop.Ry0, prow);
    else
    [Ry0block{1}, Ry0block{2}]=PIESIM_Poly2Mat_cheb(N(2), Rop.Ry0, prow);
    end
    end   % Ry0 block

%R20 block

  
     if min(size(Rop.R20))>0  

    if (flag<=1)
        prow=zeros(2,size(Rop.R20,1));
    else
        prow=p;
    end


     if (nargout==1)
     R20block{1}= PIESIM_Poly2Mat_cheb_2D(N, Rop.R20, Rop.var1, prow);
     else
     [R20block{1}, R20block{2}]= PIESIM_Poly2Mat_cheb_2D(N, Rop.R20, Rop.var1, prow);
     end
     end %R20 block

end % size(Rop.R00,2) = first column

% Discretize the second column

if (size(Rop.R0x,2)>0)

    if (mod(flag,2)==0)
    pcol=zeros(1,size(Rop.R0x,2));
    else
    pcol=px;
    end

    % R0x block

    if (size(Rop.R0x,1)>0)
        R0xblock=PIESIM_PI2Mat_opint_cheb(N(1), Rop.R0x, pcol);
    end % R0x block

      % Rxx block
      if min(size(Rop.Rxx{1}))>0 

    if (flag<=1)
        prow=zeros(1,size(Rop.Rxx{1},1));
    else
        prow=px;
    end

      Rn.R0=Rop.Rxx{1};
      Rn.R1=Rop.Rxx{2};
      Rn.R2=Rop.Rxx{3};
     if (nargout==1)
        Rxxblock{1}=PIESIM_3PI2Mat_cheb(N(1), Rn, prow, pcol);
     else
        [Rxxblock{1},Rxxblock{2}]=PIESIM_3PI2Mat_cheb(N(1), Rn, prow, pcol);
     end
     end

      % Ryx block
 if min(size(Rop.Ryx))>0 
         if (flag<=1)
        prow=zeros(1,size(Rop.Ryx,1));
    else
        prow=py;
    end

         if(nargout==1)
 Ryxblock{1}=PIESIM_1Dto1D2Mat_cheb_2D(N, Rop.Ryx, Rop.var1, prow, pcol, 'x');
         else
 [Ryxblock{1}, Ryxblock{2}]=PIESIM_1Dto1D2Mat_cheb_2D(N, Rop.Ryx, Rop.var1, prow, pcol, 'x');
         end
 end

      % R2x block
     if min(size(Rop.R2x{1,1}))>0
     if (flag<=1)
        prow=zeros(2,size(Rop.R2x{1,1},1));
    else
        prow=p;
    end
         if(nargout==1)
     R2xblock{1}=PIESIM_1Dto2D2Mat_cheb_2D(N, Rop.R2x, Rop.var1, prow, pcol, 'x');
         else
     [R2xblock{1},R2xblock{2}]=PIESIM_1Dto2D2Mat_cheb_2D(N, Rop.R2x, Rop.var1, prow, pcol,'x');
    end
 end


end % second column

% Discretize the third column

if (size(Rop.R0y,2)>0)  
    if (mod(flag,2)==0)
    pcol=zeros(1,size(Rop.R0y,2));
    else
    pcol=py;
    end


% R0y block
    if (size(Rop.R0y,1)>0)
        R0yblock=PIESIM_PI2Mat_opint_cheb(N(2), Rop.R0y, pcol);
    end % R0x block

      % Rxy block
     if min(size(Rop.Rxy,1))>0
     if (flag<=1)
        prow=zeros(1,size(Rop.Rxy,1));
    else
        prow=px;
    end
         if(nargout==1)
    [Rxyblock{1}]=PIESIM_1Dto1D2Mat_cheb_2D(N, Rop.Rxy, Rop.var1, prow, pcol, 'y');
         else
    [Rxyblock{1}, Rxyblock{2}]=PIESIM_1Dto1D2Mat_cheb_2D(N, Rop.Rxy, Rop.var1, prow, pcol, 'y');
         end
    end

      % Ryy block
      if min(size(Rop.Ryy{1}))>0
    if (flag<=1)
        prow=zeros(1,size(Rop.Ryy{1},1));
    else
        prow=py;
    end
      Rn.R0=Rop.Ryy{1};
      Rn.R1=Rop.Ryy{2};
      Rn.R2=Rop.Ryy{3};
      if (nargout==1)
      Ryyblock{1}=PIESIM_3PI2Mat_cheb(N(2), Rn, prow, pcol);
      else
      [Ryyblock{1},Ryyblock{2}]=PIESIM_3PI2Mat_cheb(N(2), Rn, prow, pcol);
      end
      end

      % R2y block
    if min(size(Rop.R2y{1}))>0
    if (flag<=1)
        prow=zeros(2,size(Rop.R2y{1},1));
    else
        prow=p;
    end
         if(nargout==1)
     R2yblock{1}=PIESIM_1Dto2D2Mat_cheb_2D(N, Rop.R2y, Rop.var1, prow, pcol, 'y');
         else
     [R2yblock{1},R2yblock{2}]=PIESIM_1Dto2D2Mat_cheb_2D(N, Rop.R2y, Rop.var1, prow, pcol,'y');
    end
    end

end % third column


% Discretize the fourth column

if (size(Rop.R02,2)>0)
    if (mod(flag,2)==0)
    pcol=zeros(2,size(Rop.R02,2));
    else
    pcol=p;
    end

    % R02 block
    if (size(Rop.R02,1)>0)
    R02block=PIESIM_2Dto0D2Mat_cheb_2D(N, Rop.R02, Rop.var1, pcol);
    end

    % Rx2 block 
    if min(size(Rop.Rx2{1}))>0
    if (flag<=1)
        prow=zeros(1,size(Rop.Rx2{1},1));
    else
        prow=px;
    end
    if(nargout==1)
     Rx2block{1}=PIESIM_2Dto1D2Mat_cheb_2D(N, Rop.Rx2, Rop.var1, prow, pcol, 'x');
         else
     [Rx2block{1},Rx2block{2}]=PIESIM_2Dto1D2Mat_cheb_2D(N, Rop.Rx2, Rop.var1, prow, pcol,'x');
    end
    end

    % Ry2 block 
    if min(size(Rop.Ry2{1}))>0 
    if (flag<=1)
        prow=zeros(1,size(Rop.Ry2{1},1));
    else
        prow=py;
    end
        if(nargout==1)
     Ry2block{1}=PIESIM_2Dto1D2Mat_cheb_2D(N, Rop.Ry2, Rop.var1, prow, pcol, 'y');
         else
     [Ry2block{1},Ry2block{2}]=PIESIM_2Dto1D2Mat_cheb_2D(N, Rop.Ry2, Rop.var1, prow, pcol,'y');
    end
    end

    % R22 block
     if min(size(Rop.R22{1,1}))>0   
    if (flag<=1)
        prow=zeros(2,size(Rop.R22{1,1},1));
    else
        prow=p;
    end
     if (nargout==1) 
     R22block{1}=PIESIM_9PI2Mat_cheb_2D(N, Rop.R22, Rop.var1, Rop.var2, prow, pcol);
     else
     [R22block{1}, R22block{2}]=PIESIM_9PI2Mat_cheb_2D(N, Rop.R22, Rop.var1, Rop.var2, prow, pcol);
     end
         end

end % fourth column

% Compose all operators

A=[];
A_2PDEstate=[];

     A=[R00block R0xblock R0yblock R02block; ...
        Rx0block{1} Rxxblock{1} Rxyblock{1} Rx2block{1};...
        Ry0block{1} Ryxblock{1} Ryyblock{1} Ry2block{1};...
        R20block{1} R2xblock{1} R2yblock{1} R22block{1}];
     if (nargout>1)
        A_2PDEstate=[R00block R0xblock R0yblock R02block; ...
        Rx0block{2} Rxxblock{2} Rxyblock{2} Rx2block{2};...
        Ry0block{2} Ryxblock{2} Ryyblock{2} Ry2block{2};...
        R20block{2} R2xblock{2} R2yblock{2} R22block{2}];
     end

