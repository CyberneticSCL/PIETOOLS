%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_fullPI2Mat_cheb_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs a discrete version of the full 2D PI operator  
%
% Inputs:
% Rop - 2D PI operator 
% psize - size of the PIE problem: all variables defining the size of the PIE problem
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
%
% NOTE: this routine also uses 1D discretization routines 
% (for operators that involve 0D-1D, 1D-0D maps or 1D-1D maps involving a single
% direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
% Initial coding YP  - 4_16_2024
=======
% Initial coding YP  -  4/16/2024
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
% DJ, 12/16/2024: Explicitly pass the variables Rop.var1 and Rop.var2 on
%                   which Rop is defined to the discretization subroutines,
%                   in case these variables are not equal to (s1,s2) and
%                   (s1_dum,s2_dum);
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
=======
% YP, 1/27/2026: Added suport for infinite-dimensional outputs
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m

function [A, A_2PDEstate]=PIESIM_fullPI2Mat_cheb_2D(Rop, psize, flag)

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
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
     R20block{1}= PIESIM_Poly2Mat_cheb_2D(N, Rop.R20, Rop.var1, p);
     else
     [R20block{1}, R20block{2}]= PIESIM_Poly2Mat_cheb_2D(N, Rop.R20, Rop.var1, p);
=======
     R20block{1}= PIESIM_Poly2Mat_cheb_2D(N, Rop.R20, Rop.var1, prow);
     else
     [R20block{1}, R20block{2}]= PIESIM_Poly2Mat_cheb_2D(N, Rop.R20, Rop.var1, prow);
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
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
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
      if sum(psize.nx)>0 
=======
      if min(size(Rop.Rxx{1}))>0 

    if (flag<=1)
        prow=zeros(1,size(Rop.Rxx{1},1));
    else
        prow=px;
    end

>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
      Rn.R0=Rop.Rxx{1};
      Rn.R1=Rop.Rxx{2};
      Rn.R2=Rop.Rxx{3};
     if (nargout==1)
        Rxxblock{1}=PIESIM_3PI2Mat_cheb(N(1), Rn, prow, pcol);
     else
        [Rxxblock{1},Rxxblock{2}]=PIESIM_3PI2Mat_cheb(N(1), Rn, prow, pcol);
     end
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
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
 Ryxblock{1}=PIESIM_PI2Mat_cheb_opint_discretize_1to1_line(N, Rop.Ryx, Rop.var1, pcol, py, 'x');
         else
 [Ryxblock{1}, Ryxblock{2}]=PIESIM_PI2Mat_cheb_opint_discretize_1to1_line(N, Rop.Ryx, Rop.var1, pcol, py, 'x');
=======
 Ryxblock{1}=PIESIM_1Dto1D2Mat_cheb_2D(N, Rop.Ryx, Rop.var1, prow, pcol, 'x');
         else
 [Ryxblock{1}, Ryxblock{2}]=PIESIM_1Dto1D2Mat_cheb_2D(N, Rop.Ryx, Rop.var1, prow, pcol, 'x');
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
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
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
     R2xblock{1}=PIESIM_1Dto2D2Mat_cheb(N, Rop.R2x, Rop.var1, p, pcol, 'x');
         else
     [R2xblock{1},R2xblock{2}]=PIESIM_1Dto2D2Mat_cheb(N, Rop.R2x, Rop.var1, p, pcol,'x');
=======
     R2xblock{1}=PIESIM_1Dto2D2Mat_cheb_2D(N, Rop.R2x, Rop.var1, prow, pcol, 'x');
         else
     [R2xblock{1},R2xblock{2}]=PIESIM_1Dto2D2Mat_cheb_2D(N, Rop.R2x, Rop.var1, prow, pcol,'x');
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
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
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
    [Rxyblock{1}]=PIESIM_PI2Mat_cheb_opint_discretize_1to1_line(N, Rop.Rxy, Rop.var1, px, pcol, 'y');
         else
    [Rxyblock{1}, Rxyblock{2}]=PIESIM_PI2Mat_cheb_opint_discretize_1to1_line(N, Rop.Rxy, Rop.var1, px, pcol, 'y');
=======
    [Rxyblock{1}]=PIESIM_1Dto1D2Mat_cheb_2D(N, Rop.Rxy, Rop.var1, prow, pcol, 'y');
         else
    [Rxyblock{1}, Rxyblock{2}]=PIESIM_1Dto1D2Mat_cheb_2D(N, Rop.Rxy, Rop.var1, prow, pcol, 'y');
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
         end
    end

      % Ryy block
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
      if sum(psize.ny)>0
=======
      if min(size(Rop.Ryy{1}))>0
    if (flag<=1)
        prow=zeros(1,size(Rop.Ryy{1},1));
    else
        prow=py;
    end
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
      Rn.R0=Rop.Ryy{1};
      Rn.R1=Rop.Ryy{2};
      Rn.R2=Rop.Ryy{3};
      if (nargout==1)
      Ryyblock{1}=PIESIM_3PI2Mat_cheb(N(2), Rn, prow, pcol);
      else
      [Ryyblock{1},Ryyblock{2}]=PIESIM_3PI2Mat_cheb(N(2), Rn, prow, pcol);
      end
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
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
     R2yblock{1}=PIESIM_1Dto2D2Mat_cheb(N, Rop.R2y, Rop.var1, p, pcol, 'y');
         else
     [R2yblock{1},R2yblock{2}]=PIESIM_1Dto2D2Mat_cheb(N, Rop.R2y, Rop.var1, p, pcol,'y');
=======
     R2yblock{1}=PIESIM_1Dto2D2Mat_cheb_2D(N, Rop.R2y, Rop.var1, prow, pcol, 'y');
         else
     [R2yblock{1},R2yblock{2}]=PIESIM_1Dto2D2Mat_cheb_2D(N, Rop.R2y, Rop.var1, prow, pcol,'y');
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
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
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
    R02block=PIESIM_PI2Mat_cheb_opint_discretize_2to0(N, Rop.R02, Rop.var1, pcol);
=======
    R02block=PIESIM_2Dto0D2Mat_cheb_2D(N, Rop.R02, Rop.var1, pcol);
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
    end

    % Rx2 block 
    if min(size(Rop.Rx2{1}))>0
    if (flag<=1)
        prow=zeros(1,size(Rop.Rx2{1},1));
    else
        prow=px;
    end
    if(nargout==1)
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
     Rx2block{1}=PIESIM_2Dto1D2Mat_cheb(N, Rop.Rx2, Rop.var1, pcol, px, 'x');
         else
     [Rx2block{1},Rx2block{2}]=PIESIM_2Dto1D2Mat_cheb(N, Rop.Rx2, Rop.var1, pcol, px,'x');
=======
     Rx2block{1}=PIESIM_2Dto1D2Mat_cheb_2D(N, Rop.Rx2, Rop.var1, prow, pcol, 'x');
         else
     [Rx2block{1},Rx2block{2}]=PIESIM_2Dto1D2Mat_cheb_2D(N, Rop.Rx2, Rop.var1, prow, pcol,'x');
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
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
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
     Ry2block{1}=PIESIM_2Dto1D2Mat_cheb(N, Rop.Ry2, Rop.var1, pcol, py, 'y');
         else
     [Ry2block{1},Ry2block{2}]=PIESIM_2Dto1D2Mat_cheb(N, Rop.Ry2, Rop.var1, pcol, py,'y');
=======
     Ry2block{1}=PIESIM_2Dto1D2Mat_cheb_2D(N, Rop.Ry2, Rop.var1, prow, pcol, 'y');
         else
     [Ry2block{1},Ry2block{2}]=PIESIM_2Dto1D2Mat_cheb_2D(N, Rop.Ry2, Rop.var1, prow, pcol,'y');
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
    end
    end

    % R22 block
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
     if (nargout==1) 
     R22block{1}=PIESIM_9PI2Mat_cheb_2D(N, Rop.R22, Rop.var1, Rop.var2, p, pcol);
     else
     [R22block{1}, R22block{2}]=PIESIM_9PI2Mat_cheb_2D(N, Rop.R22, Rop.var1, Rop.var2, p, pcol);
     end

=======
     if min(size(Rop.R22{1,1}))>0   
    if (flag<=1)
        prow=zeros(2,size(Rop.R22{1,1},1));
    else
        prow=p;
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
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

