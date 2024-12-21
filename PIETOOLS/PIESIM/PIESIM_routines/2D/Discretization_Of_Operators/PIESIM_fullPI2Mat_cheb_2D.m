%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_fullPI2Mat_cheb_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs a discrete version of the full 2D PI operator  
%
% Inputs:
% Rop - 2D PI operator 
% psize - size of the PIE problem: all variables defining the size of the PIE problem
% flag = 0 if a structure has only the first row and is acting on disturbances or control inputs (for D11, D12,
% D21, D22 operators)
% flag = 1 if a structure is a full operator acting on disturbances or control inputs (for Tw, Tu, B1 and B2
% operators)
% flag = 2 if a structure has only the first row and acts on the PDE + ODE states (for C1, C2 operators)
% flag = 3 if a structure is a full operator acting on the PDE+ODE states (for A and T)
%
% Outputs:
% A - square discretization matrix that transforms spatio-temporal PIE into
% a square system of temporal ODEs 
% A_nonsquare - nonsquare discretization matrix for transform between
% PDE and PIE solution states
% NOTE: nonsquare matrix is only required for Mcheb matrix (discrete form
% of T operator)
%
% NOTE: this routine also uses 1D discretization routines 
% (for operators that involve 0D-1D, 1D-0D maps or 1D-1D maps involving a single
% direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 12/16/2024: Explicitly pass the variables Rop.var1 and Rop.var2 on
%                   which Rop is defined to the discretization subroutines,
%                   in case these variables are not equal to (s1,s2) and
%                   (s1_dum,s2_dum);

function [A, A_nonsquare]=PIESIM_fullPI2Mat_cheb_2D(Rop, psize, flag)

N=psize.N;

% Define degree of smoothness p for 2D-1var components

psize_aux1=[1 psize.nx];
nsum=cumsum(psize.nx);
nsump1=cumsum(psize_aux1);
for i=1:length(psize.nx);
px(nsump1(i):nsum(i))=i-1;
end

psize_aux1=[1 psize.ny];
nsum=cumsum(psize.ny);
nsump1=cumsum(psize_aux1);
for i=1:length(psize.ny);
py(nsump1(i):nsum(i))=i-1;
end

% Define degree of smoothness p for 2D-2var components

psize_aux1=[1 psize.n];
nsum=cumsum(psize.n);
nsump1=cumsum(psize_aux1);
for i=1:length(psize.n);
p(nsump1(i):nsum(i))=i-1;
end

 % Iniatilize the atrix blocks of the full PI operator as empty
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

 if (mod(flag,2)==1)

% Rx0block 
    if(sum(psize.nx)>0)
    if (nargout==1)
    Rx0block{1}= PIESIM_Poly2Mat_cheb(N, Rop.Rx0, px);
    else
    [Rx0block{1}, Rx0block{2}]=PIESIM_Poly2Mat_cheb(N, Rop.Rx0, px);
    end
    end   % Rx0 block


% Ry0block 

    if(sum(psize.ny)>0)
    if (nargout==1)
    Ry0block{1}= PIESIM_Poly2Mat_cheb(N,Rop.Ry0, py);
    else
    [Ry0block{1}, Ry0block{2}]=PIESIM_Poly2Mat_cheb(N, Rop.Ry0, py);
    end
    end   % Ry0 block

%R20 block

     if(sum(psize.n)>0)
     if (nargout==1)
     R20block{1}= PIESIM_Poly2Mat_cheb_2D(N, Rop.R20, Rop.var1, p);
     else
     [R20block{1}, R20block{2}]= PIESIM_Poly2Mat_cheb_2D(N, Rop.R20, Rop.var1, p);
     end
     end %R20 block

 end %flag 

end % size(Rop.R00,2) = first column

% Discretize the second column

if (size(Rop.R0x,2)>0)

    if (flag<=1)
    pcol=zeros(size(Rop.R0x,2),1);
    else
    pcol=px;
    end

    % R0x block

    if (size(Rop.R0x,1)>0)
        R0xblock=PIESIM_PI2Mat_cheb_opint_discretize(N, Rop.R0x, pcol);
    end % R0x block

if (mod(flag,2)==1)
      % Rxx block
      Rn.R0=Rop.Rxx{1};
      Rn.R1=Rop.Rxx{2};
      Rn.R2=Rop.Rxx{3};
     if (nargout==1)
        Rxxblock{1}=PIESIM_3PI2Mat_cheb(N, Rn, px, pcol);
     else
        [Rxxblock{1},Rxxblock{2}]=PIESIM_3PI2Mat_cheb(N, Rn, px, pcol);
     end

      % Ryx block
 if sum(psize.ny)>0
         if(nargout==1)
 Ryxblock{1}=PIESIM_PI2Mat_cheb_opint_discretize_1to1_line(N, Rop.Ryx, Rop.var1, pcol, py, 'x');
         else
 [Ryxblock{1}, Ryxblock{2}]=PIESIM_PI2Mat_cheb_opint_discretize_1to1_line(N, Rop.Ryx, Rop.var1, pcol, py, 'x');
         end
 end

      % R2x block
     if sum(psize.n)>0
         if(nargout==1)
     R2xblock{1}=PIESIM_1Dto2D2Mat_cheb(N, Rop.R2x, Rop.var1, p, pcol, 'x');
         else
     [R2xblock{1},R2xblock{2}]=PIESIM_1Dto2D2Mat_cheb(N, Rop.R2x, Rop.var1, p, pcol,'x');
    end
 end

end % flag

end % second column

% Discretize the third column

if (size(Rop.R0y,2)>0)  
    if (flag<=1)
    pcol=zeros(size(Rop.R0y,2),1);
    else
    pcol=py;
    end


% R0y block
    if (size(Rop.R0y,1)>0)
        R0yblock=PIESIM_PI2Mat_cheb_opint_discretize(N, Rop.R0y, pcol);
    end % R0x block

    if (mod(flag,2)==1)
      % Rxy block
     if sum(psize.nx)>0
         if(nargout==1)
    [Rxyblock{1}]=PIESIM_PI2Mat_cheb_opint_discretize_1to1_line(N, Rop.Rxy, Rop.var1, px, pcol, 'y');
         else
    [Rxyblock{1}, Rxyblock{2}]=PIESIM_PI2Mat_cheb_opint_discretize_1to1_line(N, Rop.Rxy, Rop.var1, px, pcol, 'y');
         end
    end

      % Ryy block
      Rn.R0=Rop.Ryy{1};
      Rn.R1=Rop.Ryy{2};
      Rn.R2=Rop.Ryy{3};
      if (nargout==1)
      Ryyblock{1}=PIESIM_3PI2Mat_cheb(N, Rn, py, pcol);
      else
      [Ryyblock{1},Ryyblock{2}]=PIESIM_3PI2Mat_cheb(N, Rn, py, pcol);
      end

      % R2y block
    if sum(psize.n)>0
         if(nargout==1)
     R2yblock{1}=PIESIM_1Dto2D2Mat_cheb(N, Rop.R2y, Rop.var1, p, pcol, 'y');
         else
     [R2yblock{1},R2yblock{2}]=PIESIM_1Dto2D2Mat_cheb(N, Rop.R2y, Rop.var1, p, pcol,'y');
    end
    end

end % flag

end % third column


% Discretize the fourth column

if (size(Rop.R02,2)>0)
    if (flag<=1)
    pcol=zeros(size(Rop.R02,2),1);
    else
    pcol=p;
    end

    % R02 block
    if (size(Rop.R02,1)>0)
    R02block=PIESIM_PI2Mat_cheb_opint_discretize_2to0(N, Rop.R02, Rop.var1, pcol);
    end

    if (mod(flag,2)==1)

    % Rx2 block 
    if sum(psize.nx)>0
    if(nargout==1)
     Rx2block{1}=PIESIM_2Dto1D2Mat_cheb(N, Rop.Rx2, Rop.var1, pcol, px, 'x');
         else
     [Rx2block{1},Rx2block{2}]=PIESIM_2Dto1D2Mat_cheb(N, Rop.Rx2, Rop.var1, pcol, px,'x');
    end
    end

    % Ry2 block 
    if sum(psize.ny)>0 
        if(nargout==1)
     Ry2block{1}=PIESIM_2Dto1D2Mat_cheb(N, Rop.Ry2, Rop.var1, pcol, py, 'y');
         else
     [Ry2block{1},Ry2block{2}]=PIESIM_2Dto1D2Mat_cheb(N, Rop.Ry2, Rop.var1, pcol, py,'y');
    end
    end

    % R22 block
     if (nargout==1) 
     R22block{1}=PIESIM_9PI2Mat_cheb_2D(N, Rop.R22, Rop.var1, Rop.var2, p, pcol);
     else
     [R22block{1}, R22block{2}]=PIESIM_9PI2Mat_cheb_2D(N, Rop.R22, Rop.var1, Rop.var2, p, pcol);
     end

    end

end % fourth column

% Compose all operators

A=[];
A_nonsquare=[];


if (mod(flag,2)==1)
     A=[R00block R0xblock R0yblock R02block; ...
        Rx0block{1} Rxxblock{1} Rxyblock{1} Rx2block{1};...
        Ry0block{1} Ryxblock{1} Ryyblock{1} Ry2block{1};...
        R20block{1} R2xblock{1} R2yblock{1} R22block{1}];
     if (nargout>1)
        A_nonsquare=[R00block R0xblock R0yblock R02block; ...
        Rx0block{2} Rxxblock{2} Rxyblock{2} Rx2block{2};...
        Ry0block{2} Ryxblock{2} Ryyblock{2} Ry2block{2};...
        R20block{2} R2xblock{2} R2yblock{2} R22block{2}];
     end
else 
    A=[R00block R0xblock R0yblock R02block];
end

