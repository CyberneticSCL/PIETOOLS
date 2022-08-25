function [Pcomp] = mtimes(P1,P2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcomp] = mtimes(P1,P2) takes in two operators, 
% P1: R^m0 x L2^mx x L2^my x L2^m2 to R^k0 x L2^kx x L2^ky x L2^k2
% P2: R^k0 x L2^kx x L2^ky x L2^k2 to R^m0 x L2^mx x L2^my x L2^m2
% It returns the composition 
% P1 o P2:  R^m0 x L2^mx x L2^my x L2^m2 to R^n0 x L2^nx x L2^ny x L2^n2
% Version 1.0
% Date: 02/01/21
% 
% INPUT
% P1, P2: dopvar2d or opvar2d class objects with same inner dimension
%   - YOU CANNOT MULTIPLY A DOPVAR WITH A DOPVAR!
%
% OUTPUT
% Pcomp: the matlab structure which is the equivalent operator of
% P1 o P2, where 'o' represents composition of operators.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - mtimes
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 07_12_2021
% 02/19/2022 - DJ: Update for product of empty opvar with matrix
% 04/06/2022 - DJ: Update to increase speed
% 04/14/2022 - DJ: Update to use "int_simple"


if isa(P1,'dopvar2d') && isa(P2,'dopvar2d')
    error("Two decision opvars cannot be multiplied");
        
elseif (isa(P1,'dopvar2d') && isa(P2,'opvar2d')) || (isa(P1,'opvar2d') && isa(P2,'dopvar2d'))
    %P1.dim = P1.dim; P2.dim = P2.dim;
    if any(P1.dim(:,2)~=P2.dim(:,1))
        error("Composition requires inner dimensions of the operators to match");
    end
    if any(P1.I~=P2.I)
        error('Operators act on different intervals and cannot be composed');
    end
    
    ds = P2.var1(:);
    dt = P2.var2(:);
    I = P1.I;
    a = I(:,1);
    b = I(:,2);
    
    % Initialize the composition
    Pcomp = dopvar2d([],[P1.dim(:,1),P2.dim(:,2)],I,ds,dt);
    
    n0 = Pcomp.dim(1,1);    m0 = Pcomp.dim(1,2);
    nx = Pcomp.dim(2,1);    mx = Pcomp.dim(2,2);
    ny = Pcomp.dim(3,1);    my = Pcomp.dim(3,2);
    n2 = Pcomp.dim(4,1);    m2 = Pcomp.dim(4,2);
    
    % No need to compose if either object is empty
    if (isempty(P1) || isempty(P2)) %|| (P1==0 || P2==0)
        return
    end    
    
    % Distinguish most common case, in which opvar only maps 2D functions
    if all(all(P1.dim(1:3,:)==zeros(3,2))) && all(all(P2.dim(1:3,:)==zeros(3,2)))
        
        R22 = PL2L_compose2D(P1.R22,P2.R22,ds,dt,I);
        Pcomp.R22 = R22;
        
    else
    
    % Compute the composite functions
    % % Note the extensive use of "if" statements,
    % % to avoid use of expensive composition functions
    
    % % % Terms mapping to R % % %
    if n0>0
        
    % % % R --> R component
    if m0>0
    % R -> R -> R
    Pcomp.R00 = P1.R00*P2.R00;
    % R -> L2[x] -> R
    Pcomp.R00 = Pcomp.R00 + int_simple(P1.R0x*P2.Rx0,ds(1),a(1),b(1));
    % R -> L2[y] -> R
    Pcomp.R00 = Pcomp.R00 + int_simple(P1.R0y*P2.Ry0,ds(2),a(2),b(2));
    % R -> L2[x,y] -> R
    Pcomp.R00 = Pcomp.R00 + int_simple(int_simple(P1.R02*P2.R20,ds(1),a(1),b(1)),ds(2),a(2),b(2));
    end
    
    % % % L2[x] --> R component
    if mx>0
    % L2[x] -> R -> R
    Pcomp.R0x = P1.R00*P2.R0x;
    % L2[x] -> L2[x] -> R
    Pcomp.R0x = Pcomp.R0x + PIoMat_composeND(P1.R0x,P2.Rxx,ds,dt,I);
    % L2[x] -> L2[y] -> R
    Pcomp.R0x = Pcomp.R0x + int_simple(P1.R0y*P2.Ryx,ds(2),a(2),b(2));
    % L2[x] -> L2[x,y] -> R
    Pcomp.R0x = Pcomp.R0x + int_simple(PIoMat_composeND(P1.R02,P2.R2x,ds,dt,I),ds(2),a(2),b(2));
    end
    
    % % % L2[y] --> R component
    if my>0
    % L2[y] -> R -> R
    Pcomp.R0y = P1.R00*P2.R0y;
    % L2[y] -> L2[x] -> R
    Pcomp.R0y = Pcomp.R0y + int_simple(P1.R0x*P2.Rxy,ds(1),a(1),b(1));
    % L2[y] -> L2[y] -> R
    Pcomp.R0y = Pcomp.R0y + PIoMat_composeND(P1.R0y,P2.Ryy,ds,dt,I);
    % L2[y] -> L2[x,y] -> R
    Pcomp.R0y = Pcomp.R0y + int_simple(PIoMat_composeND(P1.R02,P2.R2y,ds,dt,I),ds(1),a(1),b(1));
    end
    
    % % % L2[x,y] --> R component
    if m2>0
    % L2[x,y] -> R -> R
    Pcomp.R02 = P1.R00*P2.R02;
    % L2[x,y] -> L2[x] -> R
    Pcomp.R02 = Pcomp.R02 + PIoMat_composeND(P1.R0x,P2.Rx2,ds,dt,I);
    % L2[x,y] -> L2[y] -> R
    Pcomp.R02 = Pcomp.R02 + PIoMat_composeND(P1.R0y,P2.Ry2,ds,dt,I);
    % L2[x,y] -> L2[x,y] -> R
    Pcomp.R02 = Pcomp.R02 + PIoMat_composeND(P1.R02,P2.R22,ds,dt,I);
    end
    
    end
    % % % Components mapping to L2[x] % % %
    if nx~=0
        
    % % % R --> L2[x] component
    if m0>0
    % R -> R -> L2[x]
    Pcomp.Rx0 = P1.Rx0*P2.R00;
    % R -> L2[x] -> L2[x]
    Pcomp.Rx0 = Pcomp.Rx0 + PIoMat_composeND(P1.Rxx,P2.Rx0,ds,dt,I);
    % R -> L2[y] -> L2[x]
    Pcomp.Rx0 = Pcomp.Rx0 + int_simple(P1.Rxy*P2.Ry0,ds(2),a(2),b(2));
    % R -> L2[x,y] -> L2[x]
    Pcomp.Rx0 = Pcomp.Rx0 + int_simple(PIoMat_composeND(P1.Rx2,P2.R20,ds,dt,I),ds(2),a(2),b(2));
    end
    
    % % % L2[x] --> L2[x] component
    if mx>0
    % L2[x] -> R -> L2[x]
    Rxx_00 = P1.Rx0*subs(P2.R0x,ds(1),dt(1));
    % L2[x] -> L2[x] -> L2[x]
    Rxx_xx = PL2L_compose1D(P1.Rxx,P2.Rxx,ds,dt,I,1);
    % L2[x] -> L2[y] -> L2[x]
    Rxx_yy = P1.Rxy*subs(P2.Ryx,ds(1),dt(1));
    % L2[x] -> L2[x,y] -> L2[x]
    Rxx_22 = PL2L_compose1D(P1.Rx2,P2.R2x,ds,dt,I,1);
    
    Pcomp.Rxx{1} = Rxx_xx{1} + int_simple(Rxx_22{1},ds(2),a(2),b(2));
    Pcomp.Rxx{2} = Rxx_00 + Rxx_xx{2} + int_simple(Rxx_yy + Rxx_22{2},ds(2),a(2),b(2));
    Pcomp.Rxx{3} = Rxx_00 + Rxx_xx{3} + int_simple(Rxx_yy + Rxx_22{3},ds(2),a(2),b(2));
    end
    
    % % % L2[y] --> L2[x] component
    if my>0
    % L2[y] -> R -> L2[x]
    Pcomp.Rxy = P1.Rx0*P2.R0y;
    % L2[y] -> L2[x] -> L2[x]
    Pcomp.Rxy = Pcomp.Rxy + PIoMat_composeND(P1.Rxx,P2.Rxy,ds,dt,I);
    % L2[y] -> L2[y] -> L2[x]
    Pcomp.Rxy = Pcomp.Rxy + PIoMat_composeND(P1.Rxy,P2.Ryy,ds,dt,I);
    % L2[y] -> L2[x,y] -> L2[x]
    Pcomp.Rxy = Pcomp.Rxy + PIoMat_composeND(P1.Rx2,P2.R2y,ds,dt,I);
    end
    
    % % % L2[x,y] --> L2[x] component
    if m2>0
    % L2[x,y] -> R -> L2[x]
    Rx2_00 = P1.Rx0*subs(P2.R02,ds(1),dt(1));
    % L2[x,y] -> L2[x] -> L2[x]
    Rx2_xx = PL2L_compose1D(P1.Rxx,P2.Rx2,ds,dt,I,1);
    P2Ry2  = {subs(P2.Ry2{1},ds(1),dt(1)),subs(P2.Ry2{2},ds(1),dt(1)),subs(P2.Ry2{3},ds(1),dt(1))};
    % L2[x,y] -> L2[y] -> L2[x]
    Rx2_yy = PIoMat_composeND(P1.Rxy,P2Ry2,ds,dt,I);
    % L2[x,y] -> L2[x,y] -> L2[x]
    Rx2_22_0 = PL2L_compose1D(P1.Rx2,P2.R22(:,1),ds,dt,I,1);
    Rx2_22_1 = PL2L_compose1D(P1.Rx2,P2.R22(:,2),ds,dt,I,1);
    Rx2_22_2 = PL2L_compose1D(P1.Rx2,P2.R22(:,3),ds,dt,I,1);
    Rx2_22 = cell(3,1);
    Rx2_22{1} = Rx2_22_0{1} + int_simple(var_swap(Rx2_22_1{1},ds(2),dt(2)),dt(2),ds(2),b(2))...
                            + int_simple(var_swap(Rx2_22_2{1},ds(2),dt(2)),dt(2),a(2),ds(2));
    Rx2_22{2} = Rx2_22_0{2} + int_simple(var_swap(Rx2_22_1{2},ds(2),dt(2)),dt(2),ds(2),b(2))...
                            + int_simple(var_swap(Rx2_22_2{2},ds(2),dt(2)),dt(2),a(2),ds(2));
    Rx2_22{3} = Rx2_22_0{3} + int_simple(var_swap(Rx2_22_1{3},ds(2),dt(2)),dt(2),ds(2),b(2))...
                            + int_simple(var_swap(Rx2_22_2{3},ds(2),dt(2)),dt(2),a(2),ds(2));
                        
    Pcomp.Rx2{1} = Rx2_xx{1} + Rx2_22{1};
    Pcomp.Rx2{2} = Rx2_00 + Rx2_yy + Rx2_xx{2} + Rx2_22{2};
    Pcomp.Rx2{3} = Rx2_00 + Rx2_yy + Rx2_xx{3} + Rx2_22{3};
    end
    
    end
    % % % Components mapping to L2[y] % % %
    if ny~=0
        
    % % % R --> L2[y] component
    if m0>0
    % R -> R -> L2[y]
    Pcomp.Ry0 = P1.Ry0*P2.R00;
    % R -> L2[x] -> L2[y]
    Pcomp.Ry0 = Pcomp.Ry0 + int_simple(P1.Ryx*P2.Rx0,ds(1),a(1),b(1));
    % R -> L2[y] -> L2[y]
    Pcomp.Ry0 = Pcomp.Ry0 + PIoMat_composeND(P1.Ryy,P2.Ry0,ds,dt,I);
    % R -> L2[x,y] -> L2[y]
    Pcomp.Ry0 = Pcomp.Ry0 + int_simple(PIoMat_composeND(P1.Ry2,P2.R20,ds,dt,I),ds(1),a(1),b(1));
    end
    
    % % % L2[x] --> L2[y] component
    if mx>0
    % L2[x] -> R -> L2[y]
    Pcomp.Ryx = P1.Ry0*P2.R0x;
    % L2[x] -> L2[x] -> L2[y]
    Pcomp.Ryx = Pcomp.Ryx + PIoMat_composeND(P1.Ryx,P2.Rxx,ds,dt,I);
    % L2[x] -> L2[y] -> L2[y]
    Pcomp.Ryx = Pcomp.Ryx + PIoMat_composeND(P1.Ryy,P2.Ryx,ds,dt,I);
    % L2[x] -> L2[x,y] -> L2[y]
    Pcomp.Ryx = Pcomp.Ryx + PIoMat_composeND(P1.Ry2,P2.R2x,ds,dt,I);
    end
    
    % % % L2[y] --> L2[y] component
    if my>0
    % L2[y] -> R -> L2[y]
    Ryy_00 = P1.Ry0*subs(P2.R0y,ds(2),dt(2));
    % L2[y] -> L2[x] -> L2[y]
    Ryy_xx = P1.Ryx*subs(P2.Rxy,ds(2),dt(2));
    % L2[y] -> L2[y] -> L2[y]
    Ryy_yy = PL2L_compose1D(P1.Ryy,P2.Ryy,ds,dt,I,2);    
    % L2[y] -> L2[x,y] -> L2[y]
    Ryy_22 = PL2L_compose1D(P1.Ry2,P2.R2y,ds,dt,I,2);
    
    Pcomp.Ryy{1} = Ryy_yy{1} + int_simple(Ryy_22{1},ds(1),a(1),b(1));
    Pcomp.Ryy{2} = Ryy_00 + Ryy_yy{2} + int_simple(Ryy_xx + Ryy_22{2},ds(1),a(1),b(1));
    Pcomp.Ryy{3} = Ryy_00 + Ryy_yy{3} + int_simple(Ryy_xx + Ryy_22{3},ds(1),a(1),b(1));
    end
    
    % % % L2[x,y] --> L2[y] component
    if m2>0
    % L2[x,y] -> R -> L2[y]
    Ry2_00 = P1.Ry0*subs(P2.R02,ds(2),dt(2));
    % L2[x,y] -> L2[x] -> L2[y]
    P2Rx2  = {subs(P2.Rx2{1},ds(2),dt(2));subs(P2.Rx2{2},ds(2),dt(2));subs(P2.Rx2{3},ds(2),dt(2))};
    Ry2_xx = PIoMat_composeND(P1.Ryx,P2Rx2,ds,dt,I);
    % L2[x,y] -> L2[y] -> L2[y]
    Ry2_yy = PL2L_compose1D(P1.Ryy,P2.Ry2,ds,dt,I,2);
    % L2[x,y] -> L2[x,y] -> L2[y]
    Ry2_22_0 = PL2L_compose1D(P1.Ry2,P2.R22(1,:),ds,dt,I,2);
    Ry2_22_1 = PL2L_compose1D(P1.Ry2,P2.R22(2,:),ds,dt,I,2);
    Ry2_22_2 = PL2L_compose1D(P1.Ry2,P2.R22(3,:),ds,dt,I,2);
    Ry2_22 = cell(1,3);
    Ry2_22{1} = Ry2_22_0{1} + int_simple(var_swap(Ry2_22_1{1},ds(1),dt(1)),dt(1),ds(1),b(1))...
                            + int_simple(var_swap(Ry2_22_2{1},ds(1),dt(1)),dt(1),a(1),ds(1));
    Ry2_22{2} = Ry2_22_0{2} + int_simple(var_swap(Ry2_22_1{2},ds(1),dt(1)),dt(1),ds(1),b(1))...
                            + int_simple(var_swap(Ry2_22_2{2},ds(1),dt(1)),dt(1),a(1),ds(1));
    Ry2_22{3} = Ry2_22_0{3} + int_simple(var_swap(Ry2_22_1{3},ds(1),dt(1)),dt(1),ds(1),b(1))...
                            + int_simple(var_swap(Ry2_22_2{3},ds(1),dt(1)),dt(1),a(1),ds(1));
                        
    Pcomp.Ry2{1} = Ry2_yy{1} + Ry2_22{1};
    Pcomp.Ry2{2} = Ry2_00 + Ry2_xx + Ry2_yy{2} + Ry2_22{2};
    Pcomp.Ry2{3} = Ry2_00 + Ry2_xx + Ry2_yy{3} + Ry2_22{3};
    end
    
    end
    % % % Components mapping to L2[x,y] % % %
    if n2~=0
        
    % % % R --> L2[x,y] component
    if m0>0
    % R -> R -> L2[x,y]
    Pcomp.R20 = P1.R20*P2.R00;
    % R -> L2[x] -> L2[x,y]
    Pcomp.R20 = Pcomp.R20 + PIoMat_composeND(P1.R2x,P2.Rx0,ds,dt,I);
    % R -> L2[y] -> L2[x,y]
    Pcomp.R20 = Pcomp.R20 + PIoMat_composeND(P1.R2y,P2.Ry0,ds,dt,I);
    % R -> L2[x,y] -> L2[x,y]
    Pcomp.R20 = Pcomp.R20 + PIoMat_composeND(P1.R22,P2.R20,ds,dt,I);
    end
    
    % % % L2[x] --> L2[x,y] component
    if mx>0
    % L2[x] -> R -> L2[x,y]
    R2x_00 = P1.R20*subs(P2.R0x,ds(1),dt(1));
    % L2[x] -> L2[x] -> L2[x,y]
    R2x_xx = PL2L_compose1D(P1.R2x,P2.Rxx,ds,dt,I,1);
    % L2[x] -> L2[y] -> L2[x,y]
    R2x_yy = PIoMat_composeND(P1.R2y,subs(P2.Ryx,ds(1),dt(1)),ds,dt,I);
    % L2[x] -> L2[x,y] -> L2[x,y]
    P2_2x = P2.R2x;
    P2_2x{1} = subs(P2_2x{1},ds(2),dt(2));
    P2_2x{2} = subs(P2_2x{2},ds(2),dt(2));
    P2_2x{3} = subs(P2_2x{3},ds(2),dt(2));
    R2x_22_0 = PL2L_compose1D(P1.R22(:,1),P2_2x,ds,dt,I,1);
    R2x_22_1 = PL2L_compose1D(P1.R22(:,2),P2_2x,ds,dt,I,1);
    R2x_22_2 = PL2L_compose1D(P1.R22(:,3),P2_2x,ds,dt,I,1);
    R2x_22 = cell(3,1);
    R2x_22{1} = subs(R2x_22_0{1},dt(2),ds(2))...
                + int_simple(R2x_22_1{1},dt(2),a(2),ds(2))...
                + int_simple(R2x_22_2{1},dt(2),ds(2),b(2));
    R2x_22{2} = subs(R2x_22_0{2},dt(2),ds(2))... 
                + int_simple(R2x_22_1{2},dt(2),a(2),ds(2))...
                + int_simple(R2x_22_2{2},dt(2),ds(2),b(2));
    R2x_22{3} = subs(R2x_22_0{3},dt(2),ds(2))... 
                + int_simple(R2x_22_1{3},dt(2),a(2),ds(2))...
                + int_simple(R2x_22_2{3},dt(2),ds(2),b(2));
            
    Pcomp.R2x{1} = R2x_xx{1} + R2x_22{1};
    Pcomp.R2x{2} = R2x_00 + R2x_yy + R2x_xx{2} + R2x_22{2};
    Pcomp.R2x{3} = R2x_00 + R2x_yy + R2x_xx{3} + R2x_22{3};
    end
    
    % % % L2[y] --> L2[x,y] component
    if my>0
    % L2[y] -> R -> L2[x,y]
    R2y_00 = P1.R20*subs(P2.R0y,ds(2),dt(2));
    % L2[y] -> L2[x] -> L2[x,y]
    R2y_xx = PIoMat_composeND(P1.R2x,subs(P2.Rxy,ds(2),dt(2)),ds,dt,I);
    % L2[y] -> L2[y] -> L2[x,y]
    R2y_yy = PL2L_compose1D(P1.R2y,P2.Ryy,ds,dt,I,2);
    % L2[y] -> L2[x,y] -> L2[x,y]
    P2_2y = P2.R2y;
    P2_2y{1} = subs(P2_2y{1},ds(1),dt(1));
    P2_2y{2} = subs(P2_2y{2},ds(1),dt(1));
    P2_2y{3} = subs(P2_2y{3},ds(1),dt(1));
    R2y_22_0 = PL2L_compose1D(P1.R22(1,:),P2_2y,ds,dt,I,2);
    R2y_22_1 = PL2L_compose1D(P1.R22(2,:),P2_2y,ds,dt,I,2);
    R2y_22_2 = PL2L_compose1D(P1.R22(3,:),P2_2y,ds,dt,I,2);
    R2y_22 = cell(1,3);
    R2y_22{1} = subs(R2y_22_0{1},dt(1),ds(1))...
                + int_simple(R2y_22_1{1},dt(1),a(1),ds(1))...
                + int_simple(R2y_22_2{1},dt(1),ds(1),b(1));
    R2y_22{2} = subs(R2y_22_0{2},dt(1),ds(1))... 
                + int_simple(R2y_22_1{2},dt(1),a(1),ds(1))...
                + int_simple(R2y_22_2{2},dt(1),ds(1),b(1));
    R2y_22{3} = subs(R2y_22_0{3},dt(1),ds(1))... 
                + int_simple(R2y_22_1{3},dt(1),a(1),ds(1))...
                + int_simple(R2y_22_2{3},dt(1),ds(1),b(1));
    
    Pcomp.R2y{1} = R2y_yy{1} + R2y_22{1};
    Pcomp.R2y{2} = R2y_00 + R2y_xx + R2y_yy{2} + R2y_22{2};
    Pcomp.R2y{3} = R2y_00 + R2y_xx + R2y_yy{3} + R2y_22{3};
    end    
    
    % % L2[x,y] --> L2[x,y] component
    if m2>0
    % L2[x,y] -> R -> L2[x,y]
    R22_00 = P1.R20*subs(P2.R02,ds,dt);
    % L2[x,y] -> L2[x] -> L2[x,y]
    P2Rx2  = {subs(P2.Rx2{1},ds(2),dt(2));subs(P2.Rx2{2},ds(2),dt(2));subs(P2.Rx2{3},ds(2),dt(2))};
    R22_xx = PL2L_compose1D(P1.R2x,P2Rx2,ds,dt,I,1);
    % L2[x,y] -> L2[y] -> L2[x,y]
    P2Ry2  = {subs(P2.Ry2{1},ds(1),dt(1)),subs(P2.Ry2{2},ds(1),dt(1)),subs(P2.Ry2{3},ds(1),dt(1))};
    R22_yy = PL2L_compose1D(P1.R2y,P2Ry2,ds,dt,I,2);
    % L2[x,y] -> L2[x,y] -> L2[x,y]
    R22_22 = PL2L_compose2D(P1.R22,P2.R22,ds,dt,I);
    
    Pcomp.R22{1,1} = R22_22{1,1};
    Pcomp.R22{2,1} = R22_yy{1} + R22_22{2,1};
    Pcomp.R22{3,1} = R22_yy{1} + R22_22{3,1};
    Pcomp.R22{1,2} = R22_xx{1} + R22_22{1,2};
    Pcomp.R22{1,3} = R22_xx{1} + R22_22{1,3};
    Pcomp.R22{2,2} = R22_00 + R22_xx{2} + R22_yy{2} + R22_22{2,2};
    Pcomp.R22{3,2} = R22_00 + R22_xx{3} + R22_yy{2} + R22_22{3,2};
    Pcomp.R22{2,3} = R22_00 + R22_xx{2} + R22_yy{3} + R22_22{2,3};
    Pcomp.R22{3,3} = R22_00 + R22_xx{3} + R22_yy{3} + R22_22{3,3};
    end
    
    end
    
    end
    

elseif ~isa(P2,'dopvar2d') && ~isa(P2,'opvar2d') %multiplication of operator times matrix
    Pcomp = P1;
    
    if all(size(P2)==[1,1]) %scalar multiplication
        Pcomp.R00 = P2*P1.R00;
        Pcomp.R0x = P2*P1.R0x;
        Pcomp.R0y = P2*P1.R0y;
        Pcomp.R02 = P2*P1.R02;
        Pcomp.Rx0 = P2*P1.Rx0;
        Pcomp.Rxy = P2*P1.Rxy;
        Pcomp.Ry0 = P2*P1.Ry0;
        Pcomp.Ryx = P2*P1.Ryx;
        Pcomp.R20 = P2*P1.R20;
        for i=1:3
            Pcomp.Rxx{i,1} = P2*P1.Rxx{i,1};
            Pcomp.Rx2{i,1} = P2*P1.Rx2{i,1};
            Pcomp.R2x{i,1} = P2*P1.R2x{i,1};
            
            Pcomp.Ryy{1,i} = P2*P1.Ryy{1,i};
            Pcomp.Ry2{1,i} = P2*P1.Ry2{1,i};
            Pcomp.R2y{1,i} = P2*P1.R2y{1,i};
            
            for j=1:3
                Pcomp.R22{i,j} = P2*P1.R22{i,j};
            end
        end

    else
        if size(P2,1)~=sum(P1.dim(:,2))
            error("Multiplication requires inner dimensions of the operators to match");
        elseif all(all(P1.dim(:,:)==0))
            warning(['Object P in dopvar2d-matrix product P*M has no rows or columns, P.dim(:,:)=0.'...
                        ' The product has ambiguous dimensions; returning the input matrix M.'])
            Pcomp = P2;
            return
        end

        mindx = cumsum(P1.dim(:,2));
        P2_0 = P2(1:mindx(1),:);            P2_x = P2(mindx(1)+1:mindx(2),:);
        P2_y = P2(mindx(2)+1:mindx(3),:);   P2_2 = P2(mindx(3)+1:end,:);
        
        Pcomp.R00 = P1.R00*P2_0;
        Pcomp.R0x = P1.R0x*P2_x;
        Pcomp.R0y = P1.R0y*P2_y;
        Pcomp.R02 = P1.R02*P2_2;
        Pcomp.Rx0 = P1.Rx0*P2_0;
        Pcomp.Rxy = P1.Rxy*P2_y;
        Pcomp.Ry0 = P1.Ry0*P2_0;
        Pcomp.Ryx = P1.Ryx*P2_x;
        Pcomp.R20 = P1.R20*P2_0;
        
        for i=1:3
            Pcomp.Rxx{i,1} = P1.Rxx{i,1}*P2_x;
            Pcomp.Rx2{i,1} = P1.Rx2{i,1}*P2_2;
            Pcomp.R2x{i,1} = P1.R2x{i,1}*P2_x;
            
            Pcomp.Ryy{1,i} = P1.Ryy{1,i}*P2_y;
            Pcomp.Ry2{1,i} = P1.Ry2{1,i}*P2_2;
            Pcomp.R2y{1,i} = P1.R2y{1,i}*P2_y;
            
            for j=1:3
                Pcomp.R22{i,j} = P1.R22{i,j}*P2_2;
            end
        end
    end
    
else %multiplication of matrix times the operator
    Pcomp = P2;
    
    if all(size(P1)==[1,1]) %scalar multiplication
        Pcomp.R00 = P1*P2.R00;
        Pcomp.R0x = P1*P2.R0x;
        Pcomp.R0y = P1*P2.R0y;
        Pcomp.R02 = P1*P2.R02;
        Pcomp.Rx0 = P1*P2.Rx0;
        Pcomp.Rxy = P1*P2.Rxy;
        Pcomp.Ry0 = P1*P2.Ry0;
        Pcomp.Ryx = P1*P2.Ryx;
        Pcomp.R20 = P1*P2.R20;
        for i=1:3
            Pcomp.Rxx{i,1} = P1*P2.Rxx{i,1};
            Pcomp.Rx2{i,1} = P1*P2.Rx2{i,1};
            Pcomp.R2x{i,1} = P1*P2.R2x{i,1};
            
            Pcomp.Ryy{1,i} = P1*P2.Ryy{1,i};
            Pcomp.Ry2{1,i} = P1*P2.Ry2{1,i};
            Pcomp.R2y{1,i} = P1*P2.R2y{1,i};
            
            for j=1:3
                Pcomp.R22{i,j} = P1*P2.R22{i,j};
            end
        end

    else
        if size(P1,2)~=sum(P2.dim(:,1))
            error("Multiplication requires inner dimensions of the operators to match");
        elseif all(all(P2.dim(:,:)==0))
            warning(['Object P in matrix-dopvar2d product M*P has no rows or columns, P.dim(:,:)=0.'...
                        ' The product has ambiguous dimensions; returning the input matrix M.'])
            Pcomp = P1;
            return
        end

        nindx = cumsum(P2.dim(:,1));
        P1_0 = P1(:,1:nindx(1));            P1_x = P1(:,nindx(1)+1:nindx(2));
        P1_y = P1(:,nindx(2)+1:nindx(3));   P1_2 = P1(:,nindx(3)+1:end);
        if isempty(P1_0)
            n1 = size(P1,1) * (P2.dim(1,2)>0);
            P1_0 = zeros(n1,P2.dim(1,1));
        end
        if isempty(P1_x)
            n1 = size(P1,1) * (P2.dim(2,2)>0);
            P1_x = zeros(n1,P2.dim(2,1));
        end
        if isempty(P1_y)
            n1 = size(P1,1) * (P2.dim(3,2)>0);
            P1_y = zeros(n1,P2.dim(3,1));
        end
        if isempty(P1_2)
            n1 = size(P1,1) * (P2.dim(4,2)>0);
            P1_2 = zeros(n1,P2.dim(4,1));
        end
        
        Pcomp.R00 = P1_0*P2.R00;
        Pcomp.R0x = P1_0*P2.R0x;
        Pcomp.R0y = P1_0*P2.R0y;
        Pcomp.R02 = P1_0*P2.R02;
        Pcomp.Rx0 = P1_x*P2.Rx0;
        Pcomp.Rxy = P1_x*P2.Rxy;
        Pcomp.Ry0 = P1_y*P2.Ry0;
        Pcomp.Ryx = P1_y*P2.Ryx;
        Pcomp.R20 = P1_2*P2.R20;
        
        for i=1:3
            Pcomp.Rxx{i,1} = P1_x*P2.Rxx{i,1};
            Pcomp.Rx2{i,1} = P1_x*P2.Rx2{i,1};
            Pcomp.R2x{i,1} = P1_2*P2.R2x{i,1};
            
            Pcomp.Ryy{1,i} = P1_y*P2.Ryy{1,i};
            Pcomp.Ry2{1,i} = P1_y*P2.Ry2{1,i};
            Pcomp.R2y{1,i} = P1_2*P2.R2y{1,i};
            
            for j=1:3
                Pcomp.R22{i,j} = P1_2*P2.R22{i,j};
            end
        end
    end
        
end

end






function [Pc] = PL2L_compose1D(P1,P2,var1,var2,I,indx)                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the composition of two operators on L2[var1(indx)] 
% and returns an equivalent operator on L2[var1(indx)]
% Version 1.0
% Date: 06/21/21
% 
% Inputs:
% P1: Cell describing second opvar object
% P2: Cell describing first opvar object
% var1: Primary variable of the polynomials
% var2: Secondary variable of the polynomials
% I: domain of the variables ds,dt (must be of same size as both ds and dt)
% indx: dimension along which the composition is performed
% 
% Output:
% Returns the composition operator structure.
%
% NOTE: This function is designed only for composition of PI operators on
% just 1 dimension, no other dimensional PI operators

% Determine the limits of the domain
a = I(indx,1);
b = I(indx,2);

% Define a dummy variable for integration
pvar bt1 ;
db = bt1;
ds = var1(indx);
dt = var2(indx);

% Initialize the composition
if indx==1
    Pc = cell(3,1);
else
    arr_dim = [ones(1,indx-1),3];
    Pc = cell(arr_dim);
end

% Describe the compositions of multipliers with integrals
Pc{1} = P1{1}*P2{1};
Pc{2} = P1{1}*P2{2} + P1{2}*subs(P2{1},ds,dt);
Pc{3} = P1{1}*P2{3} + P1{3}*subs(P2{1},ds,dt);
    
    
% Replace intermediate variable by beta
P1{2} = subs(P1{2},dt,db);
P1{3} = subs(P1{3},dt,db);

P2{2} = subs(P2{2},ds,db);
P2{3} = subs(P2{3},ds,db);
    

% Describe compositions of 1D integrals
Pc{2} = Pc{2} + int_simple(P1{2}*P2{3},db,a,dt) ...
              + int_simple(P1{2}*P2{2},db,dt,ds) ...
              + int_simple(P1{3}*P2{2},db,ds,b);
Pc{3} = Pc{3} + int_simple(P1{2}*P2{3},db,a,ds) ...
              + int_simple(P1{3}*P2{3},db,ds,dt) ...
              + int_simple(P1{3}*P2{2},db,dt,b);

end





function [Pc] = PL2L_compose2D(P1,P2,ds,dt,I)                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the composition of two operators on L2[ds(1),ds(2)] 
% and returns an equivalent operator on L2[ds(1),ds(2)]
% Version 1.0
% Date: 02/01/21
% 
% Inputs:
% P1: Cell describing second opvar object
% P2: Cell describing first opvar object
% ds: Primary variables of the polynomials
% dt: Secondary variables of the polynomials
% I: domain of the variables ds,dt (must be of equal size as both ds and dt)
% 
% Output:
% Returns the composition operator structure.
%
% NOTE: This function is designed only for composition of pure 2D PI 
% operators, no other dimensional PI operators

% Determine the limits of the domain
a = I(:,1);
b = I(:,2);

% Define a dummy variable for integration
pvar bt1 bt2;
db = [bt1; bt2];

% Initialize the composition
Pc = cell(3,3);
% Convert potential "double"s to "polynomial"
for j=1:9
    if isa(P1{j},'double')
        P1{j} = polynomial(P1{j});
    end
    if isa(P2{j},'double')
        P2{j} = polynomial(P2{j});
    end
end

% To reduce computation time, we avoid computations with parameters which
% are zero
%ismult1 = any(any(~isequal(P1{1,1},0)));
isint1_x = any(any(P1{2,1}.C~=0)) || any(any(P1{3,1}.C~=0));
isint1_y = any(any(P1{1,2}.C~=0)) || any(any(P1{1,3}.C~=0));
isint1_xy = any(any(P1{2,2}.C~=0)) || any(any(P1{3,2}.C~=0)) ...
             || any(any(P1{2,3}.C~=0)) || any(any(P1{3,3}.C~=0));
         
ismult2 = any(any(P2{1,1}.C~=0));
isint2_x = any(any(P2{2,1}.C~=0)) || any(any(P2{3,1}.C~=0));
isint2_y = any(any(P2{1,2}.C~=0)) || any(any(P2{1,3}.C~=0));
isint2_xy = any(any(P2{2,2}.C~=0)) || any(any(P2{3,2}.C~=0)) ...
             || any(any(P2{2,3}.C~=0)) || any(any(P2{3,3}.C~=0));


% % Describe the compositions of multiplier with PI operator
Pc{1,1} = P1{1,1}*P2{1,1};
Pc{2,1} = P1{1,1}*P2{2,1};
Pc{3,1} = P1{1,1}*P2{3,1};
Pc{1,2} = P1{1,1}*P2{1,2};
Pc{1,3} = P1{1,1}*P2{1,3};
Pc{2,2} = P1{1,1}*P2{2,2};
Pc{3,2} = P1{1,1}*P2{3,2};
Pc{2,3} = P1{1,1}*P2{2,3};
Pc{3,3} = P1{1,1}*P2{3,3};

% % Add compositions of PI operator with multiplier
if ismult2
    if isint1_x
        Pc{2,1} = Pc{2,1} + P1{2,1}*subs(P2{1,1},ds(1),dt(1));
        Pc{3,1} = Pc{3,1} + P1{3,1}*subs(P2{1,1},ds(1),dt(1));
    end
    if isint1_y
        Pc{1,2} = Pc{1,2} + P1{1,2}*subs(P2{1,1},ds(2),dt(2));
        Pc{1,3} = Pc{1,3} + P1{1,3}*subs(P2{1,1},ds(2),dt(2));
    end
    if isint1_xy
        Pc{2,2} = Pc{2,2} + P1{2,2}*subs(P2{1,1},ds,dt);
        Pc{3,2} = Pc{3,2} + P1{3,2}*subs(P2{1,1},ds,dt);
        Pc{2,3} = Pc{2,3} + P1{2,3}*subs(P2{1,1},ds,dt);
        Pc{3,3} = Pc{3,3} + P1{3,3}*subs(P2{1,1},ds,dt);
    end
end

% % Add compositions of 1D y-integrals with 1D x-integrals
if isint2_x
    if isint1_y
        Pc{2,2} = Pc{2,2} + P1{1,2}*subs(P2{2,1},ds(2),dt(2));
        Pc{3,2} = Pc{3,2} + P1{1,2}*subs(P2{3,1},ds(2),dt(2));
        Pc{2,3} = Pc{2,3} + P1{1,3}*subs(P2{2,1},ds(2),dt(2));
        Pc{3,3} = Pc{3,3} + P1{1,3}*subs(P2{3,1},ds(2),dt(2));
    end
    if isint1_x || isint1_xy
        % Replace intermediate variable by beta for further compositions
        P2{2,1} = subs(P2{2,1},ds(1),db(1));
        P2{3,1} = subs(P2{3,1},ds(1),db(1));
    end
end
    
% % Add compositions of 1D x-integrals with 1D y-integrals
if isint2_y
    if isint1_x
        Pc{2,2} = Pc{2,2} + P1{2,1}*subs(P2{1,2},ds(1),dt(1));
        Pc{3,2} = Pc{3,2} + P1{3,1}*subs(P2{1,2},ds(1),dt(1));
        Pc{2,3} = Pc{2,3} + P1{2,1}*subs(P2{1,3},ds(1),dt(1));
        Pc{3,3} = Pc{3,3} + P1{3,1}*subs(P2{1,3},ds(1),dt(1));
    end
    if isint1_y || isint1_xy
        % Replace intermediate variable by beta for further compositions
        P2{1,2} = subs(P2{1,2},ds(2),db(2));
        P2{1,3} = subs(P2{1,3},ds(2),db(2));
    end
end

if isint2_xy && (isint1_x || isint1_y || isint1_xy)
    % Replace intermediate variable by beta for further compositions
    P2{2,2} = subs(P2{2,2},ds,db);
    P2{3,2} = subs(P2{3,2},ds,db);
    P2{2,3} = subs(P2{2,3},ds,db);
    P2{3,3} = subs(P2{3,3},ds,db);
end
    

% % Add composition of 1D x-integral with PI operator
if isint1_x
    % Replace intermediate variable by beta for further compositions
    P1{2,1} = subs(P1{2,1},dt(1),db(1));
    P1{3,1} = subs(P1{3,1},dt(1),db(1));
    
    % Add compositions of 1D x-integrals with 1D x-integrals
    if isint2_x
        Pc{2,1} = Pc{2,1} + int_simple(P1{2,1}*P2{3,1},db(1),a(1),dt(1)) ...
                          + int_simple(P1{2,1}*P2{2,1},db(1),dt(1),ds(1)) ...
                          + int_simple(P1{3,1}*P2{2,1},db(1),ds(1),b(1));
        Pc{3,1} = Pc{3,1} + int_simple(P1{2,1}*P2{3,1},db(1),a(1),ds(1)) ...
                          + int_simple(P1{3,1}*P2{3,1},db(1),ds(1),dt(1)) ...
                          + int_simple(P1{3,1}*P2{2,1},db(1),dt(1),b(1));
    end
    % Add compositions of 1D x-integrals with 2D integrals
    if isint2_xy
        Pc{2,2} = Pc{2,2} + int_simple(P1{2,1}*subs(P2{3,2},db(2),ds(2)),db(1),a(1),dt(1)) ...
                  + int_simple(P1{2,1}*subs(P2{2,2},db(2),ds(2)),db(1),dt(1),ds(1)) ...
                  + int_simple(P1{3,1}*subs(P2{2,2},db(2),ds(2)),db(1),ds(1),b(1));
        Pc{3,2} = Pc{3,2} + int_simple(P1{2,1}*subs(P2{3,2},db(2),ds(2)),db(1),a(1),ds(1)) ...
                  + int_simple(P1{3,1}*subs(P2{3,2},db(2),ds(2)),db(1),ds(1),dt(1)) ...
                  + int_simple(P1{3,1}*subs(P2{2,2},db(2),ds(2)),db(1),dt(1),b(1));
        Pc{2,3} = Pc{2,3} + int_simple(P1{2,1}*subs(P2{3,3},db(2),ds(2)),db(1),a(1),dt(1)) ...
                  + int_simple(P1{2,1}*subs(P2{2,3},db(2),ds(2)),db(1),dt(1),ds(1)) ...
                  + int_simple(P1{3,1}*subs(P2{2,3},db(2),ds(2)),db(1),ds(1),b(1));
        Pc{3,3} = Pc{3,3} + int_simple(P1{2,1}*subs(P2{3,3},db(2),ds(2)),db(1),a(1),ds(1)) ...
                  + int_simple(P1{3,1}*subs(P2{3,3},db(2),ds(2)),db(1),ds(1),dt(1)) ...
                  + int_simple(P1{3,1}*subs(P2{2,3},db(2),ds(2)),db(1),dt(1),b(1));
    end
end

% % Add composition of 1D y-integral with PI operator
if isint1_y
    % Replace intermediate variable by beta for further compositions
    P1{1,2} = subs(P1{1,2},dt(2),db(2));
    P1{1,3} = subs(P1{1,3},dt(2),db(2));
    
    % Add compositions of 1D y-integrals with 1D y-integrals
    if isint2_y
        Pc{1,2} = Pc{1,2} + int_simple(P1{1,2}*P2{1,3},db(2),a(2),dt(2)) ...
                          + int_simple(P1{1,2}*P2{1,2},db(2),dt(2),ds(2)) ...
                          + int_simple(P1{1,3}*P2{1,2},db(2),ds(2),b(2));
        Pc{1,3} = Pc{1,3} + int_simple(P1{1,2}*P2{1,3},db(2),a(2),ds(2)) ...
                          + int_simple(P1{1,3}*P2{1,3},db(2),ds(2),dt(2)) ...
                          + int_simple(P1{1,3}*P2{1,2},db(2),dt(2),b(2));
    end
    % Add compositions of 1D y-integrals with 2D integrals
    if isint2_xy
        Pc{2,2} = Pc{2,2} + int_simple(P1{1,2}*subs(P2{2,3},db(1),ds(1)),db(2),a(2),dt(2)) ...
                  + int_simple(P1{1,2}*subs(P2{2,2},db(1),ds(1)),db(2),dt(2),ds(2)) ...
                  + int_simple(P1{1,3}*subs(P2{2,2},db(1),ds(1)),db(2),ds(2),b(2));
        Pc{3,2} = Pc{3,2} + int_simple(P1{1,2}*subs(P2{3,3},db(1),ds(1)),db(2),a(2),dt(2)) ...
                  + int_simple(P1{1,2}*subs(P2{3,2},db(1),ds(1)),db(2),dt(2),ds(2)) ...
                  + int_simple(P1{1,3}*subs(P2{3,2},db(1),ds(1)),db(2),ds(2),b(2));
        Pc{2,3} = Pc{2,3} + int_simple(P1{1,2}*subs(P2{2,3},db(1),ds(1)),db(2),a(2),ds(2)) ...
                  + int_simple(P1{1,3}*subs(P2{2,3},db(1),ds(1)),db(2),ds(2),dt(2)) ...
                  + int_simple(P1{1,3}*subs(P2{2,2},db(1),ds(1)),db(2),dt(2),b(2));
        Pc{3,3} = Pc{3,3} + int_simple(P1{1,2}*subs(P2{3,3},db(1),ds(1)),db(2),a(2),ds(2)) ...
                  + int_simple(P1{1,3}*subs(P2{3,3},db(1),ds(1)),db(2),ds(2),dt(2)) ...
                  + int_simple(P1{1,3}*subs(P2{3,2},db(1),ds(1)),db(2),dt(2),b(2));
    end
end

% % Add composition of 2D integral with PI operator
if isint1_xy
    % Replace intermediate variable by beta for further compositions
    P1{2,2} = subs(P1{2,2},dt,db);
    P1{3,2} = subs(P1{3,2},dt,db);
    P1{2,3} = subs(P1{2,3},dt,db);
    P1{3,3} = subs(P1{3,3},dt,db);
    
    % Add compositions of 2D integrals with 1D x-integrals
    if isint2_x
        Pc{2,2} = Pc{2,2} + int_simple(subs(P1{2,2},db(2),dt(2))*subs(P2{3,1},ds(2),dt(2)),db(1),a(1),dt(1)) ...
                  + int_simple(subs(P1{2,2},db(2),dt(2))*subs(P2{2,1},ds(2),dt(2)),db(1),dt(1),ds(1)) ...
                  + int_simple(subs(P1{3,2},db(2),dt(2))*subs(P2{2,1},ds(2),dt(2)),db(1),ds(1),b(1));
        Pc{3,2} = Pc{3,2} + int_simple(subs(P1{2,2},db(2),dt(2))*subs(P2{3,1},ds(2),dt(2)),db(1),a(1),ds(1)) ...
                  + int_simple(subs(P1{3,2},db(2),dt(2))*subs(P2{3,1},ds(2),dt(2)),db(1),ds(1),dt(1)) ...
                  + int_simple(subs(P1{3,2},db(2),dt(2))*subs(P2{2,1},ds(2),dt(2)),db(1),dt(1),b(1));
        Pc{2,3} = Pc{2,3} + int_simple(subs(P1{2,3},db(2),dt(2))*subs(P2{3,1},ds(2),dt(2)),db(1),a(1),dt(1)) ...
                  + int_simple(subs(P1{2,3},db(2),dt(2))*subs(P2{2,1},ds(2),dt(2)),db(1),dt(1),ds(1)) ...
                  + int_simple(subs(P1{3,3},db(2),dt(2))*subs(P2{2,1},ds(2),dt(2)),db(1),ds(1),b(1)); 
        Pc{3,3} = Pc{3,3} + int_simple(subs(P1{2,3},db(2),dt(2))*subs(P2{3,1},ds(2),dt(2)),db(1),a(1),ds(1)) ...
                  + int_simple(subs(P1{3,3},db(2),dt(2))*subs(P2{3,1},ds(2),dt(2)),db(1),ds(1),dt(1)) ...
                  + int_simple(subs(P1{3,3},db(2),dt(2))*subs(P2{2,1},ds(2),dt(2)),db(1),dt(1),b(1));
    end
    % Decribe compositions of 2D integrals with 1D y-integrals
    if isint2_y
        Pc{2,2} = Pc{2,2} + int_simple(subs(P1{2,2},db(1),dt(1))*subs(P2{1,3},ds(1),dt(1)),db(2),a(2),dt(2)) ...
                  + int_simple(subs(P1{2,2},db(1),dt(1))*subs(P2{1,2},ds(1),dt(1)),db(2),dt(2),ds(2)) ...
                  + int_simple(subs(P1{2,3},db(1),dt(1))*subs(P2{1,2},ds(1),dt(1)),db(2),ds(2),b(2));
        Pc{3,2} = Pc{3,2} + int_simple(subs(P1{3,2},db(1),dt(1))*subs(P2{1,3},ds(1),dt(1)),db(2),a(2),dt(2)) ...
                  + int_simple(subs(P1{3,2},db(1),dt(1))*subs(P2{1,2},ds(1),dt(1)),db(2),dt(2),ds(2)) ...
                  + int_simple(subs(P1{3,3},db(1),dt(1))*subs(P2{1,2},ds(1),dt(1)),db(2),ds(2),b(2));
        Pc{2,3} = Pc{2,3} + int_simple(subs(P1{2,2},db(1),dt(1))*subs(P2{1,3},ds(1),dt(1)),db(2),a(2),ds(2)) ...
                  + int_simple(subs(P1{2,3},db(1),dt(1))*subs(P2{1,3},ds(1),dt(1)),db(2),ds(2),dt(2)) ...
                  + int_simple(subs(P1{2,3},db(1),dt(1))*subs(P2{1,2},ds(1),dt(1)),db(2),dt(2),b(2));
        Pc{3,3} = Pc{3,3} + int_simple(subs(P1{3,2},db(1),dt(1))*subs(P2{1,3},ds(1),dt(1)),db(2),a(2),ds(2)) ...
                  + int_simple(subs(P1{3,3},db(1),dt(1))*subs(P2{1,3},ds(1),dt(1)),db(2),ds(2),dt(2)) ...
                  + int_simple(subs(P1{3,3},db(1),dt(1))*subs(P2{1,2},ds(1),dt(1)),db(2),dt(2),b(2));   
    end
    % Describe compositions of 2D integrals with 2D integrals
    if isint2_xy
        Pc{2,2} = Pc{2,2} + int_simple(int_simple(P1{2,2}*P2{3,3},db(1),a(1),dt(1)) ...
                        + int_simple(P1{2,2}*P2{2,3},db(1),dt(1),ds(1)) ...
                        + int_simple(P1{3,2}*P2{2,3},db(1),ds(1),b(1)),db(2),a(2),dt(2)) ...
                  + int_simple(int_simple(P1{2,2}*P2{3,2},db(1),a(1),dt(1)) ...
                        + int_simple(P1{2,2}*P2{2,2},db(1),dt(1),ds(1)) ...
                        + int_simple(P1{3,2}*P2{2,2},db(1),ds(1),b(1)),db(2),dt(2),ds(2)) ...
                  + int_simple(int_simple(P1{2,3}*P2{3,2},db(1),a(1),dt(1)) ...
                        + int_simple(P1{2,3}*P2{2,2},db(1),dt(1),ds(1)) ...
                        + int_simple(P1{3,3}*P2{2,2},db(1),ds(1),b(1)),db(2),ds(2),b(2));
        Pc{3,2} = Pc{3,2} + int_simple(int_simple(P1{2,2}*P2{3,3},db(1),a(1),ds(1)) ...
                        + int_simple(P1{3,2}*P2{3,3},db(1),ds(1),dt(1)) ...
                        + int_simple(P1{3,2}*P2{2,3},db(1),dt(1),b(1)),db(2),a(2),dt(2)) ...
                  + int_simple(int_simple(P1{2,2}*P2{3,2},db(1),a(1),ds(1)) ...
                        + int_simple(P1{3,2}*P2{3,2},db(1),ds(1),dt(1)) ...
                        + int_simple(P1{3,2}*P2{2,2},db(1),dt(1),b(1)),db(2),dt(2),ds(2)) ...
                  + int_simple(int_simple(P1{2,3}*P2{3,2},db(1),a(1),ds(1)) ...
                        + int_simple(P1{3,3}*P2{3,2},db(1),ds(1),dt(1)) ...
                        + int_simple(P1{3,3}*P2{2,2},db(1),dt(1),b(1)),db(2),ds(2),b(2));
        Pc{2,3} = Pc{2,3} + int_simple(int_simple(P1{2,2}*P2{3,3},db(1),a(1),dt(1)) ...
                        + int_simple(P1{2,2}*P2{2,3},db(1),dt(1),ds(1)) ...
                        + int_simple(P1{3,2}*P2{2,3},db(1),ds(1),b(1)),db(2),a(2),ds(2)) ...
                  + int_simple(int_simple(P1{2,3}*P2{3,3},db(1),a(1),dt(1)) ...
                        + int_simple(P1{2,3}*P2{2,3},db(1),dt(1),ds(1)) ...
                        + int_simple(P1{3,3}*P2{2,3},db(1),ds(1),b(1)),db(2),ds(2),dt(2)) ...
                  + int_simple(int_simple(P1{2,3}*P2{3,2},db(1),a(1),dt(1)) ...
                        + int_simple(P1{2,3}*P2{2,2},db(1),dt(1),ds(1)) ...
                        + int_simple(P1{3,3}*P2{2,2},db(1),ds(1),b(1)),db(2),dt(2),b(2));
        Pc{3,3} = Pc{3,3} + int_simple(int_simple(P1{2,2}*P2{3,3},db(1),a(1),ds(1)) ...
                        + int_simple(P1{3,2}*P2{3,3},db(1),ds(1),dt(1)) ...
                        + int_simple(P1{3,2}*P2{2,3},db(1),dt(1),b(1)),db(2),a(2),ds(2)) ...
                  + int_simple(int_simple(P1{2,3}*P2{3,3},db(1),a(1),ds(1)) ...
                        + int_simple(P1{3,3}*P2{3,3},db(1),ds(1),dt(1)) ...
                        + int_simple(P1{3,3}*P2{2,3},db(1),dt(1),b(1)),db(2),ds(2),dt(2)) ...
                  + int_simple(int_simple(P1{2,3}*P2{3,2},db(1),a(1),ds(1)) ...
                        + int_simple(P1{3,3}*P2{3,2},db(1),ds(1),dt(1)) ...
                        + int_simple(P1{3,3}*P2{2,2},db(1),dt(1),b(1)),db(2),dt(2),b(2));
    end
end                    

end





function [PoM] = PIoMat_composeND(P1,P2,ds,dt,I)                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function applies a PI operator to a matrix-valued function
% to return a matrix-valued function
% Version 1.0
% Date: 02/01/21
% 
% Inputs:
% P1: Polynomial or cell describing opvar
% P2: Polynomial or cell describing opvar
% ds: Primary variables of the polynomials
% dt: Secondary variables of the polynomials
% I: domain of the variables ds,dt (must be of equal size as both ds and dt)
%
% If P1 is opvar and P2 polynomial, PoM describes P1 applied to P2
% If P1 is polynomial and P2 opvar, PoM describes P2^T applied (from the
% right) to P1.
% NOTE: both operators can be opvars, if they operate along different
% dimensions. That is, in each dimension, one operator must act as opvar,
% the other as polynomial.
% 
% Output:
% Returns the composition operator structure.

    
if isa(P1,'cell') && ~isa(P2,'cell')
    vtot = ndims(P1);
    sz = size(P1);
    % Empty variable
    MT = pvar('MT');
    
    PoM = zeros(size(P1{1,1},1),size(P2,2));
    
    for i = 1:prod(sz)
    
        crd = cell(1,vtot);
        [crd{:}] = ind2sub(sz,i);
        a_indx = find([crd{:}]==2);
        b_indx = find([crd{:}]==3);
    
        P1i = P1{i}*subs(P2,[ds([a_indx,b_indx]);MT],[dt([a_indx,b_indx]);MT]);
    
        for d=a_indx
            P1i =  int_simple(P1i,dt(d),I(d,1),ds(d));
        end
        for d=b_indx
            P1i =  int_simple(P1i,dt(d),ds(d),I(d,2));
        end
        
        PoM = PoM + P1i;
    
    end
    
elseif ~isa(P1,'cell') && isa(P2,'cell')
    vtot = ndims(P2);
    sz = size(P2);
    
    PoM = zeros(size(P1,1),size(P2{1,1},2));
    
    for i = 1:prod(sz)
    
        crd = cell(1,vtot);
        [crd{:}] = ind2sub(sz,i);
        a_indx = find([crd{:}]==3);
        b_indx = find([crd{:}]==2);
    
        P1i = P1;
        P2i = P2{i};
        for d=[a_indx,b_indx]
            P1i = subs(P1i,ds(d),dt(d));
            P2i = var_swap(P2i,ds(d),dt(d));
        end
        P2i = P1i*P2i;
    
        for d=a_indx
            P2i =  int_simple(P2i,dt(d),I(d,1),ds(d));
        end
        for d=b_indx
            P2i =  int_simple(P2i,dt(d),ds(d),I(d,2));
        end
        
        PoM = PoM + P2i;
    
    end
    
elseif isa(P1,'cell') && isa(P2,'cell')
    vtot = max(ndims(P1),ndims(P2));
    sz_1 = size(P1);
    sz_2 = size(P2);
    % Empty variable
    MT = pvar('MT');
    
    PoM = zeros(size(P1{1,1},1),size(P2{1,1},2));
    
    for i1 = 1:prod(sz_1)
        crd1 = cell(1,vtot);
        [crd1{:}] = ind2sub(sz_1,i1);
        a_indx1 = find([crd1{:}]==2);
        b_indx1 = find([crd1{:}]==3);
        
        for i2 = 1:prod(sz_2)
    
            crd2 = cell(1,vtot);
            [crd2{:}] = ind2sub(sz_2,i2);
            a_indx2 = find([crd2{:}]==3);
            b_indx2 = find([crd2{:}]==2);
        
            P1i = P1{i1};
            P2i = P2{i2};
            for d=[a_indx2,b_indx2]
                P1i = subs(P1i,ds(d),dt(d));
                P2i = var_swap(P2i,ds(d),dt(d));
            end
    
            Pni = P1i*subs(P2i,[ds([a_indx1,b_indx1]);MT],[dt([a_indx1,b_indx1]);MT]);
    
            for d=[a_indx1,a_indx2]
                Pni =  int_simple(Pni,dt(d),I(d,1),ds(d));
            end
            for d=[b_indx1,b_indx2]
                Pni =  int_simple(Pni,dt(d),ds(d),I(d,2));
            end
            PoM = PoM + Pni;
        end
    end
    
end

end
