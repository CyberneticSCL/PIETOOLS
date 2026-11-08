function [prog,Zop] = lpivar_2d(prog,dim,d,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Zop] = lpivar_2d(prog,dim,d,options) declares an indefinite
% 0112-PI operator,
%
% Zop = [R00 R0x R0y R02]
%       [Rx0 Rxx Rxy Rx2]
%       [Ry0 Ryx Ryy Ry2]
%       [R20 R2x R2y R22],
%
% where R00 is a constant matrix;
%       R0x, Rx0, and Rxx{1} are defined by monomials Zxo(ss1);
%       Rxx{1}, Rxx{2} are defined by monomials Zxa(ss1,tt1), Zxb(ss1,tt1);
%       R0y, Ry0, and Ryy{1} are defined by monomials Zyo(ss2);
%       Ryy{1}, Ryy{2} are defined by monomials Zya(ss2,tt2), Zyb(ss2,tt2);
%       Rxy, Ryx, Rx2{1}, R2x{1}, Ry2{1}, R2y{1}, R22{1,1} are defined by
%           monomials Z2oo(ss1,ss2);
%       Rx2{2}, Rx2{3}, R2x{2}, R2x{3}, R22{2,1}, R22{3,1} are defined by
%           monomials Z2ao(ss1,ss2,tt1) and Z2bo(ss1,ss2,tt1);
%       Ry2{2}, Ry2{3}, R2y{2}, R2y{3}, R22{1,2}, R22{1,3} are defined by
%           monomials Z2oa(ss1,ss2,tt2) and Z2ob(ss1,ss2,tt2);
%       R22{2,2}, R22{3,2}, R22{2,3}, R22{3,3} are defined by monomials
%           Z2aa(ss1,ss2,tt1,tt2), ..., Z2bb(ss1,ss2,tt1,tt2).
%
% 
% INPUT 
%   prog: LPI program to modify.
%   dim:  4x2 array providing the dimensions of the PI operators
%
%   d: structure with fields dx, dy, d2, providing the maximal degrees of
%      the different variables in each of the monomials. In particular:
%
%   dx: 3x1-cell structure, 
%   dx{1}: 1x1 scalar providing maximal degree of x in Zxo(x)
%   dx{2}: 3x1 array providing maximal degree of [x;tt;x*tt] in Zxa(x,tt)
%   dx{3}: 3x1 array providing maximal degree of [x;tt;x*tt] in Zxb(x,tt)
%       Here the maximal degree in x*tt corresponds to the maximal joint
%       degree in variables x and tt in any monomial in Zxb
%
%   dy: 1x3-cell structure, 
%   dy{1}: 1x1 scalar providing maximal degree of y in Zyo(y)
%   dy{2}: 1x3 array providing maximal degree of [y;nu;y*nu] in Zya(y,nu)
%   dy{3}: 1x3 array providing maximal degree of [y;nu;y*nu] in Zyb(y,nu)
%
%   d2: 3x3 cell describing the maximal degrees for the parameters Z2 as
%
%       [ d2{1,1}, d2{1,2}, d2{1,3} ] ~ [ Z2oo, Z2oa, Z2ob ]
%       [ d2{2,1}, d2{2,2}, d2{2,3} ] ~ [ Z2ao, Z2aa, Z2ab ]
%       [ d2{3,1}, d2{3,2}, d2{3,3} ] ~ [ Z2bo, Z2ba, Z2bb ]
%
%   here, d2{1,1} is a 2x2 array, d2{i,1} are 4x2 arrays, d2{1,j} are 2x4
%   arrays, and d2{i,j} are 4x4 arrays for i,j in {2,3}. These arrays
%   provide maximal individual and joint degrees as
%
%                 0 |           y |          nu |         y*nu
%          ---------|-------------|-------------|-----------------
%                 x |         x*y |        x*nu |       x*y*nu
%          ---------|-------------|-------------|-----------------
%                tt |        tt*y |       tt*nu |      tt*y*nu 
%          ---------|-------------|-------------|-----------------
%              x*tt |      x*tt*y |     x*tt*nu |    x*tt*y*nu 
%
%   Note that the first element of this array should always be zero. Rows
%   and columns associated to variables that do not appear in the monomial
%   can be omitted, or should otherwise be set to zero.
%
%   options.ismultiplier: Binary value set to 1 if the PI operator should
%                       not include any integrator terms
%
%   options.sep is a length 6 binary vector where
%      options.sep(1) = 1 if Rxx{2} = Rxx{3}
%      options.sep(2) = 1 if Ryy{2} = Ryy{3}
%      options.sep(3) = 1 if R22{2,1} = R22{3,1}
%      options.sep(4) = 1 if R22{1,2} = R22{1,3}
%      options.sep(5) = 1 if R22{2,2} = R22{3,2} and R22{2,3} = R22{3,3}
%      options.sep(6) = 1 if R22{2,2} = R22{2,3} and R22{3,2} = R22{3,3}
%
%   options.exclude is a 4x4 array of logical values specifying whether to
%       exclude and of the maps R00 through R22 in the operator Zop.
%       Can also be specified as 4x4 cell, where now elements {2,2}, {3,3},
%       {2,4}, {4,2}, {3,4} and {4,3} are 1x3 cells and element {4,4} is a
%       3x3 cell, specifying whether to exclude particular components of
%       the maps Rxx, Ryy, Rx2, R2x, Ry2, R2y, and R22, respectively.
% 
% OUTPUT 
%   prog: modified LPI program with new variables and constraints
%   Zop:  dopvar2d object describing an indefinite PI operator
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - lpivar_2d
%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 07/01/2024
%                   (ss1,ss2,tt1,tt2) --> (s1,s2,s1_dum,s2_dum);
% DJ, 04/09/2026: Add 'options.exclude' to exclude certain maps;


% % First, check the spatial domain on which the program is defined.
if ~isfield(prog,'dom') || size(prog.dom,1)==0
    error('The program structure does not include a spatial domain -- please use ''lpiprogram'' to initialize your program');
else
    I = prog.dom;
end
if any(size(I)~=[2,2])
    error('For declaring a 2D PI operator variable, the spatial domain in the LPI program structure must be of of size 2x2.')
end

% Initialize default degrees
dx = {1;[1;1;1];[1;1;1]};
dy = {1,[1,1,1],[1,1,1]};
d2 = {[0,1;1,2],          [0,1,1,1;1,2,2,2], [0,1,1,1;1,2,2,2];
      [0,1,1,1;1,2,2,2]', [0,1,1,1;1,2,2,2;1,2,2,2;1,2,2,2], [0,1,1,1;1,2,2,2;1,2,2,2;1,2,2,2];
      [0,1,1,1;1,2,2,2]', [0,1,1,1;1,2,2,2;1,2,2,2;1,2,2,2], [0,1,1,1;1,2,2,2;1,2,2,2;1,2,2,2]};

% Initialize default exclusiong options
use_sep  = zeros(1,6);
use_sym = false;
excludeL = {0,0,0,0; 
            0,[0,0,0],0,[0,0,0];
            0,0,[0,0,0],[0,0,0];
            0,[0,0,0],[0,0,0],[0,0,0;0,0,0;0,0,0]};

% % Extract the input arguments
switch nargin
    case 1
        error('Not enough inputs!')
    case 2
        fprintf('\n Warning: No degrees are specified. Continuing with default values. \n')
    case 3
        
    case 4
        if isfield(options,'psatz')
            warning('There is no psatz option for indefinite PI operators.')
        end
        if isfield(options,'sep')
            if isscalar(options.sep)
                use_sep = options.sep*ones(1,6);
            else
                use_sep = options.sep;
            end
        end
        if isfield(options,'sym')
            use_sym = options.sym;
        end
        if isfield(options,'exclude')
            excludeL = options.exclude;
        end
end


% % % Check if inputs are properly specified
% % Start with the dimensions of the operator, should be a 4x1 or 4x2 array
% P\in R^(n0 x n0), Rxx\in L_2[x]^(nx x nx), Ryy\in L_2[y]^(ny x ny), 
%   R22\in L_2[x,y]^(n2 x n2),
if isscalar(dim)
    warning('Only 1 dimension is provided, assuming this to refer to input and output PDE state')
    m = [0;0;0;dim(1)]';
    n = [0;0;0;dim(1)]';
elseif all(size(dim)==[1,2])
    warning('Only 1 row and column dimension is provided, assuming this to refer to the PDE state')
    m = [0;0;0;dim(1)];
    n = [0;0;0;dim(2)];
elseif all(size(dim)==[2,2])
    warning('Only 2 row and column dimensions are provided, assuming these to refer to the ODE and PDE states')
    m = [dim(1,1);0;0;dim(2,2)];
    n = [dim(1,2);0;0;dim(2,2)];
elseif all(size(dim)==[4,2])
    m = dim(:,1);
    n = dim(:,2);
else
    error("Desired dimensions of the operator must be specified as 4x2 array.")
end
m0 = m(1);    mx = m(2);    my = m(3);    m2 = m(4);
n0 = n(1);    nx = n(2);    ny = n(3);    n2 = n(4);
if n0==0 && nx==0 && ny==0 && n2==0
    error('All input dimensions are zero')
end
if m0==0 && mx==0 && my==0 && m2==0
    error('All output dimensions are zero')
end

% Check that degrees are properly specified
if nargin>=3
    [dx,dy,d2] = check_degs(d);
end

% Check that the options are appropriately specified
% % Check if any parameters are to be excluded (set to 0)
if isnumeric(excludeL)
    if numel(excludeL)==16
        % Assume parameters to exclude are specified as per poslpivar2d
        % format
        exclude_diag = [excludeL(1);
                        all(excludeL(2:4));
                        all(excludeL(5:7));
                        all(excludeL(8:end))];
        excludeL_new = num2cell(exclude_diag | exclude_diag');
        excludeL_new{2,2} = reshape(excludeL(2:4),1,[]);
        excludeL_new{3,3} = reshape(excludeL(5:7),1,[]);
        excludeL_new{4,4} = [excludeL(8),excludeL(11),excludeL(12);
                         excludeL(9),excludeL(13),excludeL(15);
                         excludeL(10),excludeL(14),excludeL(16)];
        excludeL_new{2,4} = excludeL_new{2,2} | excludeL_new{4,4}(:,1)';
        excludeL_new{4,2} = excludeL_new{2,4};
        excludeL_new{3,4} = excludeL_new{3,2} | excludeL_new{4,4}(1,:);
        excludeL_new{4,3} = excludeL_new{3,4};
        excludeL = excludeL_new;
    elseif any(size(excludeL)~=[4,4])
        error("List of parameters to exclude should be specified as 4x4 binary array or 4x4 cell.")
    else
    excludeL = {excludeL(1,1),excludeL(1,2),excludeL(1,3),excludeL(1,4);
        excludeL(2,1),excludeL(2,2)*ones(1,3),excludeL(2,3),excludeL(2,4)*ones(1,3);
        excludeL(3,1),excludeL(3,2),excludeL(3,3)*ones(1,3),excludeL(3,4)*ones(1,3);
        excludeL(4,1),excludeL(4,2)*ones(1,3),excludeL(4,3)*ones(1,3),...
        excludeL(4,4)*ones(3,3)};
    end
elseif isa(excludeL,'cell')
    if any(size(excludeL)~=[4,4])
        error("List of parameters to exclude should be specified as 4x4 cell.")
    end
end
% Exclude empty parameters
if m(1)==0
    excludeL(1,:) = repmat({true},1,4);
end
if n(1)==0
    excludeL(:,1) = repmat({true},4,1);
end
if m(2)==0
    excludeL(2,:) = {true,true(1,3),true,true(1,3)};
end
if n(2)==0
    excludeL(:,2) = {true;true(1,3);true;true(1,3)};
end
if m(3)==0
    excludeL(3,:) = {true,true,true(1,3),true(1,3)};
end
if n(3)==0
    excludeL(:,3) = {true;true;true(1,3);true(1,3)};
end
if m(4)==0
    excludeL(4,:) = {true,true(1,3),true(1,3),true(3,3)};
end
if n(4)==0
    excludeL(:,4) = {true;true(1,3);true(1,3);true(3,3)};
end

% % Check whether the operator should be symmetric
if use_sym
    if any(n(:,1)~=n(:,2))
        error('Row and column dimensions must match for symmetric operators.')
    end
    matrixstr = 'symmetric';
    excludeL{[2;3;4;7]} = true;
    excludeL{[8;12]} = true(1,3);
    excludeL{2,2}(3) = true;
    excludeL{3,3}(3) = true;
    excludeL{4,4}([3;6;7;9]) = true;
else
    matrixstr = '';
end

% Define the primary variables
if isempty(prog.vartable) || length(prog.vartable)<4
    error("LPI program has not been properly initialized for 2D dopvar declaration; initialize the program with 'lpiprogram', declaring at least two spatial variables.")
else
    var = polynomial(prog.vartable(1:4,1));
    var1 = var((1:2)');     var2 = var((3:4)');
end
vars = [var1,var2];



% % % Build the indefinite operator

% Construct monomials given the degrees
[Zx,Zy,Z2] = build_monomials_2D(vars,dx,dy,d2);


% % Initialize the indefinite operator variable.
dopvar2d Zop;
Zop.dim = [m,n];
Zop.I = I; 
Zop.var1 = var1;    Zop.var2 = var2;

% % Define the parameters and add the decision variables to the program
% % structure

% Parameters mapping to R^m0
if ~excludeL{1,1}
    [prog,R00] = sospolymatrixvar(prog,monomials(var1,0),[m0,n0],matrixstr);
    Zop.R00 = R00;
end
if ~excludeL{1,2}
    [prog,R0x] = sospolymatrixvar(prog,Zx{1},[m0,nx]);
    Zop.R0x = R0x;
end
if ~excludeL{1,3}
    [prog,R0y] = sospolymatrixvar(prog,Zy{1},[m0,ny]);
    Zop.R0y = R0y;
end
if ~excludeL{1,4}
    [prog,R02] = sospolymatrixvar(prog,Z2{1,1},[m0,n2]);
    Zop.R02 = R02;
end

% Parameters mapping to L_2^mx[x]
if use_sym
    Zop.Rx0 = Zop.R0x';
elseif ~excludeL{2,1}
    [prog,Rx0] = sospolymatrixvar(prog,Zx{1},[mx,n0]);
    Zop.Rx0 = Rx0;
end
if ~excludeL{2,2}(1)
    [prog,Rxx_0] = sospolymatrixvar(prog,Zx{1},[mx,nx],matrixstr);
    Zop.Rxx{1} = Rxx_0;
end
if ~excludeL{2,2}(2)
    [prog,Rxx_1] = sospolymatrixvar(prog,Zx{2},[mx,nx]);
    Zop.Rxx{2} = Rxx_1;
end
if use_sep(1)
    Zop.Rxx{3} = Zop.Rxx{2};
elseif use_sym
    Zop.Rxx{3} = var_swap(Zop.Rxx{2}',var1(1),var2(1));
elseif ~excludeL{2,2}(3)
    [prog,Rxx_2] = sospolymatrixvar(prog,Zx{3},[mx,nx]);
    Zop.Rxx{3} = Rxx_2;
end
if ~excludeL{2,3}
    [prog,Rxy] = sospolymatrixvar(prog,Z2{1,1},[mx,ny]);
    Zop.Rxy = Rxy;
end
if ~excludeL{2,4}(1)
    [prog,Rx2_0] = sospolymatrixvar(prog,Z2{1,1},[mx,n2]);
    Zop.Rx2{1} = Rx2_0;
end
if ~excludeL{2,4}(2)
    [prog,Rx2_1] = sospolymatrixvar(prog,Z2{2,1},[mx,n2]);
    Zop.Rx2{2} = Rx2_1;
end
if use_sep(3)
    Zop.Rx2{3} = Zop.Rx2{2};
elseif ~excludeL{2,4}(3)
    [prog,Rx2_2] = sospolymatrixvar(prog,Z2{3,1},[mx,n2]);
    Zop.Rx2{3} = Rx2_2;
end

% Parameters mapping to L_2^my[y]
if use_sym
    Zop.Ry0 = Zop.R0y';
elseif ~excludeL{3,1}
    [prog,Ry0] = sospolymatrixvar(prog,Zy{1},[my,n0]);
    Zop.Ry0 = Ry0;
end
if use_sym
    Zop.Ryx = Zop.Rxy';
elseif ~excludeL{3,2}
    [prog,Ryx] = sospolymatrixvar(prog,Z2{1,1},[my,nx]);
    Zop.Ryx = Ryx;
end
if ~excludeL{3,3}(1)
    [prog,Ryy_0] = sospolymatrixvar(prog,Zy{1},[my,ny]);
    Zop.Ryy{1} = Ryy_0;
end
if ~excludeL{3,3}(2)
    [prog,Ryy_1] = sospolymatrixvar(prog,Zy{2},[my,ny]);
    Zop.Ryy{2} = Ryy_1;
end
if use_sep(2)
    Zop.Ryy{3} = Zop.Ryy{2};
elseif use_sym
    Zop.Ryy{3} = var_swap(Zop.Ryy{2}',var1(2),var2(2));
elseif ~excludeL{3,3}(3)
    [prog,Ryy_2] = sospolymatrixvar(prog,Zy{3},[my,ny]);
    Zop.Ryy{3} = Ryy_2;
end
if ~excludeL{3,4}(1)
    [prog,Ry2_0] = sospolymatrixvar(prog,Z2{1,1},[my,n2]);
    Zop.Ry2{1} = Ry2_0;
end
if ~excludeL{3,4}(2)
    [prog,Ry2_1] = sospolymatrixvar(prog,Z2{1,2},[my,n2]);
    Zop.Ry2{2} = Ry2_1;
end
if use_sep(4)
    Zop.Ry2{3} = Zop.Ry2{2};
elseif ~excludeL{3,4}(3)
    [prog,Ry2_2] = sospolymatrixvar(prog,Z2{1,3},[my,n2]);
    Zop.Ry2{3} = Ry2_2;
end

% Parameters mapping to L_2^m2[x,y]
if use_sym
    Zop.R20 = Zop.R02';
    Zop.R2x{1} = Zop.Rx2{1}';
    Zop.R2x{2} = var_swap(Zop.Rx2{3}',var1(1),var2(1));
    Zop.R2x{3} = var_swap(Zop.Rx2{2}',var1(1),var2(1));
    Zop.R2y{1} = Zop.Ry2{1}';
    Zop.R2y{2} = var_swap(Zop.Ry2{3}',var1(1),var2(1));
    Zop.R2y{3} = var_swap(Zop.Ry2{2}',var1(1),var2(1));
else
    if ~excludeL{4,1}
        [prog,R20] = sospolymatrixvar(prog,Z2{1,1},[m2,n0]);
        Zop.R20 = R20;
    end
    if ~excludeL{4,2}(1)
        [prog,R2x_0] = sospolymatrixvar(prog,Z2{1,1},[m2,nx]);
        Zop.R2x{1} = R2x_0;
    end
    if ~excludeL{4,2}(2)
        [prog,R2x_1] = sospolymatrixvar(prog,Z2{2,1},[m2,nx]);
        Zop.R2x{2} = R2x_1;
    end
    if use_sep(3)
        Zop.R2x{3} = Zop.R2x{2};
    elseif ~excludeL{4,2}(3)
        [prog,R2x_2] = sospolymatrixvar(prog,Z2{3,1},[m2,nx]);
        Zop.R2x{3} = R2x_2;
    end
    if ~excludeL{4,3}(1)
        [prog,R2y_0] = sospolymatrixvar(prog,Z2{1,1},[m2,ny]);
        Zop.R2y{1} = R2y_0;
    end
    if ~excludeL{4,3}(2)
        [prog,R2y_1] = sospolymatrixvar(prog,Z2{1,2},[m2,ny]);
        Zop.R2y{2} = R2y_1;
    end
    if use_sep(4)
        Zop.R2y{3} = Zop.R2y{2};
    elseif ~excludeL{4,3}(3)
        [prog,R2y_2] = sospolymatrixvar(prog,Z2{1,3},[m2,ny]);
        Zop.R2y{3} = R2y_2;
    end
end
if ~excludeL{4,4}(1,1)
    [prog,R22_00] = sospolymatrixvar(prog,Z2{1,1},[m2,n2],matrixstr);
    Zop.R22{1,1} = R22_00;
end
if ~excludeL{4,4}(2,1)
    [prog,R22_10] = sospolymatrixvar(prog,Z2{2,1},[m2,n2]);
    Zop.R22{2,1} = R22_10;
end
if ~excludeL{4,4}(1,2)
    [prog,R22_01] = sospolymatrixvar(prog,Z2{1,2},[m2,n2]);
    Zop.R22{1,2} = R22_01;
end
if ~excludeL{4,4}(2,2)
    [prog,R22_11] = sospolymatrixvar(prog,Z2{2,2},[m2,n2]);
    Zop.R22{2,2} = R22_11;
end
if use_sep(6)
    Zop.R22{2,3} = Zop.R22{2,2};
elseif ~excludeL{4,4}(2,3)
    [prog,R22_12] = sospolymatrixvar(prog,Z2{2,3},[m2,n2]);
    Zop.R22{2,3} = R22_12;
end
if use_sym
    Zop.R22{3,1} = var_swap(Zop.R22{2,1}',var1(1),var2(1));
    Zop.R22{1,3} = var_swap(Zop.R22{1,2}',var1(2),var2(2));
    Zop.R22{3,3} = var_swap(var_swap(Zop.R22{2,2}',var1(1),var2(1)),var1(2),var2(2));
    Zop.R22{3,2} = var_swap(var_swap(Zop.R22{2,3}',var1(1),var2(1)),var1(2),var2(2));
else
    if use_sep(3)
        Zop.R22{3,1} = Zop.R22{2,1};
    elseif ~excludeL{4,4}(3,1)
        [prog,R22_20] = sospolymatrixvar(prog,Z2{3,1},[m2,n2]);
        Zop.R22{3,1} = R22_20;
    end
    if use_sep(4)
        Zop.R22{1,3} = Zop.R22{1,2};
    elseif ~excludeL{4,4}(1,3)
        [prog,R22_02] = sospolymatrixvar(prog,Z2{1,3},[m2,n2]);
        Zop.R22{1,3} = R22_02;
    end
    if use_sep(5)
        Zop.R22{3,2} = Zop.R22{2,2};
        Zop.R22{3,3} = Zop.R22{2,3};
    else
        if ~excludeL{4,4}(3,2)
            [prog,R22_21] = sospolymatrixvar(prog,Z2{3,2},[m2,n2]);
            Zop.R22{3,2} = R22_21;
        end
        if ~excludeL{4,4}(3,3)
            [prog,R22_22] = sospolymatrixvar(prog,Z2{3,3},[m2,n2]);
            Zop.R22{3,3} = R22_22;
        end
    end
end

end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function [Zx,Zy,Z2] = build_monomials_2D(vars,dx,dy,d2)
% A subroutine to build monomial vectors in variables "vars" of the desired
% degrees dx, dy and d2. 

ss1 = vars(1,1);        tt1 = vars(1,2);
ss2 = vars(2,1);        tt2 = vars(2,2);

% Construct the monomials based on the specified maximal degrees
% Constructing Zxo(ss1)
Zxo = build_monoms(ss1.varname,dx{1});

% Constructing Zxa(ss1,tt1) and Zxb(ss1,tt1)
% In this implementation, Zxa will have degree dx{2}(2) in tt1 and degree
% dx{2}(1) in ss1 and max degree of ss1+tt1 is dx{2}(3). Similarly for Zxb(ss1,tt1)
st_varname = [ss1.varname; tt1.varname];
Zxa = build_monoms(st_varname,dx{2});
Zxb = build_monoms(st_varname,dx{3});

% % % % %

% Constructing Zyo(ss2)
Zyo = build_monoms(ss2.varname,dy{1});

% Constructing Zya(ss2,tt2) and Zyb(ss2,tt2)
% In this implementation, Zya will have degree dy{2}(2) in tt2 and degree
% dy{2}(1) in ss2 and max degree of ss2+tt2 is dy{2}(3). Similarly for Zyb(ss2,tt2)
st_varname = [ss2.varname; tt2.varname];
Zya = build_monoms(st_varname,dy{2});
Zyb = build_monoms(st_varname,dy{3});

% % % % %

% Constructing Z2oo(ss1,ss2)
ss_varname = [ss1.varname; ss2.varname];
Z2oo = build_monoms(ss_varname,d2{1,1});

% Constructing Z2ao(ss1,ss2,tt1) and Z2bo(ss1,ss2,tt1)
sst_varname = [ss1.varname; tt1.varname; ss2.varname];
Z2ao = build_monoms(sst_varname,d2{2,1});
Z2bo = build_monoms(sst_varname,d2{3,1});

% Constructing Z2oa(ss1,ss2,tt2) and Z2ob(ss1,ss2,tt2)
sst_varname = [ss1.varname; ss2.varname; tt2.varname];
Z2oa = build_monoms(sst_varname,d2{1,2});
Z2ob = build_monoms(sst_varname,d2{1,3});

% Constructing Z2aa(ss1,ss2,tt1,tt2) ... Z2bb(ss1,ss2,tt1,tt2)
sstt_varname = [ss1.varname; tt1.varname; ss2.varname; tt2.varname];
Z2aa = build_monoms(sstt_varname,d2{2,2});
Z2ba = build_monoms(sstt_varname,d2{3,2});
Z2ab = build_monoms(sstt_varname,d2{2,3});
Z2bb = build_monoms(sstt_varname,d2{3,3});

Zx = {Zxo; Zxa; Zxb};
Zy = {Zyo; Zya; Zyb};
Z2 = {Z2oo, Z2oa, Z2ob; Z2ao, Z2aa, Z2ab; Z2bo, Z2ba, Z2bb};

end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function Z = build_monoms(vartable,maxdegs)
% This function builds monomials in the variables "vartable" up to 
% maximal degrees "maxdegs". 
%
% INPUT
%   vartable: A (nvars x 1)-cellstr object specifying the variable names
%   maxdegs: An array of 2^nvars elements, for nvars variables. Required to
%            be a (reshaped version of a) 2x2x...x2 array, with each
%            dimension corresponding to a single variable, so that e.g.
%            element (2,1) corresponds to the degree of ss1, element (1,2)
%            corresponds to the degree of ss2, and element (2,2)
%            corresponds to the joint degree in ss1*ss2. Note that the
%            standard degree object used in poslpivar, taking the form
%
%                 0 |         ss2 |         tt2 |         ss2*tt2
%          ---------|-------------|-------------|-----------------
%               ss1 |     ss1*ss2 |     ss1*tt2 |     ss1*ss2*tt2
%          ---------|-------------|-------------|-----------------
%               tt1 |     tt1*ss2 |     tt1*tt2 |     tt1*ss2*tt2 
%          ---------|-------------|-------------|-----------------
%           ss1*tt1 | ss1*tt1*ss2 | ss1*tt1*tt2 | ss1*tt1*ss2*tt2  
%
%            can be easily reshaped to the 2x2x2x2 form. Also note that the
%            first element must always be zero (anything else is ignored).
%
% OUTPUT
%       Z: Polynomial class object corresponding to the monomials in the 
%           desired variables up to the desired maximal degrees
%

nvars = length(vartable);
if numel(maxdegs)==2^nvars-1
    maxdegs = [0;maxdegs(:)];
elseif numel(maxdegs)~=2^nvars
    error('The "maxdegs" input should contain 2^nvars elements')
end

Mdeg = max(maxdegs(:));

% Initialize
Z = monomials(vartable,0:Mdeg);
if ischar(vartable)
    Z_degmat = Z;
    Zvarname = {vartable};
elseif iscellstr(vartable)
    Z_degmat = Z;
    Zvarname = vartable;
elseif isa(vartable,'polynomial')
    Z_degmat = Z.degmat;
    Zvarname = Z.varname;
else
    error('Variable names should be specified as a cellstr or polynomial')
end

% Loop over all specified maximal (joint) degrees, and discard monomials
% that exceed the specified maximum
for i=2:2^nvars
    bindx = str2num(dec2bin(i-1,nvars)')';  
    cindx = bindx(end:-1:1)>0;                      % Logical array indicating which columns in degmat to consider
    retain = sum(Z_degmat(:,cindx),2)<=maxdegs(i);  % Logical array indicating rows with sufficiently small (joint) degree
    Z_degmat = Z_degmat(retain,:);                  % Retain only those monomials
end
% Get rid of variables that do not contribute
keep_var = any(Z_degmat,1);
Z_degmat = Z_degmat(:,keep_var);
Zvarname = Zvarname(keep_var);

% Build new monomials with maximal degrees as specified
nZ = size(Z_degmat,1);
Z = polynomial(eye(nZ),Z_degmat,Zvarname,[nZ,1]);

end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function [maxdegs,n_updates] = reduce_joint_degs(maxdegs)
% This function reduces the joint degrees presented in maxdegs to sensible
% values, making sure that the individual degrees do not exceed the joint
% degrees, nor that the sum of the individual degrees exceed the joint
% degrees.
%
% INPUT
%   maxdegs: An array of 2^nvars elements, for nvars variables. Required to
%            be a (reshaped version of a) 2x2x...x2 array, with each
%            dimension corresponding to a single variable, so that e.g.
%            element (2,1) corresponds to the degree of ss1, element (1,2)
%            corresponds to the degree of ss2, and element (2,2)
%            corresponds to the joint degree in ss1*ss2. Note that the
%            standard degree object used in e.g. poslpivar, taking the form
%
%                 0 |         ss2 |         tt2 |         ss2*tt2
%          ---------|-------------|-------------|-----------------
%               ss1 |     ss1*ss2 |     ss1*tt2 |     ss1*ss2*tt2
%          ---------|-------------|-------------|-----------------
%               tt1 |     tt1*ss2 |     tt1*tt2 |     tt1*ss2*tt2 
%          ---------|-------------|-------------|-----------------
%           ss1*tt1 | ss1*tt1*ss2 | ss1*tt1*tt2 | ss1*tt1*ss2*tt2  
%
%            can be easily reshaped to the 2x2x2x2 form. Also note that the
%            first element must always be zero (anything else is ignored).
%
% OUTPUT
%   maxdegs: An array of the same dimensions as the input, where now the
%            degrees are reduced to such an extent that:
%           - The joint degree in any set of variables (x_1,...,x_n) does
%             not exceed the joint degree in the set
%             (x_1,...,x_{i-1},x_{i+1},...,x_n) for any i\in\{1,...,n\};
%           - The joint degree in any set of variables (x_1,...,x_n) does
%             not exceed the sum of the degrees in x_i and in the set
%             (x_1,...,x_{i-1},x_{i+1},...,x_n) for any i\in\{1,...,n\};
%            Note that degrees are only reduced, never increased.
%   n_updates: Number of cycles of updates performed to obtain the final
%              maxdegs output. If the joint degrees in the input maxdegs
%              are already sensible, n_updates will be zero. Otherwise, it
%              should be just 1, but there may be particular (edge) cases
%              for which more than 1 cycle is necessary. Please share such
%              cases with us. The number of cycles is limited to 5.

% % Check how many variables we're considering
if all(size(maxdegs)==2)
    nvars = ndims(maxdegs);
else
    ndegs = numel(maxdegs);
    nvars = round(log(ndegs)/log(2));
    if ndegs==2^nvars-1
        maxdegs = [0;maxdegs(:)];
    elseif ndegs~=2^nvars
        error('Number of degrees in the input should be 2^nvars for nvars variables')
    end
end
dsize_prod = 2.^(0:nvars-1);    % Vector translating step 1 increase in dimension d to associated increase in linear index

% if maxdegs(1)~=0 && ~isnan(maxdegs(1))
% %    warning('The very first element of the maximal degree object does not correspond to any variable; the input value will be ignored')
%     maxdegs(1) = 0;
% end

% % To be safe, we perform several cycles of updates to the degrees, though
% % 1 should be sufficient (in general)
maxdegs_old = maxdegs;
max_updates = 5;
n_updates = max_updates;
for m=1:max_updates
% % Iterate over each of the (joint) degrees, making sure the joint degree
% % does not decrease upon including an additional variable, nor increase
% % with more than the maximal degree of this variable
for k=2^nvars-1:-1:2
    b_indx = str2num(dec2bin(k-1,nvars)')';  
    v_indx = b_indx(end:-1:1)>0;     % Logical array indicating which variables contribute to this (joint) degree
    for v=find(~v_indx)     % For each variable that does not contribute...
        i_indx = dsize_prod(v) + 1; % Establish the linear index associated to this individual variable      
        j_indx = dsize_prod(v) + k; % Establish the linear index associated to the joint degree with this variable
        
        % Make sure the joint degree without the variable does not exceed that with the variable
        maxdegs(k) = min(maxdegs(k),maxdegs(j_indx));
        
        % Make sure the joint degree with the variable does not exceed that 
        % of the sum of the degree without the variable and that of the 
        % individual variable
        maxdegs(j_indx) = min(maxdegs(j_indx),maxdegs(k)+maxdegs(i_indx)); 
    end
end
if ~any(maxdegs_old(:)-maxdegs(:))
    n_updates = m-1;
    break
elseif m==max_updates
    warning(['More than the allowed ',num2str(max_updates),' updates to the degree structure are necessary in reducing the joint degrees; something might be wrong...'])
else
    maxdegs_old = maxdegs;
end
end

end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function [dx,dy,d2] = check_degs(d)
% A subroutine to check that the degrees d of the monomials are
% appropriately specified.

if isnumeric(d) && all(size(d)==1)
    % Only 1 degree is specified, use this degree for all monomials.
    dval = d;
    dx = {dval;[dval;dval;dval];[dval;dval;dval]};       
    dy = {dval,[dval,dval,dval],[dval,dval,dval]};
    d2 = {dval*[0,1;1,2], dval*[0,1,1,1;1,1,1,2], dval*[0,1,1,1;1,1,1,2];
          dval*[0,1,1,1;1,1,1,2]', dval*ones(4,4), dval*ones(4,4);
          dval*[0,1,1,1;1,1,1,2]', dval*ones(4,4), dval*ones(4,4)};
elseif isa(d,'cell') && numel(d)==3
    % Degree is specified in 1D format.
    % --> fill in the gaps as per 1D approach
    if isscalar(d(:))
        d{2}=[d{1},d{1},2*d{1}];
        d{3}=d{2};
    else
        if isscalar(d{2})
            d{2} = d{2}*ones(1,3);
        elseif length(d{2})==2
            d{2}(3) = max(d{2});
        end
        if numel(d)==2
            d{3}=d{2};
        else
            if isscalar(d{3})
                d{3} = d{3}*ones(1,3);
            elseif length(d{3})==2
                d{3}(3) = max(d{3});
            end
        end
    end
    % Use specified 1D degrees along both variables.
    dx = d';        dy = d';
    d2{1,1} = [0,d{1};d{1},2*d{1}];
    d2{1,2} = [0,reshape(d{2},1,[]);d{1},d{1}+reshape(d{2},1,[])];
    d2{1,3} = d2{1,2};
    d2{2,1} = d2{1,2}';
    d2{3,1} = d2{1,3}';
    d2{2,2} = [0, reshape(d{2},1,[]); reshape(d{2},[],1), reshape(d{2},1,[])+reshape(d{2},[],1)];
elseif isa(d,'struct')
    if ~isfield(d,'dx') || ~isfield(d,'dy') || ~isfield(d,'d2')
        error("Degrees for 2D positive operator should be specified as struct with fields 'dx', 'dy', and 'd2'. Call 'help lpivar_2d' for more details.")
    else
        dx = d.dx;
        dy = d.dy;
        d2 = d.d2;
    end
end

% At this point, dx, dy and d2 should be cells.
if ~iscell(dx)
    error('Input d.dx must be a 3x1-cell structure')
end
if ~iscell(dy)
    error('Input d.dy must be a 1x3-cell structure')
end
if ~iscell(d2)
    error('Input d.d2 must be a 3x3-cell structure')
end
% % Check if degrees for x monomials are properly specified
if isscalar(dx(:))
    dx{2}=[dx{1};dx{1};2*dx{1}];
    dx{3}=dx{2};
elseif length(dx(:))==2
    if isscalar(dx{2})
        dx{2}(2)=dx{2}(1);
        dx{2}(3)=dx{2}(1);
    elseif length(dx{2})==2
        dx{2}(3) = ceil(0.5*(dx{2}(1) + dx{2}(2)));
    end
    dx{3}=dx{2};
else
    if isscalar(dx{2})
        dx{2}(2)=dx{2}(1);
        dx{2}(3)=dx{2}(1);
    elseif length(dx{2})==2
        dx{2}(3) = ceil(0.5*(dx{2}(1) + dx{2}(2)));
    end
    if isscalar(dx{3})
        dx{3}(2)=dx{3}(1);
        dx{3}(3)=dx{3}(1);
    elseif length(dx{3})==2
        dx{3}(3) = ceil(0.5*(dx{3}(1) + dx{3}(2)));
    end
end
% % Check if degrees for y monomials are properly specified
if isscalar(dy(:))
    dy{2}=[dy{1},dy{1},2*dy{1}];
    dy{3}=dy{2};
elseif length(dy(:))==2
    if isscalar(dy{2})
        dy{2}(2)=dy{2}(1);
        dy{2}(3)=dy{2}(1);
    elseif length(dy{2})==2
        dy{2}(3) = ceil(0.5*(dy{2}(1) + dy{2}(2)));
    end
    dy{3}=dy{2};
else
    if isscalar(dy{2})
        dy{2}(2)=dy{2}(1);
        dy{2}(3)=dy{2}(1);
    elseif length(dy{2})==2
        dy{2}(3) = ceil(0.5*(dy{2}(1) + dy{2}(2)));
    end
    if isscalar(dy{3})
        dy{3}(2)=dy{3}(1);
        dy{3}(3)=dy{3}(1);
    elseif length(dy{3})==2
        dy{3}(3) = ceil(0.5*(dy{3}(1) + dy{3}(2)));
    end
end
% % Check if degrees for 2D monomials are properly specified
if ~all(size(d2)==[3,3])
    if isscalar(d2) && all(size(d2{1,1})>=[2,2])
        d2 = repmat(d2,[3,3]);
    elseif all(size(d2)==[2,2]) && all(size(d2{2,1})>=[2,2]) && all(size(d2{1,2})>=[2,2]) && all(size(d2{2,2})>=[2,2])
        d2 = [d2,d2(:,2); d2(2,:),d2(2,2)];
    else
        error('Input ''d.d2'' must be a 3x3 cell');
    end
else
    % Check if sufficient degrees for R22oo(x,y) are specified
    if ~all(size(d2{1,1})>=[2,2])    % degrees in x and y
        error('Maximal degrees ''d.d2{1,1}'' must be specified using a 2x2 or 4x4 array');
    end
    % Check if sufficient degrees for R22ao(x,y,tt), R22bo(x,y,tt) are specified
    if ~all(size(d2{2,1})>=[4,2]) || ~all(size(d2{3,1})>=[4,2])
        for i=2:3
            if all(size(d2{i,1})==[2,2])    % only x and y degrees are specified
                d2{i,1} = [d2{i,1};d2{i,1}(2,:);d2{i,1}(2,:)];
            elseif all(size(d2{i,1})==[3,2]) % no joint degrees are specified
                d2{i,1} = [d2{i,1};sum(d2{i,1}(2:3,:),1)];
            else
                error(['Maximal degrees ''d.d2{',num2str(i),',1}'' must be specified using a 4x2 or 4x4 array']);
            end
        end
    end
    % Check if sufficient degrees for R22oa(x,y,nu), R22ob(x,y,nu) are specified
    if ~all(size(d2{2,1})>=[4,2]) || ~all(size(d2{3,1})>=[4,2])
        for j=2:3
            if all(size(d2{1,j})==[2,2])    % only x and y degrees are specified
                d2{1,j} = [d2{1,j},d2{1,j}(:,2),d2{1,j}(:,2)];
            elseif all(size(d2{1,j})==[2,3]) % no joint degrees are specified
                d2{1,j} = [d2{1,j},sum(d2{1,j}(:,2:3),2)];
            else
                error(['Maximal degrees ''d.d2{1,',num2str(j),'}'' must be specified using a 2x4 or 4x4 array']);
            end
        end
    end
    % Check if sufficient degrees for R22oa(x,y,nu), R22ob(x,y,nu) are specified
    if ~all(size(d2{2,1})>=[4,2]) || ~all(size(d2{3,1})>=[4,2])
        for i=2:3
          for j=2:3
            if all(size(d2{i,j})==[2,2])    % only x and y degrees are specified
                d2{i,j} = [d2{i,j},repmat(d2{i,j}(:,2),[1,2]);repmat(d2{i,j}(2,:),[2,1]),repmat(d2{i,j}(2,2),[2,2])];
            elseif all(size(d2{i,j})==[3,3]) % no joint degrees with secondary vars are specified
                d2{i,j} = [d2{i,j},sum(d2{i,j}(:,2:3),2);sum(d2{i,j}(2:3,:),1),d2{i,j}(3,3)];
            else
                error(['Maximal degrees ''d.d2{',num2str(i),',',num2str(j),'}'' must be specified using a 2x4 or 4x4 array']);
            end
          end
        end
    end
end
for k=1:9
    [d2{k},n_updates_k] = reduce_joint_degs(d2{k});
    % if n_updates_k>=1
    %     warning(['At least one of the (joint) degrees in d.d2{',num2str(k),'} is either too large or too small to make sense with all the joint degrees; reducing degrees to a sensible value'])
    % end
end
% In case insufficient degrees are provided, process the degrees
%[dx,dy,d2] = process_degrees(dx,dy,d2);

end