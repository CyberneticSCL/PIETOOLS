function [prog,Zop] = lpivar_2d(prog,dim,I,d,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Zop] = lpivar_2d(prog,dim,I,d,options) declares an indefinite
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
%   prog: SOS program to modify.
%   dim:  4x2 array providing the dimensions of the PI operators
%   I:    2x2 array [a,b;c,d] providing the domain ss1\in[a,b], ss2\in[c,d]
%
%   d, structure with fields dx, dy, d2, where:
%   dx{1}: degree of s in Zx^o(s)
%   dx{2}(1): degree of s in Zx^a(s,th), defaults to d(1)
%   dx{2}(2): degree of th in Zx^a(s,th), defaults to d(1)
%   dx{2}(3): joint degree of s,th in Zx^a(s,th), defaults to d(2,1)+d(2,2)
%   dx{3}(1): degree of s in Zx^b(s,th), defaults to d(1)
%   dx{3}(2): degree of th in Zx^b(s,th), defaults to d(1)
%   dx{3}(3): joint degree of s,th in Zx^b(s,th), defaults to d(3,1)+d(3,2)
%
%   d2{1,1}(1): degree of s1 in Z2^oo(s1,s2)
%   d2{1,1}(2): degree of s2 in Z2^oo(s1,s2)
%   d2{1,1}(3): joint degree of s1,s2 in Z2^oo(s1,s2), defaults to d{1,1}(1)+d{1,1}(2)
%
%   d2{2,1}(1): degree of s1 in Z2^ao(s1,th1,s2), defaults to d2{1,1}(1)
%   d2{2,1}(2): degree of th1 in Z2^ao(s1,th1,s2), defaults to d2{1,1}(1)
%   d2{2,1}(3): joint degree of s1,th1 in Z2^ao(s1,th1,s2), defaults to d{2,1}(1)+d{2,1}(2)
%   d2{2,1}(4): degree of s2 in Z2^ao(s1,th1,s2), defaults to d2{1,1}(2)
%   d2{2,1}(5): joint degree of s1,th1,s2 in Z2^ao(s1,th1,s2), defaults to d{2,1}(3)+d{2,1}(4)
%   d2{3,1}(1): degree of s1 in Z2^bo(s1,th1,s2), defaults to d2{1,1}(1)
%   d2{3,1}(2): degree of th1 in Z2^bo(s1,th1,s2), defaults to d2{1,1}(1)
%   d2{3,1}(3): joint degree of s1,th1 in Z2^bo(s1,th1,s2), defaults to d{3,1}(1)+d{3,1}(2)
%   d2{3,1}(4): degree of s2 in Z2^bo(s1,th1,s2), defaults to d2{1,1}(2)
%   d2{3,1}(5): joint degree of s1,th1,s2 in Z2^bo(s1,th1,s2), defaults to d{3,1}(3)+d{3,1}(4)
%
%   d2{1,2}(1): degree of s1 in Z2^oa(s1,s2,th2), defaults to d2{1,1}(1)
%   d2{1,2}(2): degree of s2 in Z2^oa(s1,s2,th2), defaults to d2{1,1}(2)
%   d2{1,2}(3): degree of th2 in Z2^oa(s1,s2,th2), defaults to d2{1,1}(2)
%   d2{1,2}(4): joint degree of s2,th2 in Z2^oa(s1,s2,th2), defaults to d2{1,2}(2)+d{1,2}(3)
%   d2{1,3}(1): degree of s1 in Z2^ob(s1,s2,th2), defaults to d2{1,1}(1)
%   d2{1,3}(2): degree of s2 in Z2^ob(s1,s2,th2), defaults to d2{1,1}(2)
%   d2{1,3}(3): degree of th2 in Z2^ob(s1,s2,th2), defaults to d2{1,1}(2)
%   d2{1,3}(4): joint degree of s2,th2 in Z2^ob(s1,s2,th2), defaults to d2{1,3}(2)+d{1,3}(3)
%
%   d2{2,2}(1,1): degree of s1 in Z2^aa(s1,th1,s2,th2), defaults to d2{1,1}(1)
%   d2{2,2}(2,1): degree of th1 in Z2^aa(s1,th1,s2,th2), defaults to d2{1,1}(1)
%   d2{2,2}(3,1): joint degree of s1,th1 in Z2^aa(s1,th1,s2,th2), defaults to d{2,2}(1,1)+d{2,2}(2,1)
%   d2{2,2}(1,2): degree of s2 in Z2^aa(s1,th1,s2,th2), defaults to d2{1,1}(2)
%   d2{2,2}(2,2): degree of th2 in Z2^aa(s1,th1,s2,th2), defaults to d2{1,1}(2)
%   d2{2,2}(3,2): joint degree of s2,th2 in Z2^aa(s1,th1,s2,th2), defaults to d{2,2}(1,2)+d{2,2}(2,2)
%
%   options.isscalar: Binary value set to 1 if the PI operator parameters
%                       should not vary in space (no dependence on s1,s2,
%                         th1,th2)
%
%   options.ismultiplier: Binary value set to 1 if the PI operator should
%                       not include any integrator terms
%
%   options.sep: Length 5 binary vector to enforce separability:
%      options.sep(1) = 1 if Rxx{2} = Rxx{3}, Rx2{2} = Rx2{3}, R2x{2} = Rx2{3},
%      options.sep(2) = 1 if Ryy{2} = Ryy{3}, Ry2{2} = Ry2{3}, R2y{2} = Ry2{3},
%      options.sep(3) = 1 if R22{2,1} = R22{3,1}
%      options.sep(4) = 1 if R22{1,2} = R22{1,3}
%      options.sep(5) = 1 if R22{2,2} = R22{3,2} = R22{2,3} = R22{3,3}
% 
% OUTPUT 
%   prog: modified SOS program with new variables and constraints
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


% % Extract the input arguments
switch nargin
    case 2
        error('Not enough inputs!')
    case 3
        dx = {1;[1;1;1];[1;1;1]};
        dy = {1,[1;1;1],[1;1;1]};
        d2 = {[1;1;2],[1;1;1;1;2],[1;1;1;1;2];
            [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1];
            [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1]};
        %options.exclude = zeros(1,16);
        %options.diag = 0;
        options.ismultiplier = 0;
        options.sep = zeros(1,6);
    case 4
        if ~isfield(d,'dx') && ~isfield(d,'dy') && ~isfield(d,'d2')
            fprintf('\n Warning: No degrees are specified. Continuing with default values. \n')
            dx = {1;[1;1;1];[1;1;1]};
            dy = {1,[1;1;1],[1;1;1]};
            d2 = {[1;1;2],[1;1;1;1;2],[1;1;1;1;2];
                [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1];
                [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1]};
        else
            if isfield(d,'dx')
                dx = d.dx;
            elseif isfield(d,'dy')
                dy = d.dy;
                dx = dy;
            else
                d2 = d.d2;
                dx = cell(3,1);
                dx{1} = d2{1,1}(1);  dx{2} = d2{2,1}(1:3);  dx{3} = d2{3,1}(1:3);
            end
            if isfield(d,'dy')
                dy = d.dy;
            elseif isfield(d,'dx')
                dy = dx;
            else
                d2 = d.d2;
                dy = cell(1,3);
                dy{1} = d2{1,1}(2);  dy{2} = d2{1,2}(2:4);  dy{3} = d2{1,3}(2:4);
            end
            if isfield(d,'d2')
                d2 = d.d2;
            else
                d2 = {[dx{1};dy{1}],[dx{1};dy{2}(:);ceil(0.5*(dx{1}+dy{2}(end)))],[dx{1};dy{3}(:);ceil(0.5*(dx{1}+dy{3}(end)))];
                    [dx{2}(:);dy{1};ceil(0.5*(dx{2}(end)+dy{1}))],[dx{2};dy{2};ceil(0.5*(dx{2}+dy{2}))],[dx{2};dy{3};ceil(0.5*(dx{2}+dy{3}))];
                    [dx{3}(:);dy{1};ceil(0.5*(dx{2}(end)+dy{1}))],[dx{3};dy{2};ceil(0.5*(dx{3}+dy{2}))],[dx{3};dy{3};ceil(0.5*(dx{3}+dy{3}))]};
            end
        end
        %options.exclude = zeros(1,16);
        %options.diag = 0;
        options.ismultiplier = 0;
        options.sep = zeros(1,6);
    case 5
        if ~isfield(d,'dx') && ~isfield(d,'dy') && ~isfield(d,'d2')
            fprintf('\n Warning: No degrees are specified. Continuing with default values. \n')
            dx = {1;[1;1;1];[1;1;1]};
            dy = {1,[1;1;1],[1;1;1]};
            d2 = {[1;1;2],[1;1;1;1;2],[1;1;1;1;2];
                [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1];
                [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1]};
        else
            if isfield(d,'dx')
                dx = d.dx;
            elseif isfield(d,'dy')
                dy = d.dy;
                dx = dy;
            else
                d2 = d.d2;
                dx = cell(3,1);
                dx{1} = d2{1,1}(1);  dx{2} = d2{2,1}(1:3);  dx{3} = d2{3,1}(1:3);
            end
            if isfield(d,'dy')
                dy = d.dy;
            elseif isfield(d,'dx')
                dy = dx;
            else
                d2 = d.d2;
                dy = cell(1,3);
                dy{1} = d2{1,1}(2);  dy{2} = d2{1,2}(2:4);  dy{3} = d2{1,3}(2:4);
            end
            if isfield(d,'d2')
                d2 = d.d2;
            else
                d2 = {[dx{1};dy{1}],[dx{1};dy{2}(:);ceil(0.5*(dx{1}+dy{2}(end)))],[dx{1};dy{3}(:);ceil(0.5*(dx{1}+dy{3}(end)))];
                    [dx{2}(:);dy{1};ceil(0.5*(dx{2}(end)+dy{1}))],[dx{2};dy{2};ceil(0.5*(dx{2}+dy{2}))],[dx{2};dy{3};ceil(0.5*(dx{2}+dy{3}))];
                    [dx{3}(:);dy{1};ceil(0.5*(dx{2}(end)+dy{1}))],[dx{3};dy{2};ceil(0.5*(dx{3}+dy{2}))],[dx{3};dy{3};ceil(0.5*(dx{3}+dy{3}))]};
            end
        end
        if isfield(options,'psatz')
            warning('There is no psatz option for indefinite PI operators.')
        end
        %         if ~isfield(options,'exclude')
        %             options.exclude = zeros(1,16);
        %         end
        if ~isfield(options,'isscalar')
            options.ismultiplier = 0;
        end
        if ~isfield(options,'ismultiplier')
            options.ismultiplier = 0;
        end
        %         if isfield(options,'diag') && options.diag==1
        %             fprintf(2,'Warning: ''diag'' option is not supported for 2D PDEs, ignoring this input.');
        %         end
        %         options.diag=0;
        if ~isfield(options,'sep')
            options.sep = zeros(1,6);
        end
end

if options.isscalar
    % In this case, we want the parameters to be uniform in space
    % --> monomial degrees 0.
    dx = {0,[0,0,0],[0,0,0]};
    dy = {0,[0,0,0],[0,0,0]};
    d2 = {zeros(2,1),zeros(5,1),zeros(5,1);
          zeros(5,1),zeros(3,3),zeros(3,3);
          zeros(5,1),zeros(3,3),zeros(3,3)};
else
    % Otherwise, check degrees are specified
    if ~iscell(dx)
        error('dx must be a 3-cell structure')
    end
    if ~iscell(dy)
        error('dy must be a 3-cell structure')
    end
    if ~iscell(d2)
        error('d2 must be a 3x3-cell structure')
    end
    % In case insufficient degrees are provided, process the degrees
    %[dx,dy,d2] = process_degrees(dx,dy,d2);
end


% Extract the size of the object: P\in R^(n0 x n0), Rxx\in L_2[x]^(nx x nx),
% Ryy\in L_2[y]^(ny x ny), R22\in L_2[x,y]^(n2 x n2),
if numel(dim)==1
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
    error('dim must be a 4x2 array')
end
m0 = m(1);    mx = m(2);    my = m(3);    m2 = m(4);
n0 = n(1);    nx = n(2);    ny = n(3);    n2 = n(4);
if n0==0 && nx==0 && ny==0 && n2==0
    error('All input dimensions are zero')
end
if m0==0 && mx==0 && my==0 && m2==0
    error('All output dimensions are zero')
end


% Specify the spatial domain (s,th \in [I(1),I(2)]
if any(size(I)~=[2,2])
    error('I must be a 2x2 array')
end
if any(I(:,1)>=I(:,2))
    error('I(:,1) must be less than I(:,2)')
end


% Define the primary variables
if isempty(prog.vartable)
    pvar ss1 ss2 tt1 tt2
    var1 = [ss1;ss2];       var2 = [tt1;tt2];
elseif length(prog.vartable)==2
    var1 = polynomial(prog.vartable);
    pvar tt1 tt2;
    var2 = [tt1;tt2];
elseif length(prog.vartable)==4
    var = polynomial(prog.vartable);
    var1 = var((1:2)');     var2 = var((3:4)');
else
    error('There should be two or four spatial variables in prog.vartable in a 2D LPI problem')
end
vars = [var1,var2];

% If we want a multiplier operator, the upper and lower diagonal integrals
% will both be zero
if options.ismultiplier==0
    options.sep = ones(1,5);
end


% % % Build the indefinite operator

% Construct monomials given the degrees
[Zx,Zy,Z2] = build_monomials_2D(vars,dx,dy,d2);


% % Define the parameters and add the decision variables to the program
% % structure

% Parameters mapping to R^m0
[prog,R00] = sospolymatrixvar(prog,1,[m0,n0]);
[prog,R0x] = sospolymatrixvar(prog,Zx{1},[m0,nx]);
[prog,R0y] = sospolymatrixvar(prog,Zy{1},[m0,ny]);
[prog,R02] = sospolymatrixvar(prog,Z2{1,1},[m0,n2]);

% Parameters mapping to L_2^mx[x]
[prog,Rx0] = sospolymatrixvar(prog,Zx{1},[mx,n0]);
[prog,Rxx_0] = sospolymatrixvar(prog,Zx{1},[mx,nx]);
if options.ismultiplier
    Rxx_1 = zeros([mx,nx]);
else
    [prog,Rxx_1] = sospolymatrixvar(prog,Zx{2},[mx,nx]);
end
if options.sep(1)
    Rxx_2 = Rxx_1;
else
    [prog,Rxx_2] = sospolymatrixvar(prog,Zx{3},[mx,nx]);
end
[prog,Rxy] = sospolymatrixvar(prog,Z2{1,1},[mx,ny]);
[prog,Rx2_0] = sospolymatrixvar(prog,Z2{1,1},[mx,n2]);
if options.ismultiplier
    Rx2_1 = zeros([mx,n2]);
else
    [prog,Rx2_1] = sospolymatrixvar(prog,Z2{2,1},[mx,n2]);
end
if options.sep(3)
    Rx2_2 = Rx2_1;
else
    [prog,Rx2_2] = sospolymatrixvar(prog,Z2{3,1},[mx,n2]);
end

% Parameters mapping to L_2^my[y]
[prog,Ry0] = sospolymatrixvar(prog,Zy{1},[my,n0]);
[prog,Ryx] = sospolymatrixvar(prog,Z2{1,1},[my,nx]);
[prog,Ryy_0] = sospolymatrixvar(prog,Zy{1},[my,ny]);
if options.ismultiplier
    Ryy_1 = zeros([my,ny]);
else
    [prog,Ryy_1] = sospolymatrixvar(prog,Zy{2},[my,ny]);
end
if options.sep(2)
    Ryy_2 = Ryy_1;
else
    [prog,Ryy_2] = sospolymatrixvar(prog,Zy{3},[my,ny]);
end
[prog,Ry2_0] = sospolymatrixvar(prog,Z2{1,1},[my,n2]);
if options.ismultiplier
    Ry2_1 = zeros([my,n2]);
else
    [prog,Ry2_1] = sospolymatrixvar(prog,Z2{1,2},[my,n2]);
end
if options.sep(4)
    Ry2_2 = Ry2_1;
else
    [prog,Ry2_2] = sospolymatrixvar(prog,Z2{1,3},[my,n2]);
end

% Parameters mapping to L_2^m2[x,y]
[prog,R20] = sospolymatrixvar(prog,Z2{1,1},[m2,n0]);
[prog,R2x_0] = sospolymatrixvar(prog,Z2{1,1},[m2,nx]);
if options.ismultiplier
    R2x_1 = zeros([m2,nx]);
else
    [prog,R2x_1] = sospolymatrixvar(prog,Z2{2,1},[m2,nx]);
end
if options.sep(3)
    R2x_2 = R2x_1;
else
    [prog,R2x_2] = sospolymatrixvar(prog,Z2{3,1},[m2,nx]);
end
[prog,R2y_0] = sospolymatrixvar(prog,Z2{1,1},[m2,ny]);
if options.ismultiplier
    R2y_1 = zeros([m2,ny]);
else
    [prog,R2y_1] = sospolymatrixvar(prog,Z2{1,2},[m2,ny]);
end
if options.sep(4)
    R2y_2 = R2y_1;
else
    [prog,R2y_2] = sospolymatrixvar(prog,Z2{1,3},[m2,ny]);
end
[prog,R22_00] = sospolymatrixvar(prog,Z2{1,1},[m2,n2]);
if options.ismultiplier
    R22_10 = zeros([m2,n2]);
else
    [prog,R22_10] = sospolymatrixvar(prog,Z2{2,1},[m2,n2]);
end
if options.sep(3)
    R22_20 = R22_10;
else
    [prog,R22_20] = sospolymatrixvar(prog,Z2{3,1},[m2,n2]);
end
if options.ismultiplier
    R22_01 = zeros([m2,n2]);
else
    [prog,R22_01] = sospolymatrixvar(prog,Z2{1,2},[m2,n2]);
end
if options.sep(4)
    R22_02 = R22_01;
else
    [prog,R22_02] = sospolymatrixvar(prog,Z2{1,3},[m2,n2]);
end
if options.ismultiplier
    R22_11 = zeros([m2,n2]);
else
    [prog,R22_11] = sospolymatrixvar(prog,Z2{2,2},[m2,n2]);
end
if options.sep(5)
    R22_12 = R22_11;
    R22_21 = R22_11;
    R22_22 = R22_11;
else
    [prog,R22_12] = sospolymatrixvar(prog,Z2{2,3},[m2,n2]);
    [prog,R22_21] = sospolymatrixvar(prog,Z2{3,2},[m2,n2]);
    [prog,R22_22] = sospolymatrixvar(prog,Z2{3,3},[m2,n2]);
end


% Construct the decision opvar2d object
Rxx = {Rxx_0; Rxx_1; Rxx_2};    Ryy = {Ryy_0, Ryy_1, Ryy_2};
Rx2 = {Rx2_0; Rx2_1; Rx2_2};    Ry2 = {Ry2_0, Ry2_1, Ry2_2};
R2x = {R2x_0; R2x_1; R2x_2};    R2y = {R2y_0, R2y_1, R2y_2};
R22 = {R22_00, R22_01, R22_02;
    R22_10, R22_11, R22_12;
    R22_20, R22_21, R22_22};

%Zop = dopvar2d([],[m,n],I,vars);
Zop = dopvar2d([],[m,n]);   Zop.I = I;
Zop.var1 = vars(:,1);       Zop.var2 = vars(:,2);

Zop.R00 = R00;  Zop.R0x = R0x;  Zop.R0y = R0y;  Zop.R02 = R02;
Zop.Rx0 = Rx0;  Zop.Rxx = Rxx;  Zop.Rxy = Rxy;  Zop.Rx2 = Rx2;
Zop.Ry0 = Ry0;  Zop.Ryx = Ryx;  Zop.Ryy = Ryy;  Zop.Ry2 = Ry2;
Zop.R20 = R20;  Zop.R2x = R2x;  Zop.R2y = R2y;  Zop.R22 = R22;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx,dy,d2] = process_degrees(dx,dy,d2)
% A subroutine to process the degrees d of the monomials specified by the
% user.

if length(dx(:))==1
    dx{2}=[dx{1},dx{1},2*dx{1}];
    dx{3}=dx{2};
elseif length(dx(:))==2
    if length(dx{2})==1
        dx{2}(2)=dx{2}(1);
        dx{2}(3)=dx{2}(1);
    elseif length(dx{2})==2
        dx{2}(3) = ceil(0.5*(dx{2}(1) + dx{2}(2)));
    end
    dx{3}=dx{2};
else
    if length(dx{2})==1
        dx{2}(2)=dx{2}(1);
        dx{2}(3)=dx{2}(1);
    elseif length(dx{2})==2
        dx{2}(3) = ceil(0.5*(dx{2}(1) + dx{2}(2)));
    end
    if length(dx{3})==1
        dx{3}(2)=dx{3}(1);
        dx{3}(3)=dx{3}(1);
    elseif length(dx{3})==2
        dx{3}(3) = ceil(0.5*(dx{3}(1) + dx{3}(2)));
    end
end

if length(dy(:))==1
    dy{2}=[dy{1},dy{1},2*dy{1}];
    dy{3}=dy{2};
elseif length(dy(:))==2
    if length(dy{2})==1
        dy{2}(2)=dy{2}(1);
        dy{2}(3)=dy{2}(1);
    elseif length(dy{2})==2
        dy{2}(3) = ceil(0.5*(dy{2}(1) + dy{2}(2)));
    end
    dy{3}=dy{2};
else
    if length(dy{2})==1
        dy{2}(2)=dy{2}(1);
        dy{2}(3)=dy{2}(1);
    elseif length(dy{2})==2
        dy{2}(3) = ceil(0.5*(dy{2}(1) + dy{2}(2)));
    end
    if length(dy{3})==1
        dy{3}(2)=dy{3}(1);
        dy{3}(3)=dy{3}(1);
    elseif length(dy{3})==2
        dy{3}(3) = ceil(0.5*(dy{3}(1) + dy{3}(2)));
    end
end

if size(d2,1)==1
    if length(d2{1,1})==1
        d2{1,1} = [d2{1,1};d2{1,1};d2{1,1}];
    elseif length(d2{1,1})==2
        d2{1,1} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(1)+d2{1,1}(2)];
    end
    
    d2{2,1} = [d2{1,1}(1);d2{1,1}(1);d2{1,1}(1);d2{1,1}(2);d2{1,1}(3)];
    d2{3,1} = d2{2,1};
        
    if size(d2,2)==1
        d2{1,2} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(2);d2{1,1}(2);d2{1,1}(3)];
        d2{1,3} = d2{1,2};
    elseif size(d2,2)==2
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        d2{1,3} = d2{1,2};
    else
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        if length(d2{1,3})==1
            d2{1,3}=[d2{1,1}(1);d2{1,3};d2{1,3};d2{1,3};ceil((d2{1,1}(1)+2*d2{1,3})/3)];
        elseif length(d2{1,3})==2
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(2);d2{1,3}(2);ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==3
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(3);ceil(0.5*(d2{1,3}(2)+d2{1,3}(3)));ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==4
            d2{1,3}=[d2{1,3}(1:4);ceil((d2{1,3}(1)+2*d2{1,3}(4))/3)];
        end
        
    end
    d2{2,2} = [d2{2,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,2}(2:4,1)))];
    d2{3,2} = [d2{3,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,2}(2:4,1)))];
    d2{2,3} = [d2{2,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,3}(2:4,1)))];
    d2{3,3} = [d2{3,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,3}(2:4,1)))];
    
    
elseif size(d2,1)==2
    if length(d2{1,1})==1
        d2{1,1}=[d2{1,1};d2{1,1};d2{1,1}];
    elseif length(d2{1,1})==2
        d2{1,1} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(1)+d2{1,1}(2)];
    end
    
    if length(d2{2,1})==1
        d2{2,1}=[d2{2,1};d2{2,1};d2{2,1};d2{2,1};ceil((2*d2{2,1}+d2{1,1}(2))/3)];
    elseif length(d2{2,1})==2
        d2{2,1}=[d2{2,1}(1);d2{2,1}(1);d2{2,1}(1);d2{2,1}(2);ceil(mean(d2{2,1}))];
    elseif length(d2{2,1})==3
        d2{2,1}=[d2{2,1}(1);d2{2,1}(2);ceil(0.5*(d2{2,1}(1)+d2{2,1}(2)));d2{2,1}(3);ceil(mean(d2{2,1}))];
    elseif length(d2{2,1})==4
        d2{2,1}=[d2{2,1}(1:4);ceil((d2{2,1}(4)+2*d2{2,1}(3))/3)];
    end
    d2{3,1} = d2{2,1};
    
    if size(d2,2)==1
        d2{1,2} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(2);d2{1,1}(2);d2{1,1}(3)];
        d2{1,3} = d2{1,2};
    elseif size(d2,2)==2
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        d2{1,3} = d2{1,2};
    else
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        if length(d2{1,3})==1
            d2{1,3}=[d2{1,1}(1);d2{1,3};d2{1,3};d2{1,3};ceil((d2{1,1}(1)+2*d2{1,3})/3)];
        elseif length(d2{1,3})==2
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(2);d2{1,3}(2);ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==3
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(3);ceil(0.5*(d2{1,3}(2)+d2{1,3}(3)));ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==4
            d2{1,3}=[d2{1,3}(1:4);ceil((d2{1,3}(1)+2*d2{1,3}(4))/3)];
        end
    end
        
    %d2{2,2} = [d2{2,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,2}(2:4,1)))];
    d2{3,2} = [d2{3,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,2}(2:4,1)))];
    %d2{2,3} = [d2{2,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,3}(2:4,1)))];
    d2{3,3} = [d2{3,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,3}(2:4,1)))];
    
else
    if length(d2{1,1})==1
        d2{1,1}=[d2{1,1};d2{1,1};d2{1,1}];
    elseif length(d2{1,1})==2
        d2{1,1} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(1)+d2{1,1}(2)];
    end
    
    if length(d2{2,1})==1
        d2{2,1}=[d2{2,1};d2{2,1};d2{2,1};d2{2,1};ceil((2*d2{2,1}+d2{1,1}(2))/3)];
    elseif length(d2{2,1})==2
        d2{2,1}=[d2{2,1}(1);d2{2,1}(1);d2{2,1}(1);d2{2,1}(2);ceil(mean(d2{2,1}))];
    elseif length(d2{2,1})==3
        d2{2,1}=[d2{2,1}(1);d2{2,1}(2);ceil(0.5*(d2{2,1}(1)+d2{2,1}(2)));d2{2,1}(3);ceil(mean(d2{2,1}))];
    elseif length(d2{2,1})==4
        d2{2,1}=[d2{2,1}(1:4);ceil((d2{2,1}(4)+2*d2{2,1}(3))/3)];
    end
    if length(d2{3,1})==1
        d2{3,1}=[d2{3,1};d2{3,1};d2{3,1};d2{3,1};ceil((2*d2{3,1}+d2{1,1}(2))/3)];
    elseif length(d2{3,1})==2
        d2{3,1}=[d2{3,1}(1);d2{3,1}(1);d2{3,1}(1);d2{3,1}(2);ceil(mean(d2{3,1}))];
    elseif length(d2{3,1})==3
        d2{3,1}=[d2{3,1}(1);d2{3,1}(2);ceil(0.5*(d2{3,1}(1)+d2{3,1}(2)));d2{3,1}(3);ceil(mean(d2{3,1}))];
    elseif length(d2{3,1})==4
        d2{3,1}=[d2{3,1}(1:4);ceil((d2{3,1}(4)+2*d2{3,1}(3))/3)];
    end
    
    if size(d2,2)==1
        d2{1,2} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(2);d2{1,1}(2);d2{1,1}(3)];
        d2{1,3} = d2{1,2};
    elseif size(d2,2)==2
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        d2{1,3} = d2{1,2};
    else
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        if length(d2{1,3})==1
            d2{1,3}=[d2{1,1}(1);d2{1,3};d2{1,3};d2{1,3};ceil((d2{1,1}(1)+2*d2{1,3})/3)];
        elseif length(d2{1,3})==2
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(2);d2{1,3}(2);ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==3
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(3);ceil(0.5*(d2{1,3}(2)+d2{1,3}(3)));ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==4
            d2{1,3}=[d2{1,3}(1:4);ceil((d2{1,3}(1)+2*d2{1,3}(4))/3)];
        end
    end    
    %d2{2,2} = [d2{2,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,2}(2:4,1)))];
    %d2{3,2} = [d2{3,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,2}(2:4,1)))];
    %d2{2,3} = [d2{2,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,3}(2:4,1)))];
    %d2{3,3} = [d2{3,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,3}(2:4,1)))];
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%%
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