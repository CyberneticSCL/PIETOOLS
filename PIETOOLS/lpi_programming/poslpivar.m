function [prog,Pop,Qmat,Zop,gs] = poslpivar(prog,n,d,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Pop,Qmat,Zop,gs] = poslpivar(prog,n,d,options) declares 
% a positive 4-PI decision operator Pop, taking the form
%
%       Pop = Zop'*((gs*Qmat)*Zop)
% 
% where
%
%   Qmat = [ Q_{11}  Q_{12} Q_{13} Q_{14}]
%          [ Q_{21}  Q_{22} Q_{23} Q_{24}] >=0
%          [ Q_{31}  Q_{32} Q_{33} Q_{34}]
%          [ Q_{41}  Q_{42} Q_{43} Q_{44}]
%
% and
%
% (Zop*[u0;u1])(s) = [u1; 
%                     Z1(s)*u2(s); 
%                     int_{a}^{s} Z2(s,th)*u2(th) dth
%                     int_{s}^{b} Z3(s,th)*u2(th) dth]
%                      
% so that, in particular,
%
%   Pop.P = Q11 int(gs,s,a,b)
%   Pop.Q(s) = g(s) Q12 Z1(s) + int( gth Q13 Z2(th,s) dth, s, b) + int( gth Q14 Z3(th,s) dth, a, s)
%   Pop.R0(s) = gs Z1(s) Q22 Z1(s)
%   Pop.R1(s,th) = gs Z1(s) Q23 Z2(s,th) + gth Z3(th,s) Q42 Z1(th) + 
%            int(geta Z2(eta,s)' Q33 Z2(eta,th) deta,s,b)+int(geta Z3(eta,s) Q43 Z2(eta,th) deta,th,s)
%           +int(geta Z3(eta,s) Q44 Z3(eta,th) deta,a,th)
%   Pop.R2(s,th) = R1(th,s)',
%
% where Zi(x)= Z_d1(x) \otimes I_n and Z_d(x) is the vector of monomials in
% variables x of degree d1 or less. Z(x,y) = Z_{d1}(x) \otimes Z_{d2}(y)
% \otimes I_n. If the application is stability of time-delay systems, d1
% will eventually be more or less the degree of the multiplier and d2 more
% or less the degree of the kernel function.
% 
% INPUT 
% - prog: 'struct' speicfying an LPI program to modify;
% - n:  2x1 array of type double specifying the dimensions of the operator.
%       Since the operator variable is symmetric, the input and output
%       dimensions will be the same. The first element n(1) should be the
%       dimension of the real part, and n(2) the dimension of the L2 part;
% - d:  1x3 cell array with each element specifying a maximal degree of
%       monomials appearing in the operator variables. Specifically:
%           d{1}: degree of s in Z1(s),                 defaults to 1
%           d{2}(1): degree of s in Z2(s,th),           defaults to d{1}
%           d{2}(2): degree of th in Z2(s,th),          defaults to d{2}(1)
%           d{2}(3): joint degree of s,th in Z2(s,th),  defaults to max(d{2}(1),d{2}(2))
%           d{3}(1): degree of s in Z3(s,th),           defaults to d{2}(1)
%           d{3}(2): degree of th in Z3(s,th),          defaults to d{2}(2)
%           d{3}(3): joint degree of s,th in Z3(s,th),  defaults to max(d{3}(1),d{3}(2))
%       If only a single degree d is specified as a scalar, a cell
%       structure will be generated as {d, [d,d,d], [d,d,d]};
% - options:    struct specifying other options for the operator, with
%               fields:
%   options.psatz   binary value, set to 1 to enforce positivity of the 
%                   operator only on the domain [a,b], using the function
%                   g(s) = (b-s)*(s-a) in the construction;
%   options.exclude 1x4 binary vector, set options.exclude(i)=1 to enforce
%                   T_{ij}=0 for j=1...4;
%   options.sep     binary value, set to 1 to enforce R1 = R2 ;
% 
% OUTPUT 
% - prog:   struct, specifying the same LPI program as the input, but now 
%           including the decision variables defining the operator Pop and
%           a constraint enforcing Pop>=0;
% - Pop:    nxn opvar representing a positive semidefinite PI operator
%           decision variable Pop = Zop'*(Qmat*gs)*Zop);
% - Qmat:   mxm object of type 'dpvar' representing the positive
%           semidefinite matrix-valued decision variable parameterizing
%           Qmat;
% - Zop:    mxn opvar object with all parameters given by monomial vectors,
%           representing the monomial operator;
% - gs:     scalar object of type 'polynomial', representing the function
%           g(s) used to enforce positivity only on the interval [a,b]. If
%           options.psatz = 1, g(s)=(s-a)*(b-s). Otherwise, g(s)=1;
% 
% NOTES:
% If the input optimization program structure "prog" corresponds to a 2D
% problem, the 2D version will be called automatically, and the returned
% operator will be an 'opvar2d' object. Call 'help poslpivar_2d' for more
% information on how to specify degrees and options for the 2D version.
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - poslpivar
%
% Copyright (C)2024 PIETOOLS Team
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
% Initial coding MMP, SS, DJ  - 09/26/2021
% DJ, 10/16/2024: Update to new 'lpiprogram' structure;
% DJ, 12/14/2024: Allow degrees to be specified as 'double' array;
% DJ, 01/23/2025: Check that 'options' are specified as struct;
% DJ, 01/25/2025: Return the monomial operator Zop and function gs;


% % % Set-up

% % First, check the spatial domain on which the program is defined.
if ~isfield(prog,'dom') || size(prog.dom,1)==0
    error('The program structure does not include a spatial domain -- please use ''lpiprogram'' to initialize your program');
else
    I = prog.dom;
end

% % Extract the input arguments
switch nargin
    case 1
        error('Not enough inputs!')
    case 2
        if all(size(I)==[2,2])
            if nargout<=2
                [prog,Pop] = poslpivar_2d(prog,n);
            elseif nargout==3
                [prog,Pop,Qmat] = poslpivar_2d(prog,n);
            elseif nargout==4
                [prog,Pop,Qmat,Zop] = poslpivar_2d(prog,n);
            else
                [prog,Pop,Qmat,Zop,gs] = poslpivar_2d(prog,n);
            end
            return
        end
        d = {1,[1,1,1],[1,1,1]};
        options.psatz=0;
        options.exclude=[0 0 0 0];
        options.sep =0;
    case 3
        if all(size(I)==[2,2])
            if nargout<=2
                [prog,Pop] = poslpivar_2d(prog,n,d);
            elseif nargout==3
                [prog,Pop,Qmat] = poslpivar_2d(prog,n,d);
            elseif nargout==4
                [prog,Pop,Qmat,Zop] = poslpivar_2d(prog,n,d);
            else
                [prog,Pop,Qmat,Zop,gs] = poslpivar_2d(prog,n,d);
            end
            return
        end
        options.psatz=0;
        options.exclude=[0 0 0 0];
        options.sep =0;
    case 4
        if all(size(I)==[2,2])
            if nargout<=2
                [prog,Pop] = poslpivar_2d(prog,n,d,options);
            elseif nargout==3
                [prog,Pop,Qmat] = poslpivar_2d(prog,n,d,options);
            elseif nargout==4
                [prog,Pop,Qmat,Zop] = poslpivar_2d(prog,n,d,options);
            else
                [prog,Pop,Qmat,Zop,gs] = poslpivar_2d(prog,n,d,options);
            end
            return
        end
        if ~isa(options,'struct')                                           % DJ, 01/23/2025
            error("Options for positive operator should be specified as struct with fields 'psatz', 'exclude', and 'sep'.")
        end
        if ~isfield(options,'psatz')
            options.psatz=0;
        end
        if ~isfield(options,'exclude')
            options.exclude=[0 0 0 0];
        end
        if ~isfield(options,'sep')
            options.sep=0;
        end
end

% % Check that the dimensions of the operator are properly specified.
if ~isnumeric(n) || any(any(n<0))
    error('Dimension of the operator must be specified as numeric array of positive integers.')
end
if numel(n)~=2
    if all(size(n)==[2,2])
        % Assume n(:,1) gives row dimension, and n(:,2) the column
        % dimension
        if any(n(:,1)~=n(:,2))
            error('Positive operator must have equal row and column dimensions.')
        else
            n = n(:,1);
        end
    else
        error('Dimension of the operator must be specified as 2x1 array.')
    end
end
% Extract the size of the object: 
%   Pop: \R^n1 xL2^n2 --> \R^n1 xL2^n2
n = n(:);
n1 = n(1);  n2 = n(2);
if all(n==0)
    % Construct an empty operator variable (for whatever reason...)
    dopvar Pop;
    Pop.I = I;
    Pop.var1 = prog.vartable(1);
    Pop.var2 = prog.vartable(2);
    return
end

% % Check that the degrees of the monomials have been properly specified
if isempty(d)
    % Use default degrees if none are specified.
    d = {1,[1,1,1],[1,1,1]};
end
if isnumeric(d)                                                             % DJ, 12/14/2024
    % Construct cell structure specifying the degrees.
    d = d(:)';
    if isscalar(d)
        d = {d};
    elseif numel(d)==2
        d = {max(d),[d,max(d)],[d,max(d)]};
    elseif numel(d)==3
        d = {max(d),d,d};
    else
        error('Degrees must be specified as 1x3 cell array.')
    end
elseif ~iscell(d)
    error('Degrees must be specified as 1x3 cell array.')
end
% Fill in any gaps in the data
if length(d(:))==1
    d{2}=[d{1},d{1},d{1}];
    d{3}=d{2};
else
    if length(d{2})==1
        d{2} = d{2}*ones(1,3);
    elseif length(d{2})==2
        d{2}(3) = max(d{2});
    end
    if numel(d)==2
        d{3}=d{2};
    elseif numel(d)==3
        if length(d{3})==1
            d{3} = d{3}*ones(1,3);
        elseif length(d{3})==2
            d{3}(3) = max(d{3});
        end
    else
        error('Degrees must be specified as 1x3 cell array.')
    end
end


% % To reduce complexity, allow certain terms to be excluded

% excludeL is a length-4 binary vector of terms to exclude
excludeL = options.exclude;

% In separable case R1=R2, so Qi3 = Qi4 and Q33 = Q44, so we may exclude
% Qi4 entirely, and use Qi3 instead.
if options.sep ==1 
    excludeL(4)=1;
end

% If there is no finite-dimensional component or infinite-dimensional
% component, we may exclude the associated term
if n1==0
    excludeL(1)=1;
end
if n2==0
    excludeL(2:4)=[1 1 1];
end
if all(excludeL)
    error('You''re creating an empty dopvar! Please change options.exclude and/or options.sep.')
end

% % This option is currently not available
if isfield(options,'diag') && options.diag
    disp('Poslpivar option ''diag'' is currently not available.')
end


% % Extract the variables defining the operator
%   R0 is function of var1, R1 and R2 of var1 and var2
var1 = prog.vartable(1);
var2 = prog.vartable(2);

% Define a dummy variable for integration (replaces eta)
pvar sss

% Define the multiplier function to be used later
if options.psatz==1 % use a function that is positive only on domain I
    gs = (var1-I(1))*(I(2)-var1);
    gth = (var2-I(1))*(I(2)-var2);
    geta = (sss-I(1))*(I(2)-sss);
else % use a function that is positive everywhere
    gs = polynomial(1);
    gth = polynomial(1);
    geta = polynomial(1);
end



% % % Construct the monomial vectors Z

% Z1(s) will have max degree d{1} in s.
Z1s = build_monoms(var1,[0;d{1}]);

% Z2(s,th) will have max degree d{2}(1) in s, d{2}(2) in th, and d{2}(3) in
% s and th combined.
Z2sth = build_monoms([var1;var2],reshape([0;d{2}(:)],[2,2]));

% Z3(s,th) will have max degree d{3}(1) in s, d{3}(2) in th, and d{3}(3) in
% s and th combined.
Z3sth = build_monoms([var1;var2],reshape([0;d{3}(:)],[2,2]));



% % % Build the matrix-valued function 
%   N(s,th,et) = ZL(s,et)'*Qmat*ZR(th,et)
% for Qmat>=0 and
%   ZL = [eye(n1), 0,                     0,                       0                      ]
%        [0,       kron(eye(n2), Z1(s)),  0,                       0                      ]
%        [0,       0,                     kron(eye(n2),Z2(et,s)),  0                      ]
%        [0,       0,                     0,                       kron(eye(n2),Z3(et,s)) ]
%   ZR = [eye(n1), 0,                     0,                       0                      ]
%        [0,       kron(eye(n2), Z1(th)), 0,                       0                      ]
%        [0,       0,                     kron(eye(n2),Z2(et,th)), 0                      ]
%        [0,       0,                     0,                       kron(eye(n2),Z3(et,th))]

% Evaluate monomials at dummy variable
Z1th = subs(Z1s,var1,var2);
Z2etath = subs(Z2sth,var1,sss);
Z2etas = subs(Z2etath,var2,var1);
Z3etath = subs(Z3sth,var1,sss);
Z3etas = subs(Z3etath,var2,var1);

% Initialize ZL and ZR
includeL = ~excludeL;
ZL = cell(1,sum(includeL));
ZR = cell(1,sum(includeL));

% Add desired monomial vectors to the block-diagonal matrix
mdim = [];      ndim = [];
Zdim = zeros(4,1);
indx = 1;
if includeL(1)
    ZL{indx} = 1;       ZR{indx} = 1;
    mdim = [mdim;n1];   ndim = [ndim;n1];
    indx = indx+1;
    Zdim(1) = n1;
end
if includeL(2)
    ZL{indx} = Z1s;     ZR{indx} = Z1th;
    mdim = [mdim;n2];   ndim = [ndim;n2];
    indx = indx+1;
    Zdim(2) = n2*size(Z1s,1);
end
if includeL(3)
    ZL{indx} = Z2etas;  ZR{indx} = Z2etath;
    mdim = [mdim;n2];   ndim = [ndim;n2];
    indx = indx+1;
    Zdim(3) = n2*size(Z2etas,1);
end
if includeL(4)
    ZL{indx} = Z3etas;  ZR{indx} = Z3etath;
    mdim = [mdim;n2];   ndim = [ndim;n2];
    Zdim(4) = n2*size(Z3etas,1);
end

% Declare the monomial operator Zop.
if nargout>=4
    opvar Zop;                                                              % DJ, 01/25/2025
    Zop.var1 = var1;    Zop.var2 = var2;
    Zop.I = I;
    Zop.dim = [0,n1;sum(Zdim),n2];
    Zop.Q2 = [eye(n1*includeL(1),n1); zeros(sum(Zdim(2:4)),n1)];
    Zop.R.R0 = [zeros(Zdim(1),n2); kron(eye(n2*includeL(2)),Z1s); zeros(sum(Zdim(3:4)),n2)];
    Zop.R.R1 = [zeros(sum(Zdim(1:2)),n2); kron(eye(n2*includeL(3)),Z2sth); zeros(Zdim(4),n2)];
    if options.sep
        Zop.R.R2 = [zeros(sum(Zdim(1:2)),n2); kron(eye(n2*includeL(3)),Z2sth)];
    else
        Zop.R.R2 = [zeros(sum(Zdim(1:3)),n2); kron(eye(n2*includeL(4)),Z3sth)];
    end
end


% Declare a matrix Qmat>=0, and matrix-valued function
%   N(s,th,et) = ZL(s,et)'*Qmat*ZR(th,et)
if nargout<=2
    [prog,N] = sosquadvar(prog,ZL,ZR,mdim,ndim,'pos');
else
    [prog,N,Qcell] = sosquadvar(prog,ZL,ZR,mdim,ndim,'pos');
    % Convert cell to 'dpvar' object.
    Qmat = dpvar([]);
    for ii=1:size(Qcell,1)
        Qmat_ii = dpvar([]);
        for jj=1:size(Qcell,2)
            Qmat_ii = [Qmat_ii,dpvar(Qcell{ii,jj})];
        end
        Qmat = [Qmat;Qmat_ii];
    end
end



% % % Build the positive dopvar object

% Since the object must be self-adjoint, we need only define 
% Pop.P = P, Pop.Q1 = QT, Pop.R.R0 = R0 and Pop.R.R1 = R1
P = dpvar(zeros(n1));
QT = dpvar(zeros(n1,n2));
R0 = dpvar(zeros(n2));
R1 = R0;

% Account for the fact that some components may be excluded.
include_indx = cumsum(includeL);

% Define P, Q1 = QT, and Q2 = QT'
if includeL(1)
    i1 = include_indx(1);     j1 = include_indx(1);
    P = P + N{i1,j1} * int(gs,var1,I(1),I(2));
    if includeL(2)
        i1 = include_indx(1);     j1 = include_indx(2);
        QT = QT + gs * subs(N{i1,j1},var2,var1);
    end
    if includeL(3) && ~options.sep==1
        i1 = include_indx(1);     j1 = include_indx(3);
        QT = QT + int(geta * subs(N{i1,j1},var2,var1),sss,var1,I(2));
    elseif includeL(3) && options.sep==1
        i1 = include_indx(1);     j1 = include_indx(3);
        QT = QT + int(geta * subs(N{i1,j1},var2,var1),sss,I(1),I(2));
    end
    if includeL(4)
        i1 = include_indx(1);     j1 = include_indx(4);
        QT = QT + int(geta * subs(N{i1,j1},var2,var1),sss,I(1),var1);
    end
end

% Define terms in R0 and R1 related to the first vector of monomials Z1
if includeL(2)
    i1 = include_indx(2);     j1 = include_indx(2);
    R0 = R0 + gs * subs(N{i1,j1},var2,var1);
    if includeL(3) && ~options.sep==1
        i1 = include_indx(2);     j1 = include_indx(3);
        R1 = R1 + gs * subs(N{i1,j1},sss,var1);
    elseif includeL(3) && options.sep==1
        i1 = include_indx(2);     j1 = include_indx(3);
        R1 = R1 + gs * subs(N{i1,j1},sss,var1) + gth * subs(N{j1,i1},sss,var2);
    end
    if includeL(4)
        i1 = include_indx(4);     j1 = include_indx(2);
        R1 = R1 + gth * subs(N{i1,j1},sss,var2);
    end
end

% Define terms in R1 related to the second vector of monomials Z2
if includeL(3) && ~options.sep==1
        i1 = include_indx(3);     j1 = include_indx(3);
        R1 = R1 + int(geta * N{i1,j1},sss,var1,I(2));
    if includeL(4)
        i1 = include_indx(4);     j1 = include_indx(3);
        R1 = R1 + int(geta * N{i1,j1},sss,var2,var1);
    end
elseif includeL(3) && options.sep==1
    i1 = include_indx(3);     j1 = include_indx(3);
    R1 = R1 + int(geta * N{i1,j1},sss,I(1),I(2));
end

% Define terms in R1 related to the third vector of monomials Z3
if includeL(4)
    i1 = include_indx(4);     j1 = include_indx(4);
    R1 = R1 + int(geta * N{i1,j1},sss,I(1),var2);
end

% Define R2
if options.sep==1
    R2 = R1;
else
    R2 = varswap(R1,var1,var2).';
end

% Construct the positive decision variable operator Pop
dopvar Pop;
Pop.P = P;
Pop.Q1 = QT;
Pop.Q2 = QT.'; 
Pop.R.R0 = R0; Pop.R.R1 = R1; Pop.R.R2 = R2;

Pop.dim = [n,n];
Pop.I = I;
Pop.var1 = var1;
Pop.var2 = var2;

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
degmat_sort = [sum(Z_degmat,2),fliplr(Z_degmat)];
degmat_sort = sortrows_integerTable(degmat_sort);       % sort monomials from low to high degree
Z_degmat = fliplr(degmat_sort(:,2:end));
Z = polynomial(eye(nZ),Z_degmat,Zvarname,[nZ,1]);

end