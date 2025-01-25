function [prog,Pop,Qmat,Zop,gss] = poslpivar_2d(prog,n,d,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Pop,Qmat,Zop,gss] = poslpivar_2d(prog,n,d,options) declares 
% a positive semidefinite 2D PI operator, 
%
%   Pop = Zop'*((Qmat*gss)*Zop)
%
% where
%
% Qmat = [ Q_{00,00} Q_{00,x1} Q_{00,x2} Q_{00,x3} Q_{00,y1} Q_{00,y2} Q_{00,y3} Q_{00,21} ... Q_{00,29}]
%        [      :         :         :         :         :         :         :         :             :   ] >= 0
%        [ Q_{29,00} Q_{29,x1} Q_{29,x2} Q_{29,x3} Q_{29,y1} Q_{29,y2} Q_{29,y3} Q_{29,21} ... Q_{29,29}]
%
% and
%
% Z*u = [                                 u0                ]
%       [                          Zxo(x)*ux(x)             ]
%       [               int_a^x Zxa(x,tt)*ux(tt) dtt        ]
%       [               int_x^b Zxb(x,tt)*ux(tt) dtt        ]
%       [                          Zyo(x)*uy(y)             ]
%       [               int_c^y Zya(y,nu)*uy(nu) dnu        ]
%       [               int_y^d Zyb(y,nu)*uy(nu) dnu        ]
%       [                       Z2oo(x,y)*u2(x,y)           ]
%       [            int_a^x Z2ao(x,y,tt)*u2(tt,y) dtt      ]
%       [            int_x^b Z2bo(x,y,tt)*u2(tt,y) dtt      ]
%       [            int_c^y Z2oa(x,y,nu)*u2(x,nu) dnu      ]
%       [            int_y^d Z2ob(x,y,nu)*u2(x,nu) dnu      ]
%       [ int_a^x int_c^y Z2aa(x,y,tt,nu)*u2(tt,nu) dnu dtt ]
%       [ int_a^x int_y^d Z2ab(x,y,tt,nu)*u2(tt,nu) dnu dtt ]
%       [ int_x^b int_c^y Z2ba(x,y,tt,nu)*u2(tt,nu) dnu dtt ]
%       [ int_x^b int_y^d Z2bb(x,y,tt,nu)*u2(tt,nu) dnu dtt ]
%
% where each function Z.. =  Z_d(x,y,tt,nu) \otimes I_n.. with Z_d(x) is a 
% vector of monomials in variables x,y,tt,ss of total degree d or less. 
% 
% INPUT 
%   prog: SOS program to modify.
%   n(1): dimension of real part
%   n(2): dimension of L2[x] part
%   n(3): dimension of L2[y] part
%   n(4): dimension of L2[x,y] part
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
%   options.psatz=1 if this is a psatz term. options.psatz=0 otherwise
%   options.exclude is a length 16 binary vector where 
%      options.exclude(i)=1 if we want to set $T_{ij}=0$ for j=1...16
%   options.sep is a length 6 binary vector where
%      options.sep(1) = 1 if Rxx{2} = Rxx{3}
%      options.sep(2) = 1 if Ryy{2} = Ryy{3}
%      options.sep(3) = 1 if R22{2,1} = R22{3,1}
%      options.sep(4) = 1 if R22{1,2} = R22{1,3}
%      options.sep(5) = 1 if R22{2,2} = R22{3,2} and R22{2,3} = R22{3,3}
%      options.sep(6) = 1 if R22{2,2} = R22{2,3} and R22{3,2} = R22{3,3}
% 
% OUTPUT 
%   prog:   struct, specifying the same LPI program as the input, but now 
%           including the decision variables defining the operator Pop and
%           a constraint enforcing Pop>=0;
%   Pop:    nxn 'opvar2d' object representing the positive semidefinite
%           PI operator decision variable;
% - Qmat:   mxm object of type 'dpvar' representing the positive
%           semidefinite matrix-valued decision variable parameterizing
%           Qmat;
% - Zop:    mxn opvar object with all parameters given by monomial vectors,
%           representing the monomial operator;
% - gss:    scalar object of type 'polynomial', representing the function
%           g(s) used to enforce positivity only on bounded domain. If
%           options.psatz = 1, gss(x,y)=(x-a)*(b-x)*(y-c)*(d-y). If
%           options.psatz = 2, gss(x,y)=(R^2 - (x-x_0)^2 - (y-y_0)^2),
%           where x_0=(b+a)/2, y_0=(d+c)/2, and R=norm([x-x_0;y-y_0]);
%           Otherwise, gss(x,y)=1;
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - poslpivar_2d
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
% Initial coding DJ - 09_28_2021
% DJ, 04/14/2022: Update to account for new degree data forma
% DJ, 12/15/2024: Bugfix in case no degrees are specified;
% DJ, 01/23/2025: Check that 'options' are specified as struct;
% DJ, 01/25/2025: Return the monomial operator Zop and function gs;

% % % Set-up % % %

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
use_monomials = 0;      % assume monomial degrees rather than actual monomials are specified
dx = {1;[1;1;1];[1;1;1]};
dy = {1,[1,1,1],[1,1,1]};
d2 = {[0,1;1,2],          [0,1,1,1;1,2,2,2], [0,1,1,1;1,2,2,2];
      [0,1,1,1;1,2,2,2]', [0,1,1,1;1,2,2,2;1,2,2,2;1,2,2,2], [0,1,1,1;1,2,2,2;1,2,2,2;1,2,2,2];
      [0,1,1,1;1,2,2,2]', [0,1,1,1;1,2,2,2;1,2,2,2;1,2,2,2], [0,1,1,1;1,2,2,2;1,2,2,2;1,2,2,2]};

% Initialize default exclude options
is_sep = ones(1,6);    % Use fully separable operator by default
% Exclude only single integrals along 2D states by default.
excludeL = [0, 0,0,0, 0,0,0, 0,1,1,1,1,0,0,0,0];
% Do not use psatz option by default
psatz = 0;
center_monomials = 0;

% Extract the input arguments
switch nargin
    case 1
        error('Not enough inputs!')
    case 2
        fprintf('\n Warning: No degrees are specified. Continuing with default values. \n')
    case 4
        if ~isa(options,'struct')                                           % DJ, 01/23/2025
            error("Options for positive operator should be specified as struct with fields 'psatz', 'exclude', and 'sep'.")
        end
        if isfield(options,'psatz')
            psatz = abs(options.psatz);
        end
        if isfield(options,'exclude')
            % Exclude only single integrals along 2D states by default.
            excludeL = options.exclude;
        end
        if isfield(options,'diag') && options.diag==1
            fprintf(2,'Warning: ''diag'' option is not supported for 2D PDEs, ignoring this input.'); 
        end
        options.diag=0;
        if isfield(options,'sep')
            is_sep = options.sep;
        end
        if isfield(options,'center_monomials')
            center_monomials = options.center_monomials;
        end
end

if nargout>5
    error('At most 5 outputs are produced');
end

% % % Check if inputs are properly specified
% % Start with the dimensions of the operator, should be a 4x1 or 4x2 array
if all(size(n)==[4,2])
    if all(n(:,1)==n(:,2))
        n = n(:,1)';
    else
        error('Positive operator must have equal row and column dimensions.')
    end
elseif numel(n)>4
    error('Operator dimensions n must be specified as 4x2 or 4x1 array')
end
n = n(:);
if length(n)==1
    warning('Only 1 dimension is provided, assuming this to refer to the 2D PDE state')
    n = [0;0;0;n(1)];
elseif length(n)==2
    warning('Only 2 dimensions are provided, assuming these to refer to the ODE and 2D PDE states')
    n = [n(1);0;0;n(2)];
elseif length(n)~=4
    error('Operator dimensions n must be specified as 4x2 or 4x1 array')
end
% Extract the size of the object: P\in R^(n0 x n0), Rxx\in L_2^(nx x nx)[x],
% Ryy\in L_2^(ny x ny)[y], R22\in L_2^(n2 x n2)[x,y],
n0=n(1);    nx=n(2);    ny=n(3);    n2=n(4);
if all(n==0)
    error('Error in poslpivar: All dimensions are zero')
end

% % Next, check the monomial degrees
if nargin>=3
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
    if length(d(:))==1
        d{2}=[d{1},d{1},2*d{1}];
        d{3}=d{2};
    else
        if length(d{2})==1
            d{2} = d{2}*ones(1,3);
        elseif length(d{2})==2
            d{2}(3) = max(d{2});
        end
        if numel(d)==2
            d{3}=d{2};
        else
            if length(d{3})==1
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
    if isa(d,'struct') && isfield(d,'Zx') && isfield(d,'Zy') && isfield(d,'Z2')
        % Monomials are specified directly rather than through degrees.
        use_monomials = 1;
        dx = d.Zx;
        dy = d.Zy;
        d2 = d.Z2;
    elseif isfield(d,'Zx') && (isfield(d,'dy') || isfield(d,'d2')) || ...
            isfield(d,'Zy') && (isfield(d,'dx') || isfield(d,'d2')) || ...
             isfield(d,'Z2') && (isfield(d,'dx') || isfield(d,'dy'))
        error('Specifying both monomials and maximal degrees is currently not supported')
    elseif (~isfield(d,'use_monomials') || ~d.use_monomials) ...
        && ~isfield(d,'dx') && ~isfield(d,'dy') && ~isfield(d,'d2')
        fprintf('\n Warning: No degrees are specified. Continuing with default values. \n')
    else
        if isfield(d,'use_monomials')
            use_monomials = d.use_monomials;
        end
        if ~isfield(d,'dx') || ~isfield(d,'dy') || ~isfield(d,'d2')
            error("Degrees for 2D positive operator should be specified as struct with fields 'dx', 'dy', and 'd2'. Call 'help posplivar_2d' for more details.")
        end
        dx = d.dx;
        dy = d.dy;
        d2 = d.d2;
    end
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
if ~use_monomials
if length(dx(:))==1
    dx{2}=[dx{1};dx{1};2*dx{1}];
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
% % Check if degrees for y monomials are properly specified
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
% % Check if degrees for 2D monomials are properly specified
if ~all(size(d2)==[3,3])
    if numel(d2)==1 && all(size(d2{1,1})>=[2,2])
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
    if n_updates_k>=5
        warning(['At least one of the (joint) degrees in d.d2{',num2str(k),'} is either too large or too small to make sense with all the joint degrees; reducing degrees to a sensible value'])
    end
end
end

% % To reduce complexity, allow certain terms to be excluded

% Check that list of parameters to exclude is properly specified.
excludeL = excludeL(:)';
if numel(excludeL)==4
    % Parameters to exclude are specified in 1D format
    % --> augment to 2D;
    excludeL = [excludeL(1),...
                excludeL(2:4),...
                excludeL(2:4),...
                reshape((excludeL').*excludeL,1,[])];
elseif numel(excludeL)~=16
    error('For 2D operator, options ''exclude'' should be specified as 1x16 array of binary values.')
end

% Check that option to enforce separable operator is properly specified.
if numel(is_sep)==1
    % Assume only full integral operators are used.
    is_sep = is_sep*ones(1,6);
elseif numel(is_sep)~=6
    error('For 2D operator, options ''sep'' should be specified as 1x6 array of binary values.')
end
% In separable case 1, set Ti3 = Ti4 and Zxa=Zxb (so that Rxx{2}=Rxx{3})
if is_sep(1)==1 
    excludeL(4) = 1;
end
% In separable case 2, set Ti6 = Ti7 and Zya=Zyb (so that Ryy{2}=Ryy{3})
if is_sep(2)==1 
    excludeL(7) = 1;
end
% In separable case 3, set Ti9 = Ti10 and Z2ao=Z2bo (so that R22{2,1}=R22{3,1})
if is_sep(3)==1 
    excludeL(10) = 1;
end
% In separable case 4, set Ti11 = Ti12 and Z2oa=Z2ob (so that R22{1,2}=R22{1,3})
if is_sep(4)==1 
    excludeL(12) = 1;
end
% In separable case 5, set Ti13 = Ti14 and Z2aa=Z2ba, and Ti15 = Ti14 and Z2ab=Z2bb  (so that R22{2,2}=R22{3,2} and R22{2,3}=R22{3,3})
% In separable case 6, set Ti13 = Ti15 and Z2aa=Z2ab, and Ti14 = Ti16 and Z2ba=Z2bb  (so that R22{2,2}=R22{2,3} and R22{3,2}=R22{3,3})
if is_sep(5)==1 && is_sep(6) 
    excludeL(14:16) = 1;
elseif is_sep(5)==1
    excludeL([14,16]) = 1;
elseif is_sep(6)==1
    excludeL([15,16]) = 1;
end

if n0==0
    excludeL(1)=1;
end
if nx==0
    excludeL(2:4)=[1 1 1];
end
if ny==0
    excludeL(5:7) = [1 1 1];
end
if n2==0
    excludeL(8:16) = [1 1 1 1 1 1 1 1 1];
end
    
if all(excludeL)
    error('You''re creating an empty dopvar! Please change options.exclude and/or options.sep and try again.')
end

%center_monomials = 0;
if use_monomials && center_monomials
    warning(['Centering monomials that are provided by the user is not recommended. '...
        'If an error is encountered, try running without centering, or try specifying the monomials using maximal degrees.'])
end


% % Define the variables and multiplier function

% Define the primary variables
if isempty(prog.vartable) || length(prog.vartable)<4
    error('Insufficient spatial variables have been specified in the LPI program structure.')
else
    var = polynomial(prog.vartable(1:4));
    var1 = var((1:2)');     var2 = var((3:4)');
    ss1 = var1(1);        ss2 = var1(2);
    tt1 = var2(1);        tt2 = var2(2);
end


% Define dummy variables for integration (replace eta)
rr1 = polynomial(1,1,{'rr1'},[1 1]);
rr2 = polynomial(1,1,{'rr2'},[1 1]);
rr = [rr1;rr2];

% Define the multiplier function to be used later
if psatz==0
    gss=polynomial(1);
elseif psatz==1
    gss=(ss1-I(1,1))*(I(1,2)-ss1)*(ss2-I(2,1))*(I(2,2)-ss2);
    
elseif psatz==2
    cntr = [mean(I(1,:)),mean(I(2,:))];
    rds = norm([I(1,2),I(2,2)] - cntr);
    
    gss = (rds^2 - (ss1-cntr(1))^2 - (ss2-cntr(2))^2);
    
else
    error('options.psatz can only assume values 0, 1, and 2')
end

% ONLY FOR TESTING PURPOSES: Take square of psatz polynomial
% if sign(options.psatz)==-1
%     gss = gss^2;  
% end

% Introduce variations with adjusted variables for later purposes
gtt = subs(gss,var1,var2);
grr = subs(gss,var1,rr);
gts = subs(gss,ss1,tt1);
gst = subs(gss,ss2,tt2);
grs = subs(gss,ss1,rr1);
gsr = subs(gss,ss2,rr2);
grt = subs(gtt,tt1,rr1);
gtr = subs(gtt,tt2,rr2);

% Introduce integrated versions for later purposes
gis = int(gss,ss1,I(1,1),I(1,2));
gsi = int(gss,ss2,I(2,1),I(2,2));
git = subs(gis,ss2,tt2);
gti = subs(gsi,ss1,tt1);
gir = subs(gis,ss2,rr2);
gri = subs(gsi,ss1,rr1);

gii = int(gis,ss2,I(2,1),I(2,2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Construct the monomial bases Z % % %

% % Build the monomials

if ~use_monomials
    % Construct the monomials based on the specified maximal degrees
    % Constructing Zxo(ss1)
    Zxo_s = build_monoms(ss1.varname,dx{1});
    
    % Constructing Zxa(ss1,tt1) and Zxb(ss1,tt1)
    % In this implementation, Zxa will have degree dx{2}(2) in tt1 and degree
    % dx{2}(1) in ss1 and max degree of ss1+tt1 is dx{2}(3). Similarly for Zxb(ss1,tt1)
    st_varname = [ss1.varname; tt1.varname];
    Zxa_st = build_monoms(st_varname,dx{2});
    Zxb_st = build_monoms(st_varname,dx{3});
    
    % % % % %
    
    % Constructing Zyo(ss2)
    Zyo_s = build_monoms(ss2.varname,dy{1});
    
    % Constructing Zya(ss2,tt2) and Zyb(ss2,tt2)
    % In this implementation, Zya will have degree dy{2}(2) in tt2 and degree
    % dy{2}(1) in ss2 and max degree of ss2+tt2 is dy{2}(3). Similarly for Zyb(ss2,tt2)
    st_varname = [ss2.varname; tt2.varname];
    Zya_st = build_monoms(st_varname,dy{2});
    Zyb_st = build_monoms(st_varname,dy{3});
    
    % % % % %
    
    % Constructing Z2oo(ss1,ss2)
    ss_varname = [ss1.varname; ss2.varname];
    Z2oo_ss = build_monoms(ss_varname,d2{1,1});
    
    % Constructing Z2ao(ss1,ss2,tt1) and Z2bo(ss1,ss2,tt1)
    sst_varname = [ss1.varname; tt1.varname; ss2.varname];
    Z2ao_sst = build_monoms(sst_varname,d2{2,1});
    Z2bo_sst = build_monoms(sst_varname,d2{3,1});
    
    % Constructing Z2oa(ss1,ss2,tt2) and Z2ob(ss1,ss2,tt2)
    sst_varname = [ss1.varname; ss2.varname; tt2.varname];
    Z2oa_sst = build_monoms(sst_varname,d2{1,2});
    Z2ob_sst = build_monoms(sst_varname,d2{1,3});
    
    % Constructing Z2aa(ss1,ss2,tt1,tt2) ... Z2bb(ss1,ss2,tt1,tt2)
    sstt_varname = [ss1.varname; tt1.varname; ss2.varname; tt2.varname];
    Z2aa_sstt = build_monoms(sstt_varname,d2{2,2});
    Z2ba_sstt = build_monoms(sstt_varname,d2{3,2});
    Z2ab_sstt = build_monoms(sstt_varname,d2{2,3});
    Z2bb_sstt = build_monoms(sstt_varname,d2{3,3});

else
    % In this case, dx{i}, dy{j} and d2{k} are all monomial degree tables.
    
    % Constructing Zxo(ss1), Zxa(ss1,tt1) and Zxb(ss1,tt1)
    st_varname = [ss1.varname; tt1.varname];
    if ~isa(dx{1},'polynomial')
        Zxo_s = polynomial(eye(size(dx{1},1)),dx{1},ss1.varname,[size(dx{1},1),1]);
    end
    if ~isa(dx{2},'polynomial')
        Zxa_st = polynomial(eye(size(dx{2},1)),dx{2},st_varname,[size(dx{2},1),1]);
    end
    if ~isa(dx{3},'polynomial')
        Zxb_st = polynomial(eye(size(dx{3},1)),dx{3},st_varname,[size(dx{3},1),1]);
    end
    
    % Constructing Zyo(ss2), Zya(ss2,tt2) and Zyb(ss2,tt2)
    st_varname = [ss2.varname; tt2.varname];
    if ~isa(dy{1},'polynomial')
        Zyo_s = polynomial(eye(size(dy{1},1)),dy{1},ss2.varname,[size(dy{1},1),1]);
    end
    if ~isa(dy{2},'polynomial')
        Zya_st = polynomial(eye(size(dy{2},1)),dy{2},st_varname,[size(dy{2},1),1]);
    end
    if ~isa(dy{3},'polynomial')
        Zyb_st = polynomial(eye(size(dy{3},1)),dy{3},st_varname,[size(dy{2},1),1]);
    end
    
    % Constructing Z2oo(ss1,ss2)
    ss_varname = [ss1.varname; ss2.varname];
    if ~isa(d2{1,1},'polynomial')
        Z2oo_ss = polynomial(eye(size(d2{1,1},1)),d2{1,1},ss_varname,[size(d2{1,1},1),1]);
    end

    % Constructing Z2ao(ss1,ss2,tt1) and Z2bo(ss1,ss2,tt1)
    sst_varname = [ss1.varname; tt1.varname; ss2.varname];
    if ~isa(d2{2,1},'polynomial')
        Z2ao_sst = polynomial(eye(size(d2{2,1},1)),d2{2,1},sst_varname,[size(d2{2,1},1),1]);
    end
    if ~isa(d2{3,1},'polynomial')
        Z2bo_sst = polynomial(eye(size(d2{3,1},1)),d2{3,1},sst_varname,[size(d2{3,1},1),1]);
    end

    % Constructing Z2oa(ss1,ss2,tt2) and Z2ob(ss1,ss2,tt2)
    sst_varname = [ss1.varname; ss2.varname; tt2.varname];
    if ~isa(d2{1,2},'polynomial')
        Z2oa_sst = polynomial(eye(size(d2{1,2},1)),d2{1,2},sst_varname,[size(d2{1,2},1),1]);
    end
    if ~isa(d2{1,3},'polynomial')
        Z2ob_sst = polynomial(eye(size(d2{1,3},1)),d2{1,3},sst_varname,[size(d2{1,3},1),1]);
    end

    % Constructing Z2aa(ss1,ss2,tt1,tt2) ... Z2bb(ss1,ss2,tt1,tt2)
    sstt_varname = [ss1.varname; tt1.varname; ss2.varname; tt2.varname];
    if ~isa(d2{2,2},'polynomial')
        Z2aa_sstt = polynomial(eye(size(d2{2,2},1)),d2{2,2},sstt_varname,[size(d2{2,2},1),1]);
    end
    if ~isa(d2{3,2},'polynomial')
        Z2ba_sstt = polynomial(eye(size(d2{3,2},1)),d2{3,2},sstt_varname,[size(d2{3,2},1),1]);
    end
    if ~isa(d2{2,3},'polynomial')
        Z2ab_sstt = polynomial(eye(size(d2{2,3},1)),d2{2,3},sstt_varname,[size(d2{2,3},1),1]);
    end
    if ~isa(d2{3,3},'polynomial')
        Z2bb_sstt = polynomial(eye(size(d2{3,3},1)),d2{3,3},sstt_varname,[size(d2{3,3},1),1]);
    end
    
end

if center_monomials
% Center all the monomials around the origin
cntr = mean(I,2);
domsz = I(:,2)-I(:,1);

ss1_new = (ss1-cntr(1))*2/domsz(1);     tt1_new = (tt1-cntr(1))*2/domsz(1);
ss2_new = (ss2-cntr(2))*2/domsz(2);     tt2_new = (tt2-cntr(2))*2/domsz(2);

% % For each of the monomials Z, establish coefficients Cf s.t. 
% % Cf2_bb*Z = Z_new, where the monomials Z_new are centered around 0

% Zxo(ss1), Zxa(ss1,tt1), Zxb(ss1,tt1)
Zxo_s = polynomial(subs(Zxo_s,ss1,ss1_new));
Zxa_st = polynomial(subs(Zxa_st,[ss1;tt1],[ss1_new;tt1_new]));
Zxb_st = polynomial(subs(Zxb_st,[ss1;tt1],[ss1_new;tt1_new]));

% Zyo(ss2), Zya(ss2,tt2), Zyb(ss2,tt2)
Zyo_s = polynomial(subs(Zyo_s,ss2,ss2_new));
Zya_st = polynomial(subs(Zya_st,[ss2;tt2],[ss2_new;tt2_new]));
Zyb_st = polynomial(subs(Zyb_st,[ss2;tt2],[ss2_new;tt2_new]));


% Zoo(ss1,ss2),     Zoa(ss1,ss2,tt2),     Zob(ss1,ss2,tt2),
% Zao(ss1,tt1,ss2), Zaa(ss1,tt1,ss2,tt2), Zab(ss1,tt1,ss2,tt2)
% Zbo(ss1,tt1,ss2), Zba(ss1,tt1,ss2,tt2), Zbb(ss1,tt1,ss2,tt2)
Z2oo_ss = polynomial(subs(Z2oo_ss,[ss1;ss2],[ss1_new;ss2_new]));
Z2ao_sst = polynomial(subs(Z2ao_sst,[ss1;tt1;ss2],[ss1_new;tt1_new;ss2_new]));
Z2bo_sst = polynomial(subs(Z2bo_sst,[ss1;tt1;ss2],[ss1_new;tt1_new;ss2_new]));
Z2oa_sst = polynomial(subs(Z2oa_sst,[ss1;ss2;tt2],[ss1_new;ss2_new;tt2_new]));
Z2ob_sst = polynomial(subs(Z2ob_sst,[ss1;ss2;tt2],[ss1_new;ss2_new;tt2_new]));
Z2aa_sstt = polynomial(subs(Z2aa_sstt,[ss1;tt1;ss2;tt2],[ss1_new;tt1_new;ss2_new;tt2_new]));
Z2ba_sstt = polynomial(subs(Z2ba_sstt,[ss1;tt1;ss2;tt2],[ss1_new;tt1_new;ss2_new;tt2_new]));
Z2ab_sstt = polynomial(subs(Z2ab_sstt,[ss1;tt1;ss2;tt2],[ss1_new;tt1_new;ss2_new;tt2_new]));
Z2bb_sstt = polynomial(subs(Z2bb_sstt,[ss1;tt1;ss2;tt2],[ss1_new;tt1_new;ss2_new;tt2_new]));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Construct the left and right monomial bases ZL and ZR % % %

% Introduce monomials with adjusted variables
Zxo_t = subs(Zxo_s,ss1,tt1);
Zxa_rt = subs(Zxa_st,ss1,rr1);
Zxa_rs = subs(Zxa_rt,tt1,ss1);
Zxb_rt = subs(Zxb_st,ss1,rr1);
Zxb_rs = subs(Zxb_rt,tt1,ss1);

Zyo_t = subs(Zyo_s,ss2,tt2);
Zya_rt = subs(Zya_st,ss2,rr2);
Zya_rs = subs(Zya_rt,tt2,ss2);
Zyb_rt = subs(Zyb_st,ss2,rr2);
Zyb_rs = subs(Zyb_rt,tt2,ss2);

Z2oo_tt = subs(Z2oo_ss,var1,var2);

Z2ao_rst = subs(Z2ao_sst,ss1,rr1);
Z2ao_rtt = subs(Z2ao_rst,ss2,tt2);
Z2ao_rss = subs(Z2ao_rst,tt1,ss1);
Z2bo_rst = subs(Z2bo_sst,ss1,rr1);
Z2bo_rtt = subs(Z2bo_rst,ss2,tt2);
Z2bo_rss = subs(Z2bo_rst,tt1,ss1);

Z2oa_srt = subs(Z2oa_sst,ss2,rr2);
Z2oa_trt = subs(Z2oa_srt,ss1,tt1);
Z2oa_srs = subs(Z2oa_srt,tt2,ss2);
Z2ob_srt = subs(Z2ob_sst,ss2,rr2);
Z2ob_trt = subs(Z2ob_srt,ss1,tt1);
Z2ob_srs = subs(Z2ob_srt,tt2,ss2);

Z2aa_rrtt = subs(Z2aa_sstt,var1,rr);
Z2aa_rrss = subs(Z2aa_rrtt,var2,var1);
Z2ba_rrtt = subs(Z2ba_sstt,var1,rr);
Z2ba_rrss = subs(Z2ba_rrtt,var2,var1);
Z2ab_rrtt = subs(Z2ab_sstt,var1,rr);
Z2ab_rrss = subs(Z2ab_rrtt,var2,var1);
Z2bb_rrtt = subs(Z2bb_sstt,var1,rr);
Z2bb_rrss = subs(Z2bb_rrtt,var2,var1);


% Collect the different monomials in cells ZL and ZR
includeL = ~excludeL;
ZL = cell(1,sum(includeL));
ZR = cell(1,sum(includeL));

mdim = [];
ndim = [];
Zdim = zeros(16,1);
indx = 1;
if includeL(1)
    ZL{indx} = 1;
    ZR{indx} = 1;
    mdim = [mdim;n0];
    ndim = [ndim;n0];
    indx = indx+1;
    Zdim(1) = n0;
end
if includeL(2)
    ZL{indx} = Zxo_s;
    ZR{indx} = Zxo_t;
    mdim = [mdim;nx];
    ndim = [ndim;nx];
    indx = indx+1;
    Zdim(2) = nx*size(Zxo_s,1);
end
if includeL(3)
    ZL{indx} = Zxa_rs;
    ZR{indx} = Zxa_rt;
    mdim = [mdim;nx];
    ndim = [ndim;nx];
    indx = indx+1;
    Zdim(3) = nx*size(Zxa_st,1);
end
if includeL(4)
    ZL{indx} = Zxb_rs;
    ZR{indx} = Zxb_rt;
    mdim = [mdim;nx];
    ndim = [ndim;nx];
    indx = indx+1;
    Zdim(4) = nx*size(Zxb_st,1);
end
if includeL(5)
    ZL{indx} = Zyo_s;
    ZR{indx} = Zyo_t;
    mdim = [mdim;ny];
    ndim = [ndim;ny];
    indx = indx+1;
    Zdim(5) = ny*size(Zyo_s,1);
end
if includeL(6)
    ZL{indx} = Zya_rs;
    ZR{indx} = Zya_rt;
    mdim = [mdim;ny];
    ndim = [ndim;ny];
    indx = indx+1;
    Zdim(6) = ny*size(Zya_st,1);
end
if includeL(7)
    ZL{indx} = Zyb_rs;
    ZR{indx} = Zyb_rt;
    mdim = [mdim;ny];
    ndim = [ndim;ny];
    indx = indx+1;
    Zdim(7) = ny*size(Zyb_st,1);
end
if includeL(8)
    ZL{indx} = Z2oo_ss;
    ZR{indx} = Z2oo_tt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
    Zdim(8) = n2*size(Z2oo_ss,1);
end
if includeL(9)
    ZL{indx} = Z2ao_rss;
    ZR{indx} = Z2ao_rtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
    Zdim(9) = n2*size(Z2ao_sst,1);
end
if includeL(10)
    ZL{indx} = Z2bo_rss;
    ZR{indx} = Z2bo_rtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
    Zdim(10) = n2*size(Z2bo_sst,1);
end
if includeL(11)
    ZL{indx} = Z2oa_srs;
    ZR{indx} = Z2oa_trt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
    Zdim(11) = n2*size(Z2oa_sst,1);
end
if includeL(12)
    ZL{indx} = Z2ob_srs;
    ZR{indx} = Z2ob_trt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
    Zdim(12) = n2*size(Z2ob_sst,1);
end
if includeL(13)
    ZL{indx} = Z2aa_rrss;
    ZR{indx} = Z2aa_rrtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
    Zdim(13) = n2*size(Z2aa_sstt,1);
end
if includeL(14)
    ZL{indx} = Z2ba_rrss;
    ZR{indx} = Z2ba_rrtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
    Zdim(14) = n2*size(Z2ba_sstt,1);
end
if includeL(15)
    ZL{indx} = Z2ab_rrss;
    ZR{indx} = Z2ab_rrtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
    Zdim(15) = n2*size(Z2ab_sstt,1);
end
if includeL(16)
    ZL{indx} = Z2bb_rrss;
    ZR{indx} = Z2bb_rrtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    Zdim(16) = n2*size(Z2bb_sstt,1);
    %indx = indx+1;
end

% Declare the monomial operator Zop.
if nargout>=4
    Zop = opvar2d([0,n0;0,nx;0,ny;sum(Zdim),n2],I,var1,var2);
    Zop.R20 = [eye(n0*includeL(1),n0); zeros(sum(Zdim(2:end)),n0)];
    Zop.R2x{1} = [zeros(Zdim(1),nx); kron(eye(nx*includeL(2)),Zxo_s); zeros(sum(Zdim(3:end)),nx)];
    Zop.R2x{2} = [zeros(sum(Zdim(1:2)),nx); kron(eye(nx*includeL(3)),Zxa_st); zeros(sum(Zdim(4:end)),nx)];
    if is_sep(1)
        Zop.R2x{3} = Zop.R2x{2};
    else
        Zop.R2x{3} = [zeros(sum(Zdim(1:3)),nx); kron(eye(nx*includeL(4)),Zxb_st); zeros(sum(Zdim(5:end)),nx)];
    end
    Zop.R2y{1} = [zeros(sum(Zdim(1:4)),ny); kron(eye(ny*includeL(5)),Zyo_s); zeros(sum(Zdim(6:end)),ny)];
    Zop.R2y{2} = [zeros(sum(Zdim(1:5)),ny); kron(eye(ny*includeL(6)),Zya_st); zeros(sum(Zdim(7:end)),ny)];
    if is_sep(2)
        Zop.R2y{3} = Zop.R2y{2};
    else
        Zop.R2y{3} = [zeros(sum(Zdim(1:6)),ny); kron(eye(ny*includeL(7)),Zyb_st); zeros(sum(Zdim(8:end)),ny)];
    end
    Zop.R22{1,1} = [zeros(sum(Zdim(1:7)),n2); kron(eye(n2*includeL(8)),Z2oo_ss); zeros(sum(Zdim(9:end)),n2)];
    Zop.R22{2,1} = [zeros(sum(Zdim(1:8)),n2); kron(eye(n2*includeL(9)),Z2ao_sst); zeros(sum(Zdim(10:end)),n2)];
    if is_sep(3)
        Zop.R22{3,1} = Zop.R22{2,1};
    else
        Zop.R22{3,1} = [zeros(sum(Zdim(1:9)),n2); kron(eye(n2*includeL(10)),Z2bo_sst); zeros(sum(Zdim(11:end)),n2)];
    end
    Zop.R22{1,2} = [zeros(sum(Zdim(1:10)),n2); kron(eye(n2*includeL(11)),Z2oa_sst); zeros(sum(Zdim(12:end)),n2)];
    if is_sep(4)
        Zop.R22{1,3} = Zop.R22{1,2};
    else
        Zop.R22{1,3} = [zeros(sum(Zdim(1:11)),n2); kron(eye(n2*includeL(12)),Z2ob_sst); zeros(sum(Zdim(13:end)),n2)];
    end
    Zop.R22{2,2} = [zeros(sum(Zdim(1:12)),n2); kron(eye(n2*includeL(13)),Z2aa_sstt); zeros(sum(Zdim(14:end)),n2)];
    if is_sep(5)
        Zop.R22{3,2} = Zop.R22{2,2};
    else
        Zop.R22{3,2} = [zeros(sum(Zdim(1:13)),n2); kron(eye(n2*includeL(14)),Z2ba_sstt); zeros(sum(Zdim(15:end)),n2)];
    end
    if is_sep(6)
        Zop.R22{2,3} = Zop.R22{2,2};
        Zop.R22{3,3} = Zop.R22{3,2};
    else
        Zop.R22{2,3} = [zeros(sum(Zdim(1:14)),n2); kron(eye(n2*includeL(15)),Z2ab_sstt); zeros(sum(Zdim(16:end)),n2)];
        if is_sep(5)
            Zop.R22{3,3} = Zop.R22{2,3};
        else
            Zop.R22{3,3} = [zeros(sum(Zdim(1:15)),n2); kron(eye(n2*includeL(15)),Z2bb_sstt)];
        end
    end
end


% Declare a decision matrix Qmat>=0, and compute the matrix-valued function
%   N = ZL'*Qmat*ZR
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


ij = cumsum(includeL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Build the positive decision operator % % %

Pop = dopvar2d([],[n0,n0; nx,nx; ny,ny; n2,n2],I,var1,var2);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Build R00, R0x, R0y and R02
if includeL(1)
    %Pop.R00 = Pop.R00 + gii*Q{1,1};
    Pop.R00 = Pop.R00 + gii*N{ij(1),ij(1)};
    if includeL(2)
        Pop.R0x = Pop.R0x + gsi*subs(N{ij(1),ij(2)},tt1,ss1);
    end
    if includeL(3) && ~is_sep(1)==1
        Pop.R0x = Pop.R0x + int_simple(gri*subs(N{ij(1),ij(3)},tt1,ss1),rr1,ss1,I(1,2));
    elseif includeL(3) && is_sep(1)==1
        Pop.R0x = Pop.R0x + int_simple(gri*subs(N{ij(1),ij(3)},tt1,ss1),rr1,I(1,1),I(1,2));
    end
    if includeL(4)
        Pop.R0x = Pop.R0x + int_simple(gri*subs(N{ij(1),ij(4)},tt1,ss1),rr1,I(1,1),ss1);
    end
    
    if includeL(5)
        Pop.R0y = Pop.R0y + gis*subs(N{ij(1),ij(5)},tt2,ss2);
    end
    if includeL(6) && ~is_sep(2)==1
        Pop.R0y = Pop.R0y + int_simple(gir*subs(N{ij(1),ij(6)},tt2,ss2),rr2,ss2,I(2,2));
    elseif includeL(6) && is_sep(2)==1
        Pop.R0y = Pop.R0y + int_simple(gir*subs(N{ij(1),ij(6)},tt2,ss2),rr2,I(2,1),I(2,2));
    end
    if includeL(7)
        Pop.R0y = Pop.R0y + int_simple(gir*subs(N{ij(1),ij(7)},tt2,ss2),rr2,I(2,1),ss2);
    end
    
    if includeL(8)
        Pop.R02 = Pop.R02 + gss*subs(N{ij(1),ij(8)},var2,var1);
    end
    if includeL(9) && ~is_sep(3)==1
        Pop.R02 = Pop.R02 + int_simple(grs*subs(N{ij(1),ij(9)},var2,var1),rr1,ss1,I(1,2));
    elseif includeL(9) && is_sep(3)==1
        Pop.R02 = Pop.R02 + int_simple(grs*subs(N{ij(1),ij(9)},var2,var1),rr1,I(1,1),I(1,2));
    end
    if includeL(10)
        Pop.R02 = Pop.R02 + int_simple(grs*subs(N{ij(1),ij(10)},var2,var1),rr1,I(1,1),ss1);
    end
    if includeL(11) && ~is_sep(4)==1
        Pop.R02 = Pop.R02 + int_simple(gsr*subs(N{ij(1),ij(11)},var2,var1),rr2,ss2,I(2,2));
    elseif includeL(11) && is_sep(4)==1
        Pop.R02 = Pop.R02 + int_simple(gsr*subs(N{ij(1),ij(11)},var2,var1),rr2,I(2,1),I(2,2));
    end
    if includeL(12)
        Pop.R02 = Pop.R02 + int_simple(gsr*subs(N{ij(1),ij(12)},var2,var1),rr2,I(2,1),ss2);
    end
    if includeL(13) && ~is_sep(5)==1 && ~is_sep(6)
        Pop.R02 = Pop.R02 + int_simple(int_simple(grr*subs(N{ij(1),ij(13)},var2,var1),rr1,ss1,I(1,2)),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(6)
        %
        Pop.R02 = Pop.R02 + int_simple(int_simple(grr*subs(N{ij(1),ij(13)},var2,var1),rr1,I(1,1),I(1,2)),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(5)
        %
        Pop.R02 = Pop.R02 + int_simple(int_simple(grr*subs(N{ij(1),ij(13)},var2,var1),rr1,ss1,I(1,2)),rr2,I(2,1),I(2,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        Pop.R02 = Pop.R02 + int_simple(int_simple(grr*subs(N{ij(1),ij(13)},var2,var1),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
    end
    if includeL(14) && ~is_sep(6)
        Pop.R02 = Pop.R02 + int_simple(int_simple(grr*subs(N{ij(1),ij(14)},var2,var1),rr1,I(1,1),ss1),rr2,ss2,I(2,2));
    elseif includeL(14) && is_sep(6)
        %
        Pop.R02 = Pop.R02 + int_simple(int_simple(grr*subs(N{ij(1),ij(14)},var2,var1),rr1,I(1,1),ss1),rr2,I(2,1),I(2,2));
    end
    if includeL(15) && ~is_sep(5)
        Pop.R02 = Pop.R02 + int_simple(int_simple(grr*subs(N{ij(1),ij(15)},var2,var1),rr1,ss1,I(1,2)),rr2,I(2,1),ss2);
    elseif includeL(15) && is_sep(5)
        %
        Pop.R02 = Pop.R02 + int_simple(int_simple(grr*subs(N{ij(1),ij(15)},var2,var1),rr1,I(1,1),I(1,2)),rr2,I(2,1),ss2);
    end
    if includeL(16) 
        Pop.R02 = Pop.R02 + int_simple(int_simple(grr*subs(N{ij(1),ij(16)},var2,var1),rr1,I(1,1),ss1),rr2,I(2,1),ss2);
    end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Build Rxx, Rxy and Rx2
if includeL(2)
    Pop.Rxx{1} = gsi * subs(N{ij(2),ij(2)},tt1,ss1); 
    if includeL(3) && ~is_sep(1)==1
        Pop.Rxx{2} = Pop.Rxx{2} + gsi * subs(N{ij(2),ij(3)},rr1,ss1); 
    elseif includeL(3) && is_sep(1)==1
        Pop.Rxx{2} = Pop.Rxx{2} + gsi * subs(N{ij(2),ij(3)},rr1,ss1) + gti * subs(N{ij(3),ij(2)},rr1,tt1); 
    end
    if includeL(4)
        Pop.Rxx{2} = Pop.Rxx{2} + gti * subs(N{ij(4),ij(2)},rr1,tt1); 
    end
    
    if includeL(5)
        Pop.Rxy = Pop.Rxy + gss * subs(N{ij(2),ij(5)},tt2,ss2);
    end
    if includeL(6) && ~is_sep(2)==1
        Pop.Rxy = Pop.Rxy + int_simple(gsr * subs(N{ij(2),ij(6)},tt2,ss2),rr2,ss2,I(2,2));
    elseif includeL(6) && is_sep(2)==1
        Pop.Rxy = Pop.Rxy + int_simple(gsr * subs(N{ij(2),ij(6)},tt2,ss2),rr2,I(2,1),I(2,2));
    end
    if includeL(7)
        Pop.Rxy = Pop.Rxy + int_simple(gsr * subs(N{ij(2),ij(7)},tt2,ss2),rr2,I(2,1),ss2);
    end
    
    if includeL(8)
        Pop.Rx2{1} = Pop.Rx2{1} + gss * subs(N{ij(2),ij(8)},var2,var1);
    end
    if includeL(11) && ~is_sep(4)==1
        Pop.Rx2{1} = Pop.Rx2{1} + int_simple(gsr * subs(N{ij(2),ij(11)},var2,var1),rr2,ss2,I(2,2));
    elseif includeL(11) && is_sep(4)==1
        Pop.Rx2{1} = Pop.Rx2{1} + int_simple(gsr * subs(N{ij(2),ij(11)},var2,var1),rr2,I(2,1),I(2,2));
    end
    if includeL(12)
        Pop.Rx2{1} = Pop.Rx2{1} + int_simple(gsr * subs(N{ij(2),ij(12)},var2,var1),rr2,I(2,1),ss2);
    end
    
    if includeL(9) && ~is_sep(3)
        Pop.Rx2{2} = Pop.Rx2{2} + gss * subs(N{ij(2),ij(9)},[rr1;tt2],var1);
    elseif includeL(9) && is_sep(3)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + gss * subs(N{ij(2),ij(9)},[rr1;tt2],var1);
        Pop.Rx2{3} = Pop.Rx2{3} + gss * subs(N{ij(2),ij(9)},[rr1;tt2],var1);
    end
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gsr * subs(N{ij(2),ij(13)},[rr1;tt2],var1),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(6)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gsr * subs(N{ij(2),ij(13)},[rr1;tt2],var1),rr2,ss2,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gsr * subs(N{ij(2),ij(13)},[rr1;tt2],var1),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(5)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gsr * subs(N{ij(2),ij(13)},[rr1;tt2],var1),rr2,I(2,1),I(2,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gsr * subs(N{ij(2),ij(13)},[rr1;tt2],var1),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gsr * subs(N{ij(2),ij(13)},[rr1;tt2],var1),rr2,I(2,1),I(2,2));
    end
    if includeL(15) && ~is_sep(6)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gsr * subs(N{ij(2),ij(15)},[rr1;tt2],var1),rr2,I(2,1),ss2);
    elseif includeL(15) && is_sep(6)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gsr * subs(N{ij(2),ij(15)},[rr1;tt2],var1),rr2,I(2,1),ss2);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gsr * subs(N{ij(2),ij(15)},[rr1;tt2],var1),rr2,I(2,1),ss2);
    end
    
    if includeL(10)
        Pop.Rx2{3} = Pop.Rx2{3} + gss * subs(N{ij(2),ij(10)},[rr1;tt2],var1);
    end
    if includeL(14) && ~is_sep(5)==1
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gsr * subs(N{ij(2),ij(14)},[rr1;tt2],var1),rr2,ss2,I(2,2));
    elseif includeL(14) && is_sep(5)==1
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gsr * subs(N{ij(2),ij(14)},[rr1;tt2],var1),rr2,I(2,1),I(2,2));
    end
    if includeL(16)
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gsr * subs(N{ij(2),ij(16)},[rr1;tt2],var1),rr2,I(2,1),ss2);
    end
end

if includeL(3) && ~is_sep(1)==1
    Pop.Rxx{2} = Pop.Rxx{2} + int_simple(gri * N{ij(3),ij(3)},rr1,ss1,I(1,2));
    if includeL(4)
        Pop.Rxx{2} = Pop.Rxx{2} + int_simple(gri * N{ij(4),ij(3)},rr1,tt1,ss1);
    end
    
    if includeL(5)
        Pop.Rxy = Pop.Rxy + int_simple(grs * subs(N{ij(3),ij(5)},tt2,ss2),rr1,ss1,I(1,2));
    end
    if includeL(6) && ~is_sep(2)==1
        Pop.Rxy = Pop.Rxy + int_simple(int_simple(grr * subs(N{ij(3),ij(6)},tt2,ss2),rr1,ss1,I(1,2)),rr2,ss2,I(2,2));
    elseif includeL(6) && is_sep(2)==1
        Pop.Rxy = Pop.Rxy + int_simple(int_simple(grr * subs(N{ij(3),ij(6)},tt2,ss2),rr1,ss1,I(1,2)),rr2,I(2,1),I(2,2));
    end
    if includeL(7)
        Pop.Rxy = Pop.Rxy + int_simple(int_simple(grr * subs(N{ij(3),ij(7)},tt2,ss2),rr1,ss1,I(1,2)),rr2,I(2,1),ss2);
    end
    
    if includeL(8)
        Pop.Rx2{3} = Pop.Rx2{3} + gts * subs(N{ij(3),ij(8)},[rr1;tt2],[tt1;ss2]);
    end
    if includeL(11) && ~is_sep(4)==1
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gtr * subs(N{ij(3),ij(11)},[rr1;tt2],[tt1;ss2]),rr2,ss2,I(2,2));
    elseif includeL(11) && is_sep(4)==1
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gtr * subs(N{ij(3),ij(11)},[rr1;tt2],[tt1;ss2]),rr2,I(2,1),I(2,2));
    end
    if includeL(12)
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gtr * subs(N{ij(3),ij(12)},[rr1;tt2],[tt1;ss2]),rr2,I(2,1),ss2);
    end
    
    if includeL(9) && ~is_sep(3)==1
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(grs * subs(N{ij(3),ij(9)},tt2,ss2),rr1,ss1,I(1,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(grs * subs(N{ij(3),ij(9)},tt2,ss2),rr1,tt1,I(1,2));
    elseif includeL(9) && is_sep(3)==1
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(grs * subs(N{ij(3),ij(9)},tt2,ss2),rr1,ss1,I(1,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(grs * subs(N{ij(3),ij(9)},tt2,ss2),rr1,ss1,I(1,2));
    end
    if includeL(13) && ~is_sep(5)==1 && ~is_sep(6)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,ss1,I(1,2)),rr2,ss2,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,tt1,I(1,2)),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(6)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,ss1,I(1,2)),rr2,ss2,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,ss1,I(1,2)),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(5)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,ss1,I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,tt1,I(1,2)),rr2,I(2,1),I(2,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,ss1,I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,ss1,I(1,2)),rr2,I(2,1),I(2,2));
    end
    if includeL(15) && ~is_sep(5)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(15)},tt2,ss2),rr1,ss1,I(1,2)),rr2,I(2,1),ss2);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(15)},tt2,ss2),rr1,tt1,I(1,2)),rr2,I(2,1),ss2);
    elseif includeL(15) && is_sep(5)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(15)},tt2,ss2),rr1,ss1,I(1,2)),rr2,I(2,1),ss2);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(15)},tt2,ss2),rr1,ss1,I(1,2)),rr2,I(2,1),ss2);
    end  
    
    if includeL(10)
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(grs * subs(N{ij(3),ij(10)},tt2,ss2),rr1,ss1,tt1);
    end
    if includeL(14) && ~is_sep(6)==1
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(14)},tt2,ss2),rr1,ss1,tt1),rr2,ss2,I(2,2));
    elseif includeL(14) && is_sep(6)==1
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(14)},tt2,ss2),rr1,ss1,tt1),rr2,I(2,1),I(2,2));
    end
    if includeL(16)
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(16)},tt2,ss2),rr1,ss1,tt1),rr2,I(2,1),ss2);
    end       
    
elseif includeL(3) && is_sep(1)==1
    Pop.Rxx{2} = Pop.Rxx{2} + int_simple(gri * N{ij(3),ij(3)},rr1,I(1,1),I(1,2));
    if includeL(5)
        Pop.Rxy = Pop.Rxy + int_simple(grs * subs(N{ij(3),ij(5)},tt2,ss2),rr1,I(1,1),I(1,2));
    end
    if includeL(6) && ~is_sep(2)==1
        Pop.Rxy = Pop.Rxy + int_simple(int_simple(grr * subs(N{ij(3),ij(6)},tt2,ss2),rr1,I(1,1),I(1,2)),rr2,ss2,I(2,2));
    elseif includeL(6) && is_sep(2)==1
        Pop.Rxy = Pop.Rxy + int_simple(int_simple(grr * subs(N{ij(3),ij(6)},tt2,ss2),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
    end
    if includeL(7)
        Pop.Rxy = Pop.Rxy + int_simple(int_simple(grr * subs(N{ij(3),ij(7)},tt2,ss2),rr1,I(1,1),I(1,2)),rr2,I(2,1),ss2);
    end
    
    if includeL(8)
        Pop.Rx2{3} = Pop.Rx2{3} + gts * subs(N{ij(3),ij(8)},[rr1;tt2],[tt1;ss2]);
        Pop.Rx2{2} = Pop.Rx2{2} + gts * subs(N{ij(3),ij(8)},[rr1;tt2],[tt1;ss2]);
    end
    if includeL(11) && ~is_sep(4)==1
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gtr * subs(N{ij(3),ij(11)},[rr1;tt2],[tt1;ss2]),rr2,ss2,I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gtr * subs(N{ij(3),ij(11)},[rr1;tt2],[tt1;ss2]),rr2,ss2,I(2,2));
    elseif includeL(11) && is_sep(4)==1
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gtr * subs(N{ij(3),ij(11)},[rr1;tt2],[tt1;ss2]),rr2,I(2,1),I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gtr * subs(N{ij(3),ij(11)},[rr1;tt2],[tt1;ss2]),rr2,I(2,1),I(2,2));
    end
    if includeL(12)
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(gtr * subs(N{ij(3),ij(12)},[rr1;tt2],[tt1;ss2]),rr2,I(2,1),ss2);
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gtr * subs(N{ij(3),ij(12)},[rr1;tt2],[tt1;ss2]),rr2,I(2,1),ss2);
    end
    
    if includeL(9) && ~is_sep(3)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(grs * subs(N{ij(3),ij(9)},tt2,ss2),rr1,tt1,I(1,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(grs * subs(N{ij(3),ij(9)},tt2,ss2),rr1,tt1,I(1,2));
    elseif includeL(9) && is_sep(3)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(grs * subs(N{ij(3),ij(9)},tt2,ss2),rr1,I(1,1),I(1,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(grs * subs(N{ij(3),ij(9)},tt2,ss2),rr1,I(1,1),I(1,2));
    end
    if includeL(13) && ~is_sep(5)==1 && ~is_sep(6)==1
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,tt1,I(1,2)),rr2,ss2,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,tt1,I(1,2)),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(6)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,I(1,1),I(1,2)),rr2,ss2,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,I(1,1),I(1,2)),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(5)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,tt1,I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,tt1,I(1,2)),rr2,I(2,1),I(2,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(13)},tt2,ss2),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));        
    end
    if includeL(15) && ~is_sep(5)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(15)},tt2,ss2),rr1,tt1,I(1,2)),rr2,I(2,1),ss2);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(15)},tt2,ss2),rr1,tt1,I(1,2)),rr2,I(2,1),ss2);
    elseif includeL(15) && is_sep(5)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(15)},tt2,ss2),rr1,I(1,1),I(1,2)),rr2,I(2,1),ss2);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(15)},tt2,ss2),rr1,I(1,1),I(1,2)),rr2,I(2,1),ss2);
    end  
    
    if includeL(10)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(grs * subs(N{ij(3),ij(10)},tt2,ss2),rr1,I(1,1),tt1);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(grs * subs(N{ij(3),ij(10)},tt2,ss2),rr1,I(1,1),tt1);
    end
    if includeL(14) && ~is_sep(6)==1
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(14)},tt2,ss2),rr1,I(1,1),tt1),rr2,ss2,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(14)},tt2,ss2),rr1,I(1,1),tt1),rr2,ss2,I(2,2));
    elseif includeL(14) && is_sep(6)==1
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(14)},tt2,ss2),rr1,I(1,1),tt1),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(14)},tt2,ss2),rr1,I(1,1),tt1),rr2,I(2,1),I(2,2));
    end
    if includeL(16)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(3),ij(16)},tt2,ss2),rr1,I(1,1),tt1),rr2,I(2,1),ss2);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(3),ij(16)},tt2,ss2),rr1,I(1,1),tt1),rr2,I(2,1),ss2);
    end  
    
end
if includeL(4)
    Pop.Rxx{2} = Pop.Rxx{2} + int_simple(gri * N{ij(4),ij(4)},rr1,I(1,1),tt1);
    if includeL(5)
        Pop.Rxy = Pop.Rxy + int_simple(grs * subs(N{ij(4),ij(5)},tt2,ss2),rr1,I(1,1),ss1);
    end
    if includeL(6) && ~is_sep(2)==1
        Pop.Rxy = Pop.Rxy + int_simple(int_simple(grr * subs(N{ij(4),ij(6)},tt2,ss2),rr1,I(1,1),ss1),rr2,ss2,I(2,2));
    elseif includeL(6) && is_sep(2)==1
        Pop.Rxy = Pop.Rxy + int_simple(int_simple(grr * subs(N{ij(4),ij(6)},tt2,ss2),rr1,I(1,1),ss1),rr2,I(2,1),I(2,2));
    end
    if includeL(7)
        Pop.Rxy = Pop.Rxy + int_simple(int_simple(grr * subs(N{ij(4),ij(7)},tt2,ss2),rr1,I(1,1),ss1),rr2,I(2,1),ss2);
    end
    
    if includeL(8)
        Pop.Rx2{2} = Pop.Rx2{2} + gts * subs(N{ij(4),ij(8)},[rr1;tt2],[tt1;ss2]);
    end
    if includeL(11) && ~is_sep(4)==1
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gtr * subs(N{ij(4),ij(11)},[rr1;tt2],[tt1;ss2]),rr2,ss2,I(2,2));
    elseif includeL(11) && is_sep(4)==1
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gtr * subs(N{ij(4),ij(11)},[rr1;tt2],[tt1;ss2]),rr2,I(2,1),I(2,2));
    end
    if includeL(12)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(gtr * subs(N{ij(4),ij(12)},[rr1;tt2],[tt1;ss2]),rr2,I(2,1),ss2);
    end
    
    if includeL(9) && ~is_sep(3)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(grs * subs(N{ij(4),ij(9)},tt2,ss2),rr1,tt1,ss1);
    elseif includeL(9) && is_sep(3)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(grs * subs(N{ij(4),ij(9)},tt2,ss2),rr1,I(1,1),ss1);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(grs * subs(N{ij(4),ij(9)},tt2,ss2),rr1,I(1,1),ss1);
    end
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(4),ij(13)},tt2,ss2),rr1,tt1,ss1),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(6)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(4),ij(13)},tt2,ss2),rr1,I(1,1),ss1),rr2,ss2,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(4),ij(13)},tt2,ss2),rr1,I(1,1),ss1),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(5)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(4),ij(13)},tt2,ss2),rr1,tt1,ss1),rr2,I(2,1),I(2,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(4),ij(13)},tt2,ss2),rr1,I(1,1),ss1),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(4),ij(13)},tt2,ss2),rr1,I(1,1),ss1),rr2,I(2,1),I(2,2));
    end
    if includeL(15) && ~is_sep(5)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(4),ij(15)},tt2,ss2),rr1,tt1,ss1),rr2,I(2,1),ss2);
    elseif includeL(15) && is_sep(5)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(4),ij(15)},tt2,ss2),rr1,I(1,1),ss1),rr2,I(2,1),ss2);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(4),ij(15)},tt2,ss2),rr1,I(1,1),ss1),rr2,I(2,1),ss2);
    end     
    
    if includeL(10)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(grs * subs(N{ij(4),ij(10)},tt2,ss2),rr1,I(1,1),tt1);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(grs * subs(N{ij(4),ij(10)},tt2,ss2),rr1,I(1,1),ss1);
    end
    if includeL(14) && ~is_sep(6)==1
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(4),ij(14)},tt2,ss2),rr1,I(1,1),tt1),rr2,ss2,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(4),ij(14)},tt2,ss2),rr1,I(1,1),ss1),rr2,ss2,I(2,2));
    elseif includeL(14) && is_sep(6)==1
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(4),ij(14)},tt2,ss2),rr1,I(1,1),tt1),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(4),ij(14)},tt2,ss2),rr1,I(1,1),ss1),rr2,I(2,1),I(2,2));
    end
    if includeL(16)
        Pop.Rx2{2} = Pop.Rx2{2} + int_simple(int_simple(grr * subs(N{ij(4),ij(16)},tt2,ss2),rr1,I(1,1),tt1),rr2,I(2,1),ss2);
        Pop.Rx2{3} = Pop.Rx2{3} + int_simple(int_simple(grr * subs(N{ij(4),ij(16)},tt2,ss2),rr1,I(1,1),ss1),rr2,I(2,1),ss2);
    end        
    
end

if is_sep(1)==1
    Pop.Rxx{3} = Pop.Rxx{2};
else
    Pop.Rxx{3} = var_swap(Pop.Rxx{2},ss1,tt1).';
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Build Ryy and Ry2
if includeL(5)
    Pop.Ryy{1} = gis * subs(N{ij(5),ij(5)},tt2,ss2); 
    if includeL(6) && ~is_sep(2)==1
        Pop.Ryy{2} = Pop.Ryy{2} + gis * subs(N{ij(5),ij(6)},rr2,ss2); 
    elseif includeL(6) && is_sep(1)==1
        Pop.Ryy{2} = Pop.Ryy{2} + gis * subs(N{ij(5),ij(6)},rr2,ss2) + git * subs(N{ij(6),ij(5)},rr2,tt2); 
    end
    if includeL(7)
        Pop.Ryy{2} = Pop.Ryy{2} + git * subs(N{ij(7),ij(5)},rr2,tt2); 
    end
    
    if includeL(8)
        Pop.Ry2{1} = Pop.Ry2{1} + gss * subs(N{ij(5),ij(8)},var2,var1); 
    end
    if includeL(9) && ~is_sep(3)==1
        Pop.Ry2{1} = Pop.Ry2{1} + int_simple(grs * subs(N{ij(5),ij(9)},var2,var1),rr1,ss1,I(1,2)); 
    elseif includeL(9) && is_sep(3)==1
        Pop.Ry2{1} = Pop.Ry2{1} + int_simple(grs * subs(N{ij(5),ij(9)},var2,var1),rr1,I(1,1),I(1,2)); 
    end
    if includeL(10)
        Pop.Ry2{1} = Pop.Ry2{1} + int_simple(grs * subs(N{ij(5),ij(10)},var2,var1),rr1,I(1,1),ss1); 
    end
    
    if includeL(11) && ~is_sep(4)
        %
        Pop.Ry2{2} = Pop.Ry2{2} + gss * subs(N{ij(5),ij(11)},[tt1;rr2],var1);
    elseif includeL(11) && is_sep(4)
        Pop.Ry2{2} = Pop.Ry2{2} + gss * subs(N{ij(5),ij(11)},[tt1;rr2],var1);
        Pop.Ry2{3} = Pop.Ry2{3} + gss * subs(N{ij(5),ij(11)},[tt1;rr2],var1);
    end
    if includeL(13) && ~is_sep(5)==1 && ~is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grs * subs(N{ij(5),ij(13)},[tt1;rr2],var1),rr1,ss1,I(1,2));
    elseif includeL(13) && ~is_sep(6)
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grs * subs(N{ij(5),ij(13)},[tt1;rr2],var1),rr1,I(1,1),I(1,2));
    elseif includeL(13) && ~is_sep(5)
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grs * subs(N{ij(5),ij(13)},[tt1;rr2],var1),rr1,ss1,I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grs * subs(N{ij(5),ij(13)},[tt1;rr2],var1),rr1,ss1,I(1,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grs * subs(N{ij(5),ij(13)},[tt1;rr2],var1),rr1,I(1,1),I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grs * subs(N{ij(5),ij(13)},[tt1;rr2],var1),rr1,I(1,1),I(1,2));
    end
    if includeL(14) && ~is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grs * subs(N{ij(5),ij(14)},[tt1;rr2],var1),rr1,I(1,1),ss1);
    elseif includeL(14) && is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grs * subs(N{ij(5),ij(14)},[tt1;rr2],var1),rr1,I(1,1),ss1);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grs * subs(N{ij(5),ij(14)},[tt1;rr2],var1),rr1,I(1,1),ss1);
    end
    
    if includeL(12)
        Pop.Ry2{3} = Pop.Ry2{3} + gss * subs(N{ij(5),ij(12)},[tt1;rr2],var1);
    end
    if includeL(15) && ~is_sep(5)==1
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grs * subs(N{ij(5),ij(15)},[tt1;rr2],var1),rr1,ss1,I(1,2));
    elseif includeL(15) && is_sep(5)==1
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grs * subs(N{ij(5),ij(15)},[tt1;rr2],var1),rr1,I(1,1),I(1,2));
    end
    if includeL(16)
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grs * subs(N{ij(5),ij(16)},[tt1;rr2],var1),rr1,I(1,1),ss1);
    end 
end

if includeL(6) && ~is_sep(2)==1
    Pop.Ryy{2} = Pop.Ryy{2} + int_simple(gir * N{ij(6),ij(6)},rr2,ss2,I(2,2));
    if includeL(7)
        Pop.Ryy{2} = Pop.Ryy{2} + int_simple(gir * N{ij(7),ij(6)},rr2,tt2,ss2);
    end
    
    if includeL(8)
        Pop.Ry2{3} = Pop.Ry2{3} + gst * subs(N{ij(6),ij(8)},[tt1;rr2],[ss1;tt2]);
    end
    if includeL(9) && ~is_sep(3)==1
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grt * subs(N{ij(6),ij(9)},[tt1;rr2],[ss1;tt2]),rr1,ss1,I(1,2));
    elseif includeL(9) && is_sep(3)==1
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grt * subs(N{ij(6),ij(9)},[tt1;rr2],[ss1;tt2]),rr1,I(1,1),I(1,2));
    end
    if includeL(10)
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grt * subs(N{ij(6),ij(10)},[tt1;rr2],[ss1;tt2]),rr1,I(1,1),ss1);
    end
    
    if includeL(11) && ~is_sep(4)==1
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(gsr * subs(N{ij(6),ij(11)},tt1,ss1),rr2,ss2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(gsr * subs(N{ij(6),ij(11)},tt1,ss1),rr2,tt2,I(2,2));
    elseif includeL(11) && is_sep(4)==1
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(gsr * subs(N{ij(6),ij(11)},tt1,ss1),rr2,ss2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(gsr * subs(N{ij(6),ij(11)},tt1,ss1),rr2,ss2,I(2,2));
    end   
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,ss2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,tt2,I(2,2));
    elseif includeL(13) && ~is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,ss2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,tt2,I(2,2));
    elseif includeL(13) && ~is_sep(5)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,ss2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,ss2,I(2,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,ss2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,ss2,I(2,2));
    end
    if includeL(14) && ~is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,ss2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,tt2,I(2,2));
    elseif includeL(14) && is_sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,ss2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,ss2,I(2,2));
    end
    
    if includeL(12)
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(gsr * subs(N{ij(6),ij(12)},tt1,ss1),rr2,ss2,tt2);
    end
    if includeL(15) && ~is_sep(5)==1
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(15)},tt1,ss1),rr1,ss1,I(1,2)),rr2,ss2,tt2);
    elseif includeL(15) && is_sep(5)==1
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(15)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,ss2,tt2);
    end
    if includeL(16)
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(16)},tt1,ss1),rr1,I(1,1),ss1),rr2,ss2,tt2);
    end       
   
elseif includeL(6) && is_sep(2)==1
    Pop.Ryy{2} = Pop.Ryy{2} + int_simple(gir * N{ij(6),ij(6)},rr2,I(2,1),I(2,2));
    
    if includeL(8)
        Pop.Ry2{2} = Pop.Ry2{2} + gst * subs(N{ij(6),ij(8)},[tt1;rr2],[ss1;tt2]);
        Pop.Ry2{3} = Pop.Ry2{3} + gst * subs(N{ij(6),ij(8)},[tt1;rr2],[ss1;tt2]);
    end
    if includeL(9) && ~is_sep(3)==1
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grt * subs(N{ij(6),ij(9)},[tt1;rr2],[ss1;tt2]),rr1,ss1,I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grt * subs(N{ij(6),ij(9)},[tt1;rr2],[ss1;tt2]),rr1,ss1,I(1,2));
    elseif includeL(9) && is_sep(3)==1
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grt * subs(N{ij(6),ij(9)},[tt1;rr2],[ss1;tt2]),rr1,I(1,1),I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grt * subs(N{ij(6),ij(9)},[tt1;rr2],[ss1;tt2]),rr1,I(1,1),I(1,2));
    end
    if includeL(10)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grt * subs(N{ij(6),ij(10)},[tt1;rr2],[ss1;tt2]),rr1,I(1,1),ss1);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(grt * subs(N{ij(6),ij(10)},[tt1;rr2],[ss1;tt2]),rr1,I(1,1),ss1);
    end
    
    if includeL(11) && ~is_sep(4)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(gsr * subs(N{ij(6),ij(11)},tt1,ss1),rr2,tt2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(gsr * subs(N{ij(6),ij(11)},tt1,ss1),rr2,tt2,I(2,2));
    elseif includeL(11) && is_sep(4)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(gsr * subs(N{ij(6),ij(11)},tt1,ss1),rr2,I(2,1),I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(gsr * subs(N{ij(6),ij(11)},tt1,ss1),rr2,I(2,1),I(2,2));
    end  
    
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,tt2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,tt2,I(2,2));
    elseif includeL(13) && ~is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,tt2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,tt2,I(2,2));
    elseif includeL(13) && ~is_sep(5)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,I(2,1),I(2,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
    end
    if includeL(14) && ~is_sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,tt2,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,tt2,I(2,2));
    elseif includeL(14) && is_sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,I(2,1),I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,I(2,1),I(2,2));
    end
    
    if includeL(12)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(gsr * subs(N{ij(6),ij(12)},tt1,ss1),rr2,I(2,1),tt2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(gsr * subs(N{ij(6),ij(12)},tt1,ss1),rr2,I(2,1),tt2);
    end
    if includeL(15) && ~is_sep(5)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(15)},tt1,ss1),rr1,ss1,I(1,2)),rr2,I(2,1),tt2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(15)},tt1,ss1),rr1,ss1,I(1,2)),rr2,I(2,1),tt2);
    elseif includeL(15) && is_sep(5)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(15)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,I(2,1),tt2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(15)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,I(2,1),tt2);
    end
    if includeL(16)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(6),ij(16)},tt1,ss1),rr1,I(1,1),ss1),rr2,I(2,1),tt2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(6),ij(16)},tt1,ss1),rr1,I(1,1),ss1),rr2,I(2,1),tt2);
    end   
end

if includeL(7)
    Pop.Ryy{2} = Pop.Ryy{2} + int_simple(gir * N{ij(7),ij(7)},rr2,I(2,1),tt2);
    
    if includeL(8)
        Pop.Ry2{2} = Pop.Ry2{2} + gst * subs(N{ij(7),ij(8)},[tt1;rr2],[ss1;tt2]);
    end
    if includeL(9) && ~is_sep(3)==1
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grt * subs(N{ij(7),ij(9)},[tt1;rr2],[ss1;tt2]),rr1,ss1,I(1,2));
    elseif includeL(9) && is_sep(3)==1
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grt * subs(N{ij(7),ij(9)},[tt1;rr2],[ss1;tt2]),rr1,I(1,1),I(1,2));
    end
    if includeL(10)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(grt * subs(N{ij(7),ij(10)},[tt1;rr2],[ss1;tt2]),rr1,I(1,1),ss1);
    end
    
    if includeL(11) && ~is_sep(4)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(gsr * subs(N{ij(7),ij(11)},tt1,ss1),rr2,tt2,ss2);
    elseif includeL(11) && is_sep(4)
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(gsr * subs(N{ij(7),ij(11)},tt1,ss1),rr2,I(1,1),ss2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(gsr * subs(N{ij(7),ij(11)},tt1,ss1),rr2,I(1,1),ss2);
    end    
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(7),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,tt2,ss2);
    elseif includeL(13) && ~is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(7),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,tt2,ss2);
    elseif includeL(13) && ~is_sep(5)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(7),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,I(2,1),ss2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(7),ij(13)},tt1,ss1),rr1,ss1,I(1,2)),rr2,I(2,1),ss2);
    elseif includeL(13) && is_sep(5) && is_sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(7),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,I(2,1),ss2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(7),ij(13)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,I(2,1),ss2);
    end
    if includeL(14) && ~is_sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(7),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,tt2,ss2);
    elseif includeL(14) && is_sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(7),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,I(2,1),ss2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(7),ij(14)},tt1,ss1),rr1,I(1,1),ss1),rr2,I(2,1),ss2);
    end
    
    if includeL(12)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(gsr * subs(N{ij(7),ij(12)},tt1,ss1),rr2,I(2,1),tt2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(gsr * subs(N{ij(7),ij(12)},tt1,ss1),rr2,I(2,1),ss2);
    end
    if includeL(15) && ~is_sep(5)==1
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(7),ij(15)},tt1,ss1),rr1,ss1,I(1,2)),rr2,I(2,1),tt2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(7),ij(15)},tt1,ss1),rr1,ss1,I(1,2)),rr2,I(2,1),ss2);
    elseif includeL(15) && is_sep(5)==1
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(7),ij(15)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,I(2,1),tt2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(7),ij(15)},tt1,ss1),rr1,I(1,1),I(1,2)),rr2,I(2,1),ss2);
    end
    if includeL(16)
        Pop.Ry2{2} = Pop.Ry2{2} + int_simple(int_simple(grr * subs(N{ij(7),ij(16)},tt1,ss1),rr1,I(1,1),ss1),rr2,I(2,1),tt2);
        Pop.Ry2{3} = Pop.Ry2{3} + int_simple(int_simple(grr * subs(N{ij(7),ij(16)},tt1,ss1),rr1,I(1,1),ss1),rr2,I(2,1),ss2);
    end        
end

if is_sep(2)==1
    Pop.Ryy{3} = Pop.Ryy{2};
else
    Pop.Ryy{3} = var_swap(Pop.Ryy{2},ss2,tt2).';
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Build R22
if includeL(8)
    Pop.R22{1,1} = Pop.R22{1,1} + gss * subs(N{ij(8),ij(8)},var2,var1);
    
    if includeL(9) && ~is_sep(3)
        Pop.R22{2,1} = Pop.R22{2,1} + gss * subs(N{ij(8),ij(9)},[rr1,tt2],var1);
    elseif includeL(9) && is_sep(3)
        Pop.R22{2,1} = Pop.R22{2,1} + gss * subs(N{ij(8),ij(9)},[rr1,tt2],var1) + gts * subs(N{ij(9),ij(8)},[rr1,tt2],[tt1;ss2]);
    end
    if includeL(10)
        Pop.R22{2,1} = Pop.R22{2,1} + gts * subs(N{ij(10),ij(8)},[rr1,tt2],[tt1;ss2]);
    end
    
    if includeL(11) && ~is_sep(4)
        Pop.R22{1,2} = Pop.R22{1,2} + gss * subs(N{ij(8),ij(11)},[tt1,rr2],var1);
    elseif includeL(11) && is_sep(4)
        Pop.R22{1,2} = Pop.R22{1,2} + gss * subs(N{ij(8),ij(11)},[tt1,rr2],var1) + gst * subs(N{ij(11),ij(8)},[tt1,rr2],[ss1;tt2]);
    end
    if includeL(12)
        Pop.R22{1,2} = Pop.R22{1,2} + gst * subs(N{ij(12),ij(8)},[tt1,rr2],[ss1;tt2]);
    end
    
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + gss * subs(N{ij(8),ij(13)},rr,var1);
    elseif includeL(13) && ~is_sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + gss * subs(N{ij(8),ij(13)},rr,var1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + gss * subs(N{ij(8),ij(13)},rr,var1);
    elseif includeL(13) && ~is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + gss * subs(N{ij(8),ij(13)},rr,var1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + gtt * subs(N{ij(13),ij(8)},rr,var2);
    elseif includeL(13) && is_sep(5) && is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + gss * subs(N{ij(8),ij(13)},rr,var1) + gtt * subs(N{ij(13),ij(8)},rr,var2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + gss * subs(N{ij(8),ij(13)},rr,var1) + gtt * subs(N{ij(13),ij(8)},rr,var2);
    end
    if includeL(14) && ~is_sep(6)
        Pop.R22{3,2} = Pop.R22{3,2} + gss * subs(N{ij(8),ij(14)},rr,var1);
    elseif includeL(14) && is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + gtt * subs(N{ij(14),ij(8)},rr,var2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + gss * subs(N{ij(8),ij(14)},rr,var1);
    end
    if includeL(15) && ~is_sep(5)
        Pop.R22{3,2} = Pop.R22{3,2} + gtt * subs(N{ij(15),ij(8)},rr,var2);
    elseif includeL(15) && is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + gtt * subs(N{ij(15),ij(8)},rr,var2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + gtt * subs(N{ij(15),ij(8)},rr,var2);
    end
    if includeL(16)
        Pop.R22{2,2} = Pop.R22{2,2} + gtt * subs(N{ij(16),ij(8)},rr,var2);
    end
    
end

if includeL(9) && ~is_sep(3)
    Pop.R22{2,1} = Pop.R22{2,1} + int_simple(grs * subs(N{ij(9),ij(9)},tt2,ss2),rr1,ss1,I(1,2));
    if includeL(10)
        Pop.R22{2,1} = Pop.R22{2,1} + int_simple(grs * subs(N{ij(10),ij(9)},tt2,ss2),rr1,tt1,ss1);
    end
    
    if includeL(11) && ~is_sep(4)
        Pop.R22{3,2} = Pop.R22{3,2} + gts * subs(N{ij(9),ij(11)},rr,[tt1;ss2]);
    elseif includeL(11) && is_sep(4)
        Pop.R22{2,2} = Pop.R22{2,2} + gst * subs(N{ij(11),ij(9)},rr,[ss1;tt2]);
        Pop.R22{3,2} = Pop.R22{3,2} + gts * subs(N{ij(9),ij(11)},rr,[tt1;ss2]);
    end
    if includeL(12)
        Pop.R22{2,2} = Pop.R22{2,2} + gst * subs(N{ij(12),ij(9)},rr,[ss1;tt2]);
    end
    
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(9),ij(13)},rr2,ss2),rr1,ss1,I(1,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(grs * subs(N{ij(9),ij(13)},rr2,ss2),rr1,tt1,I(1,2));
    elseif includeL(13) && ~is_sep(6)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(9),ij(13)},rr2,ss2),rr1,ss1,I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var11,I(1,2));
    elseif includeL(13) && ~is_sep(5)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(9),ij(13)},rr2,ss2),rr1,ss1,I(1,2))...
                                    + int_simple(grt * subs(N{ij(13),ij(9)},rr2,tt2),rr1,ss1,I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var12,I(1,2))...
        % %                             + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,var12,I(1,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(9),ij(13)},rr2,ss2),rr1,ss1,I(1,2))...
                                    + int_simple(grt * subs(N{ij(13),ij(9)},rr2,tt2),rr1,tt1,I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var11,I(1,2))...
        % %                             + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,var12,I(1,2));
    end
    if includeL(14) && ~is_sep(6)
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(grs * subs(N{ij(9),ij(14)},rr2,ss2),rr1,ss1,tt1);
    elseif includeL(14) && is_sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grt * subs(N{ij(14),ij(9)},rr2,tt2),rr1,tt1,ss1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(14)},rr2,var21),rr1,var11,var12);
    end
    if includeL(15) && ~is_sep(5)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grt * subs(N{ij(15),ij(9)},rr2,tt2),rr1,ss1,I(1,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(grt * subs(N{ij(15),ij(9)},rr2,tt2),rr1,tt1,I(1,2));
    elseif includeL(15) && is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grt * subs(N{ij(15),ij(9)},rr2,tt2),rr1,tt1,I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(15),ij(9)},rr2,var22),rr1,var12,I(1,2));
    end
    if includeL(16)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grt * subs(N{ij(16),ij(9)},rr2,tt2),rr1,tt1,ss1);
    end
    
elseif includeL(9) && is_sep(3)
    Pop.R22{2,1} = Pop.R22{2,1} + int_simple(grs * subs(N{ij(9),ij(9)},tt2,ss2),rr1,I(1,1),I(1,2));
    
    if includeL(11) && ~is_sep(4)
        Pop.R22{2,2} = Pop.R22{2,2} + gts * subs(N{ij(9),ij(11)},rr,[tt1;ss2]);
        Pop.R22{3,2} = Pop.R22{3,2} + gts * subs(N{ij(9),ij(11)},rr,[tt1;ss2]);
    elseif includeL(11) && is_sep(4)
        Pop.R22{2,2} = Pop.R22{2,2} + gst * subs(N{ij(11),ij(9)},rr,[ss1;tt2]) + gts * subs(N{ij(9),ij(11)},rr,[tt1;ss2]);
        Pop.R22{3,2} = Pop.R22{3,2} + gst * subs(N{ij(11),ij(9)},rr,[ss1;tt2]) + gts * subs(N{ij(9),ij(11)},rr,[tt1;ss2]);
    end
    if includeL(12)
        Pop.R22{2,2} = Pop.R22{2,2} + gst * subs(N{ij(12),ij(9)},rr,[ss1;tt2]);
        Pop.R22{3,2} = Pop.R22{3,2} + gst * subs(N{ij(12),ij(9)},rr,[ss1;tt2]);
    end
    
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(9),ij(13)},rr2,ss2),rr1,tt1,I(1,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(grs * subs(N{ij(9),ij(13)},rr2,ss2),rr1,tt1,I(1,2));
    elseif includeL(13) && ~is_sep(6)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(9),ij(13)},rr2,ss2),rr1,I(1,1),I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,I(1,1),I(1,2));
    elseif includeL(13) && ~is_sep(5)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(9),ij(13)},rr2,ss2),rr1,tt1,I(1,2))...
                                    + int_simple(grt * subs(N{ij(13),ij(9)},rr2,tt2),rr1,ss1,I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var12,I(1,2))...
        % %                             + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,var11,I(1,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(9),ij(13)},rr2,ss2),rr1,I(1,1),I(1,2))...
                                    + int_simple(grt * subs(N{ij(13),ij(9)},rr2,tt2),rr1,I(1,1),I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,I(1,1),I(1,2))...
        % %                             + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,I(1,1),I(1,2));
    end
    if includeL(14) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(9),ij(14)},rr2,ss2),rr1,I(1,1),tt1);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(grs * subs(N{ij(9),ij(14)},rr2,ss2),rr1,I(1,1),tt1);
    elseif includeL(14) && is_sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(9),ij(14)},rr2,ss2),rr1,I(1,1),tt1)...
                                    + int_simple(grt * subs(N{ij(14),ij(9)},rr2,tt2),rr1,I(1,1),ss1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(14)},rr2,var21),rr1,I(1,1),var12)...
        % %                             + int(grt * subs(N{ij(14),ij(9)},rr2,var22),rr1,I(1,1),var11);
    end
    if includeL(15) && ~is_sep(5)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grt * subs(N{ij(15),ij(9)},rr2,tt2),rr1,ss1,I(1,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(grt * subs(N{ij(15),ij(9)},rr2,tt2),rr1,ss1,I(1,2));
    elseif includeL(15) && is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grt * subs(N{ij(15),ij(9)},rr2,tt2),rr1,I(1,1),I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(15),ij(9)},rr2,var22),rr1,I(1,1),I(1,2));
    end
    if includeL(16)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grt * subs(N{ij(16),ij(9)},rr2,tt2),rr1,I(1,1),ss1);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(grt * subs(N{ij(16),ij(9)},rr2,tt2),rr1,I(1,1),ss1);
    end
    
end

if includeL(10)
    Pop.R22{2,1} = Pop.R22{2,1} + int_simple(grs * subs(N{ij(10),ij(10)},tt2,ss2),rr1,I(1,1),tt1);
    
    if includeL(11) && ~is_sep(4)
        Pop.R22{2,2} = Pop.R22{2,2} + gts * subs(N{ij(10),ij(11)},rr,[tt1;ss2]);
    elseif includeL(11) && is_sep(4)
        Pop.R22{2,2} = Pop.R22{2,2} + gts * subs(N{ij(10),ij(11)},rr,[tt1;ss2]);
        Pop.R22{3,2} = Pop.R22{3,2} + gst * subs(N{ij(11),ij(10)},rr,[ss1;tt2]);
    end
    if includeL(12)
        Pop.R22{3,2} = Pop.R22{3,2} + gst * subs(N{ij(12),ij(10)},rr,[ss1;tt2]);
    end
    
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(10),ij(13)},rr2,ss2),rr1,tt1,ss1);
    elseif includeL(13) && ~is_sep(6)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(10),ij(13)},rr2,ss2),rr1,I(1,1),ss1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(10),ij(13)},rr2,var21),rr1,I(1,1),var11);
    elseif includeL(13) && ~is_sep(5)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(10),ij(13)},rr2,ss2),rr1,tt1,ss1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(13),ij(10)},rr2,var22),rr1,var11,var12);
    elseif includeL(13) && is_sep(5) && is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(10),ij(13)},rr2,ss2),rr1,I(1,1),ss1)...
                                    + int_simple(grt * subs(N{ij(13),ij(10)},rr2,tt2),rr1,I(1,1),tt1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(10),ij(13)},rr2,var21),rr1,I(1,1),var11)...
        % %                             + int(grt * subs(N{ij(13),ij(10)},rr2,var22),rr1,I(1,1),var12);
    end
    if includeL(14) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(10),ij(14)},rr2,ss2),rr1,I(1,1),tt1);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(grs * subs(N{ij(10),ij(14)},rr2,ss2),rr1,I(1,1),ss1);
    elseif includeL(14) && is_sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grs * subs(N{ij(10),ij(14)},rr2,ss2),rr1,I(1,1),tt1)...
                                    + int_simple(grt * subs(N{ij(14),ij(10)},rr2,tt2),rr1,I(1,1),tt1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(10),ij(14)},rr2,var21),rr1,I(1,1),var11)...
        % %                             + int(grt * subs(N{ij(14),ij(10)},rr2,var22),rr1,I(1,1),var11);
    end
    if includeL(15) && ~is_sep(5)
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(grt * subs(N{ij(15),ij(10)},rr2,tt2),rr1,ss1,tt1);
    elseif includeL(15) && is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grt * subs(N{ij(15),ij(10)},rr2,tt2),rr1,I(1,1),tt1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(15),ij(10)},rr2,var22),rr1,I(1,1),var12);
    end
    if includeL(16)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(grt * subs(N{ij(16),ij(10)},rr2,tt2),rr1,I(1,1),tt1);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(grt * subs(N{ij(16),ij(10)},rr2,tt2),rr1,I(1,1),ss1);
    end
    
end

if includeL(11) && ~is_sep(4)
    Pop.R22{1,2} = Pop.R22{1,2} + int_simple(gsr * subs(N{ij(11),ij(11)},tt1,ss1),rr2,ss2,I(2,2));
    if includeL(12)
        Pop.R22{1,2} = Pop.R22{1,2} + int_simple(gsr * subs(N{ij(12),ij(11)},tt1,ss1),rr2,tt2,ss2);
    end
    
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(11),ij(13)},rr1,ss1),rr2,ss2,I(2,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(gtr * subs(N{ij(13),ij(11)},rr1,tt1),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(11),ij(13)},rr1,ss1),rr2,ss2,I(2,2))...
                                    + int_simple(gtr * subs(N{ij(13),ij(11)},rr1,tt1),rr2,ss2,I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var21,I(2,2))...
        % %                             + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var21,I(2,2));
    elseif includeL(13) && ~is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(11),ij(13)},rr1,ss1),rr2,ss2,I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var22,I(2,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(11),ij(13)},rr1,ss1),rr2,ss2,I(2,2))...
                                    + int_simple(gtr * subs(N{ij(13),ij(11)},rr1,tt1),rr2,tt2,I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var22,I(2,2))...
        % %                             + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var21,I(2,2));
    end
    if includeL(14) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gtr * subs(N{ij(14),ij(11)},rr1,tt1),rr2,ss2,I(2,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(gsr * subs(N{ij(11),ij(14)},rr1,ss1),rr2,ss2,I(2,2));
    elseif includeL(14) && is_sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gtr * subs(N{ij(14),ij(11)},rr1,tt1),rr2,tt2,I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(11),ij(14)},rr1,var11),rr2,var21,I(2,2));
    end
    if includeL(15) && ~is_sep(5)
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(gtr * subs(N{ij(15),ij(11)},rr1,tt1),rr2,tt2,ss2);
    elseif includeL(15) && is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gtr * subs(N{ij(15),ij(11)},rr1,tt1),rr2,tt2,ss2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(15),ij(11)},rr1,var12),rr2,var22,var21);
    end
    if includeL(16)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gtr * subs(N{ij(16),ij(11)},rr1,tt1),rr2,tt2,ss2);
    end
    
elseif includeL(11) && is_sep(4)
    Pop.R22{1,2} = Pop.R22{1,2} + int_simple(gsr * subs(N{ij(11),ij(11)},tt1,ss1),rr2,I(2,1),I(2,2));
    
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(11),ij(13)},rr1,ss1),rr2,tt2,I(2,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(gtr * subs(N{ij(13),ij(11)},rr1,tt1),rr2,ss2,I(2,2));
    elseif includeL(13) && ~is_sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(11),ij(13)},rr1,ss1),rr2,tt2,I(2,2))...
                                    + int_simple(gtr * subs(N{ij(13),ij(11)},rr1,tt1),rr2,ss2,I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var21,I(2,2))...
        % %                             + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var22,I(2,2));
    elseif includeL(13) && ~is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(11),ij(13)},rr1,ss1),rr2,I(2,1),I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,I(2,1),I(2,2));
    elseif includeL(13) && is_sep(5) && is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(11),ij(13)},rr1,ss1),rr2,I(2,1),I(2,2))...
                                    + int_simple(gtr * subs(N{ij(13),ij(11)},rr1,tt1),rr2,I(2,1),I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,I(2,1),I(2,2))...
        % %                             + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,I(2,1),I(2,2));
    end    
    if includeL(14) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gtr * subs(N{ij(14),ij(11)},rr1,tt1),rr2,ss2,I(2,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(gsr * subs(N{ij(11),ij(14)},rr1,ss1),rr2,tt2,I(2,2));
    elseif includeL(14) && is_sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gtr * subs(N{ij(14),ij(11)},rr1,tt1),rr2,I(2,1),I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(11),ij(14)},rr1,var11),rr2,I(2,1),I(2,2));
    end
    if includeL(15) && ~is_sep(5)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(11),ij(15)},rr1,ss1),rr2,I(2,1),tt2);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(gtr * subs(N{ij(15),ij(11)},rr1,tt1),rr2,I(2,1),ss2);
    elseif includeL(15) && is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(11),ij(15)},rr1,ss1),rr2,I(2,1),tt2)...
                                    + int_simple(gtr * subs(N{ij(15),ij(11)},rr1,tt1),rr2,I(2,1),ss2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(15),ij(11)},rr1,var12),rr2,I(2,1),var21)...
        % %                             + int(gsr * subs(N{ij(11),ij(15)},rr1,var11),rr2,I(2,1),var22);
    end
    if includeL(16)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gtr * subs(N{ij(16),ij(11)},rr1,tt1),rr2,I(2,1),ss2);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(gsr * subs(N{ij(11),ij(16)},rr1,ss1),rr2,I(2,1),tt2);
    end
    
end

if includeL(12)
    Pop.R22{1,2} = Pop.R22{1,2} + int_simple(gsr * subs(N{ij(12),ij(12)},tt1,ss1),rr2,I(2,1),tt2);
    
    if includeL(13) && ~is_sep(5) && ~is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(12),ij(13)},rr1,ss1),rr2,tt2,ss2);
    elseif includeL(13) && ~is_sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(12),ij(13)},rr1,ss1),rr2,tt2,ss2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(12),ij(13)},rr1,var11),rr2,var22,var21);
    elseif includeL(13) && ~is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(12),ij(13)},rr1,ss1),rr2,I(2,1),ss2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(12)},rr1,var12),rr2,I(2,1),var22);
    elseif includeL(13) && is_sep(5) && is_sep(6)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(12),ij(13)},rr1,ss1),rr2,I(2,1),ss2)...
                                    + int_simple(gtr * subs(N{ij(13),ij(12)},rr1,tt1),rr2,I(2,1),tt2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(12)},rr1,var12),rr2,I(2,1),var22)...
        % %                             + int(gsr * subs(N{ij(12),ij(13)},rr1,var11),rr2,I(2,1),var21);
    end 
    
    if includeL(14) && ~is_sep(6)
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(gsr * subs(N{ij(12),ij(14)},rr1,ss1),rr2,tt2,ss2);
    elseif includeL(14) && is_sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gtr * subs(N{ij(14),ij(12)},rr1,tt1),rr2,I(2,1),tt2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(12),ij(14)},rr1,var11),rr2,I(2,1),var21);
    end
    if includeL(15) && ~is_sep(5)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(12),ij(15)},rr1,ss1),rr2,I(2,1),tt2);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(gtr * subs(N{ij(15),ij(12)},rr1,tt1),rr2,I(2,1),tt2);
    elseif includeL(15) && is_sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gsr * subs(N{ij(12),ij(15)},rr1,ss1),rr2,I(2,1),tt2)...
                                    + int_simple(gtr * subs(N{ij(15),ij(12)},rr1,tt1),rr2,I(2,1),tt2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(15),ij(12)},rr1,var12),rr2,I(2,1),var22)...
        % %                             + int(gsr * subs(N{ij(12),ij(15)},rr1,var11),rr2,I(2,1),var22);
    end
    if includeL(16)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(gtr * subs(N{ij(16),ij(12)},rr1,tt1),rr2,I(2,1),tt2);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(gsr * subs(N{ij(12),ij(16)},rr1,ss1),rr2,I(2,1),tt2);
    end    
end

if includeL(13) && ~is_sep(5) && ~is_sep(6)
    Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(13),ij(13)},rr1,ss1,I(1,2)),rr2,ss2,I(2,2));
    Pop.R22{3,2} = Pop.R22{3,2} + int_simple(int_simple(grr * N{ij(13),ij(13)},rr1,tt1,I(1,2)),rr2,ss2,I(2,2));
    
    if includeL(14)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(14),ij(13)},rr1,tt1,ss1),rr2,ss2,I(2,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(int_simple(grr * N{ij(13),ij(14)},rr1,ss1,tt1),rr2,ss2,I(2,2));
    end
    if includeL(15)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(15),ij(13)},rr1,ss1,I(1,2)),rr2,tt2,ss2);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(int_simple(grr * N{ij(15),ij(13)},rr1,tt1,I(1,2)),rr2,tt2,ss2);
    end
    if includeL(16)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(16),ij(13)},rr1,tt1,ss1),rr2,tt2,ss2);
    end
    
elseif includeL(13) && ~is_sep(6)
    %
    Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(13),ij(13)},rr1,I(1,1),I(1,2)),rr2,ss2,I(2,2));
    % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(13),ij(13)},rr1,I(1,1),I(1,2)),rr2,var21,I(2,2));
    
    if includeL(15)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(15),ij(13)},rr1,I(1,1),I(1,2)),rr2,tt2,ss2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(15),ij(13)},rr1,I(1,1),I(1,2)),rr2,var22,var21);
    end
    
elseif includeL(13) && ~is_sep(5)
    %
    Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(13),ij(13)},rr1,ss1,I(1,2)),rr2,I(2,1),I(2,2));
    % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(13),ij(13)},rr1,var12,I(1,2)),rr2,I(2,1),I(2,2));
    
    if includeL(14)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(14),ij(13)},rr1,tt1,ss1),rr2,I(2,1),I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(13),ij(14)},rr1,var11,var12),rr2,I(2,1),I(2,2));
    end

elseif includeL(13) && is_sep(5) && is_sep(6)
    Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(13),ij(13)},rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
    % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(13),ij(13)},rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
end

if includeL(14) && ~is_sep(6)
    Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(14),ij(14)},rr1,I(1,1),tt1),rr2,ss2,I(2,2));
    Pop.R22{3,2} = Pop.R22{3,2} + int_simple(int_simple(grr * N{ij(14),ij(14)},rr1,I(1,1),ss1),rr2,ss2,I(2,2));
    
    if includeL(15)
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(int_simple(grr * N{ij(15),ij(14)},rr1,ss1,tt1),rr2,tt2,ss2);
    end
    if includeL(16)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(16),ij(14)},rr1,I(1,1),tt1),rr2,tt2,ss2);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(int_simple(grr * N{ij(16),ij(14)},rr1,I(1,1),ss1),rr2,tt2,ss2);
    end
    
elseif includeL(14) && is_sep(6)
    %
    Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(14),ij(14)},rr1,I(1,1),tt1),rr2,I(2,1),I(2,2));
    % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(14),ij(14)},rr1,I(1,1),var11),rr2,I(2,1),I(2,2));
    
    if includeL(15)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(14),ij(15)},rr1,tt1,ss1),rr2,I(2,1),tt2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(15),ij(14)},rr1,var11,var12),rr2,I(2,1),var21);
    end
    
end

if includeL(15) && ~is_sep(5)
    Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(15),ij(15)},rr1,ss1,I(1,2)),rr2,I(2,1),tt2);
    Pop.R22{3,2} = Pop.R22{3,2} + int_simple(int_simple(grr * N{ij(15),ij(15)},rr1,tt1,I(1,2)),rr2,I(2,1),tt2);
    
    if includeL(16)
        Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(16),ij(15)},rr1,tt1,ss1),rr2,I(2,1),tt2);
        Pop.R22{3,2} = Pop.R22{3,2} + int_simple(int_simple(grr * N{ij(15),ij(16)},rr1,ss1,tt1),rr2,I(2,1),tt2);
    end
    
elseif includeL(15) && is_sep(5)
    %
    Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(15),ij(15)},rr1,I(1,1),I(1,2)),rr2,I(2,1),tt2);
    % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(15),ij(15)},rr1,I(1,1),I(1,2)),rr2,I(2,1),var22);
    
end

if includeL(16)
    Pop.R22{2,2} = Pop.R22{2,2} + int_simple(int_simple(grr * N{ij(16),ij(16)},rr1,I(1,1),tt1),rr2,I(2,1),tt2);
    Pop.R22{3,2} = Pop.R22{3,2} + int_simple(int_simple(grr * N{ij(16),ij(16)},rr1,I(1,1),ss1),rr2,I(2,1),tt2);
end

if is_sep(3)==1
    Pop.R22{3,1} = Pop.R22{2,1};
else
    Pop.R22{3,1} = var_swap(Pop.R22{2,1},ss1,tt1).';
end
if is_sep(4)==1
    Pop.R22{1,3} = Pop.R22{1,2};
else
    Pop.R22{1,3} = var_swap(Pop.R22{1,2},ss2,tt2).';
end
if is_sep(5) && is_sep(6)
    Pop.R22{3,2} = Pop.R22{2,2};
    Pop.R22{2,3} = Pop.R22{2,2};
    Pop.R22{3,3} = Pop.R22{2,2};
elseif is_sep(5)
    Pop.R22{3,2} = Pop.R22{2,2};
    Pop.R22{2,3} = var_swap(var_swap(Pop.R22{3,2},ss1,tt1),ss2,tt2).';
    Pop.R22{3,3} = Pop.R22{2,3};
elseif is_sep(6)
    Pop.R22{2,3} = Pop.R22{2,2};
    Pop.R22{3,2} = var_swap(var_swap(Pop.R22{2,3},ss1,tt1),ss2,tt2).';
    Pop.R22{3,3} = Pop.R22{3,2};
else
    Pop.R22{3,3} = var_swap(var_swap(Pop.R22{2,2},ss1,tt1),ss2,tt2).';
    Pop.R22{2,3} = var_swap(var_swap(Pop.R22{3,2},ss1,tt1),ss2,tt2).';
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

Pop.Rx0 = Pop.R0x.';
Pop.Ry0 = Pop.R0y.';
Pop.R20 = Pop.R02.';

Pop.Ryx = Pop.Rxy.';
Pop.R2x{1,1} = Pop.Rx2{1,1}.';
Pop.R2x{2,1} = var_swap(Pop.Rx2{3,1}.',ss1,tt1);
Pop.R2x{3,1} = var_swap(Pop.Rx2{2,1}.',ss1,tt1);
Pop.R2y{1,1} = Pop.Ry2{1,1}.';
Pop.R2y{1,2} = var_swap(Pop.Ry2{1,3}.',ss2,tt2);
Pop.R2y{1,3} = var_swap(Pop.Ry2{1,2}.',ss2,tt2);


end


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
degmat_sort = sortrows_integerTable(degmat_sort);   % sort monomials from low to high degree
Z_degmat = fliplr(degmat_sort(:,2:end));
Z = polynomial(eye(nZ),Z_degmat,Zvarname,[nZ,1]);

end



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