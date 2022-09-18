function [sos,Deop] = lpi_ineq_2d(sos,Pop, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sos = lpi_ineq_2d(sos,P, options) sets up inequality constraints such that P>=0
% 
% INPUT
%   prog: SOS program to modify.
%   P: PI dopvar2d variable
%   -Optional INPUTS
%   options: determines if the operator needs to be pure or to include psatz
%   term. 
%   If options.pure(1)=1, then Rxx^o term is excluded. 
%   If options.pure(2)=1, then Ryy^o term is excluded. 
%   If options.pure(3)=1, then R22^oo term is excluded. 
%   If options.pure(4)=1, then R22^ao and R22^bo terms are excluded. 
%   If options.pure(5)=1, then R22^oa and R22^ob terms are excluded. 
%   If options.psatz=1, then the psatz multiplier is included.
% 
% OUTPUT 
%   sos: SOS program
%   Dop: dopvar_2d object defining the LPI equality Pop = Dop, where Dop>=0
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - lpi_ineq_2d
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 07_21_2021 
% 02/21/2022 - DJ: Update to allow multiple psatz terms

% Extract the inputs
switch nargin
    case 1
        error(['Not enough inputs!'])
    case 2
        % If no options are specified, use no special options
        options.psatz = 0;
        options.pure = zeros(1,5);
        options.sep = zeros(1,6);
    case 3
        % If a certain option is not specified, do no use it
        if ~isfield(options,'psatz')
            options.psatz=0;
        end
        if ~isfield(options,'pure')
            options.pure = zeros(1,5);
        end
        if ~isfield(options,'sep')
            options.sep = zeros(1,6);
        end
end

dim = Pop.dim;
if dim(:,1)~=dim(:,2)
    error('Non-symmetric Operators cannot be sign definite. Unable to set the inequality');
end

% We're going to enforce inequality P>=0, by building an operator Deop>=0,
% and setting P=Deop. If psatz is used, build an additional Deop2>=0, and
% set P=Deop+Deop2.

% First, extract the dimensions and domain from P
n0 = dim(1,1);
nx = dim(2,1);
ny = dim(3,1);
n2 = dim(4,1);
dom = Pop.I;

% Next, exclude certain elements from Deop (options2) and Deop2 (options3)
% if specified
options2.exclude = zeros(1,16);
options3.exclude = zeros(1,16);

options2.exclude(2) = options.pure(1);
options3.exclude(2) = options.pure(1); 

options2.exclude(5) = options.pure(2);
options3.exclude(5) = options.pure(2);

tol = 1e-14;
if (isa(Pop.R22{1,1},'double') && max(max(Pop.R22{1,1}))<tol) || all(max(max(Pop.R22{1,1}.C))<tol)
    options2.exclude(8) = 1;
    options3.exclude(8) = 1;
else
    options2.exclude(8) = max([options.pure(3),options.pure(4),options.pure(5)]);
    options3.exclude(8) = max([options.pure(3),options.pure(4),options.pure(5)]);
end
if (isa(Pop.R22{2,1},'double') && max(max(Pop.R22{2,1}))<tol) || all(max(max(Pop.R22{2,1}.C))<tol)
    options2.exclude(9) = 1;
    options3.exclude(9) = 1;
else
    options2.exclude(9) = options.pure(4);
    options3.exclude(9) = options.pure(4); 
end
if (isa(Pop.R22{3,1},'double') && max(max(Pop.R22{3,1}))<tol) || all(max(max(Pop.R22{3,1}.C))<tol)
    options2.exclude(10) = 1;
    options3.exclude(10) = 1;
else
    options2.exclude(10) = options.pure(4);
    options3.exclude(10) = options.pure(4); 
end
if (isa(Pop.R22{1,2},'double') && max(max(Pop.R22{1,2}))<tol) || all(max(max(Pop.R22{1,2}.C))<tol)
    options2.exclude(11) = 1;
    options3.exclude(11) = 1;
else
    options2.exclude(11) = options.pure(5);
    options3.exclude(11) = options.pure(5); 
end
if (isa(Pop.R22{1,3},'double') && max(max(Pop.R22{1,3}))<tol) || all(max(max(Pop.R22{1,3}.C))<tol)
    options2.exclude(12) = 1;
    options3.exclude(12) = 1;
else
    options2.exclude(12) = options.pure(5);
    options3.exclude(12) = options.pure(5); 
end

% Now, determine appropriate degrees for Deop and Deop2 to match degrees
% for Pop. Two options:
% 1. NOT RECOMMENDED
%   Establish maximal degrees of the monomials defining Deop such that the
%   maximal degrees of the monomials appearing in Deop matches the maximal
%   degree of the monomials in Pop. This option does not guarantee all the
%   necessary monomials appear, so the maximal degree may be increased
%   later on. Slightly safer, but more expensive!
% 2. 
%   Establish a specific set of monomials used to define Deop. In this
%   case, we decompose Deop into four terms: Deop = Deop_xy + Deop_oy +
%   Deop_xo + Deop_oo. Here, Deop_oy has separable kernels in the
%   x-integrals (e.g. Deop_oy.Rxx{2}=Deop_oy.Rxx{3}), Deop_oy has separable
%   kernels in the y-integrals (e.g. Deop_xo.Ryy{2}=Deop_xo.Ryy{3}), and
%   Deop_oo has separable kernels in both the x and y integrals (e.g. 
%   Deop_oo.R22{2,2}=Deop_oo.R22{3,2}=Deop_oo.R22{2,3}=Deop_oo.R22{3,3}).
toggle = 2;
if toggle==1
    degs = degbalance(Pop,options2);
elseif toggle==2
    [degs,degs_oy,degs_xo,degs_oo] = getmonomials_lpi_ineq(Pop,options2);
end

% Finally, we can construct the positive operators Deop and Deop2, and use
% them to enforce positivity of P
[sosD, Deop] = poslpivar_2d(sos, [n0, nx, ny, n2],dom,degs,options2);

% % Add aditional degrees of freedom/monomials to obtain a (hopefully)
% % feasible problem
if toggle==1    % Using maximal degrees
    % % Check that the parameters of Deop indeed contain all monomials that appear
    % % in the parameters of Pop
    p_indx = [2;3;4;6;7;8;11;12;16];    % Check only lower-triangular parameters
    [isgood_Deop,isgood_Dpar,degs] = checkdeg_lpi_eq_2d(Pop,Deop,degs,p_indx);

    % If Deop is missing monomials, keep increasing the degrees until Deop has
    % all the necessary monomials
    while ~isgood_Deop
        warning('The specified options for ''lpi_ineq'' do not allow sufficient freedom to enforce the inequality constraint. Additional monomials are being added.')

        % Construct a new positive operator Deop with greater degrees
        [sosD, Deop] = poslpivar_2d(sos, [n0, nx, ny, n2],dom,degs,options2);

        % Check that now Deop has all the necessary monomials
        p_indx = find(~(isgood_Dpar(:)));   % Indices of parameters we still need to verify are okay
        [isgood_Deop,isgood_Dpar,degs] = checkdeg_lpi_eq_2d(Pop,Deop,degs,p_indx);    
    end
elseif toggle==2    % Using predefined monomials
    % % Add addtional terms with separable kernels
    if any(~options2.exclude([3,4,9,10,13,14,15,16]))
        options2.sep = [1,0,1,0,1,0];
        [sosD, Deop_x] = poslpivar_2d(sosD, [n0, nx, ny, n2],dom,degs_oy,options2);
    end
    if any(~options2.exclude([6,7,11,12,13,14,15,16]))
        options2.sep = [0,1,0,1,0,1];
        [sosD, Deop_y] = poslpivar_2d(sosD, [n0, nx, ny, n2],dom,degs_xo,options2);
    end
    if any(~options2.exclude([3,4,6,7,9,10,11,12,13,14,15,16]))
        options2.sep = [1,1,1,1,1,1];
        [sosD, Deop_xy] = poslpivar_2d(sosD, [n0, nx, ny, n2],dom,degs_oo,options2);
    end
    Deop = Deop + Deop_x + Deop_y + Deop_xy;
end
sos = sosD; % Make sure the SOS program contains the right operator Deop
    
% If desired, add a Psatz term to enforce only local positivity
for j=1:length(options.psatz)
    if options.psatz(j)~=0
        options3.psatz = options.psatz(j);
        [sos, De2op] = poslpivar_2d(sos, [n0 ,nx, ny, n2],dom,degs,options3);
        Deop = Deop+De2op; 
    end
end
% Enforce the constraint
sos = lpi_eq_2d(sos,Deop-Pop);  %Pop==Deop

end




function [degs_xy,degs_oy,degs_xo,degs_oo] = getmonomials_lpi_ineq(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [degs_xy,degs_oy,degs_xo,degs_oo] = getmonomials_lpi_ineq(P,opts) 
% provides monomial degrees for "poslpivar" to construct dopvar2d objects
% Q_xy, Q_oy, Q_xo and Q_oo such that the sum Q = Q_xy + Q_oy + Q_xo + Q_oo 
% is defined by similar monomials as the dopvar2d object P.
% 
% INPUT
%   P: PI dopvar2d class object, for which we want to impose some LPI
%   inequality P = Q with Q = Q_xy + Q_oy + Q_xo + Q_oo >=0
%   opts: struct specifying restrictions on the opvar2d object Q.
%            In particular, opts should have a field "exclude", specifying
%            which monomials can be excluded in the construction of the
%            positive operator Q. See also the header for "poslpivar_2d".
% 
% OUTPUT 
%   degs_xy: A struct with fields Zx, Zy and Z2, each of which is a cell
%       object, describing which monomials to use in constructing the
%       positive semidefinite operator Q_xy.
%   degs_oy: A struct with same fields as degs_xy, specifying the degrees
%       of the monomials used to construct Q_oy. Here, Q_oy is constructed
%       with separable kernels for the x-integrals, so that e.g.
%       Q_oy.Rxx{2} = Q_oy.Rxx{3} and Q_oy.R22{2,2} = Q_oy.R22{3,2}.
%   degs_xo: A struct with same fields as degs_xy, specifying the degrees
%       of the monomials used to construct Q_xo. Here, Q_xo is constructed
%       with separable kernels for the y-integrals, so that e.g.
%       Q_xo.Ryy{2} = Q_xo.Ryy{3} and Q_xo.R22{2,2} = Q_xo.R22{2,3}.
%   degs_oo: A struct with same fields as degs_xy, specifying the degrees
%       of the monomials used to construct Q_oo. Here, Q_oo is constructed
%       with separable kernels for both the x- and y-integrals, so that e.g.
%       Q_oo.R22{2,2} = Q_oo.R22{2,3} = Q_oo.R22{3,2} = Q_oo.R22{3,3}.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - getmonomials_lpi_ineq
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
% Initial coding DJ - 04_28_2022

% Process the inputs
if nargin==1
    Zexclude = zeros(1,16);
elseif nargin==2
    if isfield(opts,'exclude')
        Zexclude = opts.exclude;
    else
        Zexclude = zeros(1,16);
    end
end

% Extract the dimensions of the operator P
n = P.dim(:,1);
if any(P.dim(:,1)~=P.dim(:,2))
    error('The input operator must be square (symmetric) in order to impose LPI (in)equality constraints')
end
% If the operator does not map to/from certain states, there is no need to
% implement associated monomials
if n(1)==0
    Zexclude(1)=1;
end
if n(2)==0
    Zexclude(2:4)=1;
end
if n(3)==0
    Zexclude(5:7)=1;
end
if n(4)==0
    Zexclude(8:16)=1;
end

% % Some other settings:
% We will construct the monomials based on joint degrees in x*tt, and y*nu.
% Then, we will allow all combinations of degrees in e.g. x and tt that add
% up to the allowed max degree in x*tt. The following values specify a
% minimal degree in each of the variables in the different monomials,
% including only monomials with a degree in each variable that is as least
% as big as specified. Will decrease the freedom if not 0, and thus reduce
% complexity, but is not recommended.
xdeg_min_Zx = 0;    % Minimal degree of x in Zx
ttdeg_min_Zx = 0;   % Minimal degree of tt in Zx
ydeg_min_Zy = 0;    % Minimal degree of y in Zy
nudeg_min_Zy = 0;   % Minimal degree of nu in Zy
xdeg_min_Z2 = 0;    % Minimal degree of x in Z2
ttdeg_min_Z2 = 0;   % Minimal degree of tt in Z2
ydeg_min_Z2 = 0;    % Minimal degree of y in Z2
nudeg_min_Z2 = 0;   % Minimal degree of nu in Z2

% For the operators Q_oy, Q_xo and Q_xy, with separable kernels, the
% separability reduces the freedom in the choice of the operator, freedom
% which we tend to need. To account for this, we scale the number of
% monomials defining e.g. Q_oy.Rxx{2}=Q_oy.Rxx{3} with a factor of
% "fctr_x", aiming to retrieve the lost freedom. Default values are
% sqrt(2), based on the idea that we wish to double number of decision 
% variables, and the number of decision variables corresponds roughly to 
% the number of monomials squared. Slightly heuristic though, you can
% increase for security or decrease for speed.
fctr_x_Zx = sqrt(2);
fctr_y_Zy = sqrt(2);
fctr_x_Z2 = sqrt(2);
fctr_y_Z2 = sqrt(2);

% Extract the degrees of the variables appearing in each monomial
% in the operator
Pdegs = getdeg(P);

degs_xy = struct();     degs_xo = struct();
degs_oy = struct();     degs_oo = struct();
degs_xy.use_monomials = 1;  degs_xo.use_monomials = 1;
degs_oy.use_monomials = 1;  degs_oo.use_monomials = 1;

% % % =============================================================== % % %
% % % Compute monomials for Zx
Zx_x = {[];[];[]};      Zx_o = {[];[];[]};
xdeg_min = xdeg_min_Zx;
ttdeg_min = ttdeg_min_Zx;

if ~Zexclude(2)
% % Start with multiplier: P.Rxx{1} = Zx{1}(x)' * Qxx{1,1} * Zx{1}(x)
% % Degrees of Zx{1} are just halve degree of P.Rxx{1}
mindeg_x = floor(min(Pdegs.Rxx{1}(:,1))/2);
maxdeg_x = ceil(max(Pdegs.Rxx{1}(:,1))/2);
Zx_x{1} = makesparse((mindeg_x : maxdeg_x)');
Zx_o{1} = Zx_x{1};
else
    Zx_x{1} = zeros(0,1);
    Zx_o{1} = zeros(0,1);
end

if ~Zexclude(3) || ~Zexclude(4)
    % % For integrators, focus on terms involving product of monomial vectors
    % % with themselves:
    %
    % P.Rxx{2} = ... + int_x^b Zx{2}(eta,x)' * Qxx{2,2} * Zx{2}(eta,tt) d eta
    %           +... + int_a^tt Zx{3}(eta,x)' * Qxx{3,3} * Zx{3}(eta,tt) d eta
    %
    % % We note that each integral produces two terms: a limit that varies in x
    % % or tt, and a constant limit. We focus first on the varying limits.
    
    % For the varying limits, the joint degrees in x*tt of the result are given
    % by the sum of the joint degrees in x*tt of the monomials Zx{i}(x,tt),
    % plus 1 as a result of the integration. Extract first these joint degrees:
    Pdmat = [Pdegs.Rxx{2}(:,[1,2]);Pdegs.Rxx{3}(:,[1,2])];    % All degrees in [x, tt]
    joint_degs = sum(Pdmat,2)-1;                            % Joint degrees in x*tt, subtracting contribution from integral
    
    % Next, establish the necessary joint degrees in x*tt of the monomials
    % Zx{i}, so that their product returns all the necessary joint degrees:
    mindeg_xtt = max(floor(min(joint_degs)/2),0);
    maxdeg_xtt = floor(max(joint_degs)/2);
    Zd_joint = (mindeg_xtt : maxdeg_xtt)';
    
    % Finally, for each joint degree, establish all combinations of degrees in
    % x and tt that add up to this joint degree:
    nZd_x = sum(max(Zd_joint+1-xdeg_min-ttdeg_min,1));
    Zd_x = zeros(nZd_x,2);
    indx_n = 0;
    for k=1:size(Zd_joint,1)
        Zdk_xtt = Zd_joint(k);                                  % Required joint degree in x*tt
        i1x = min(xdeg_min,Zdk_xtt);
        i1tt = min(ttdeg_min,Zdk_xtt-i1x);
        Zddk = [(i1x:Zdk_xtt-i1tt)',(Zdk_xtt-i1x:-1:i1tt)'];    % Degrees in x,tt that add up to Zdk_xtt
        
        % Add to the list of monomials
        indx = indx_n;
        indx_n = indx + size(Zddk,1);
        Zd_x(indx+1:indx_n,:) = Zddk;
    end
    Zd_x = makesparse(Zd_x);
    Zx_x{2} = Zd_x;     Zx_x{3} = Zd_x;
    
    % % With these monomials Zx_x, all the monomials that appear in Pop.Rx{2}
    % % and Pop.Rx{3} will be produced as a result of the varying limits of the
    % % integrals
    % P.Rxx{2} = ... + int_x^b Zx_x{2}(eta,x)' * Qxx_x{2,2} * Zx_x{2}(eta,tt) d eta
    %           +... + int_a^tt Zx_x{3}(eta,x)' * Qxx_x{3,3} * Zx_x{3}(eta,tt) d eta
    % % However, evaluating these integrals at the lower limits introduces
    % % additional terms, imposing undesirable constraints on our decision
    % % variables. We have to introduce additional freedom so that these
    % % additional terms can be canceled if necessary, or conversely, the
    % % varying limits of the integrals can be discarded. To this end, we add a
    % % term to the operator with separable kernels:
    % P.Rxx{2} = ... + int_a^b Zx_o{2}(eta,x)' * Qxx_o{2,2} * Zx_o{2}(eta,tt) d eta
    % % Note that the degrees in x of the monomials Zx_o{2} do not contribute
    % % to the degrees of P.Rxx{2}. As such, we can include as many monomials
    % % as needed to get the required freedom.
    
    % First, establish what degrees of tt we want. Note that the degrees in tt
    % of the monomials define the joint degrees in x*tt of the integral of the
    % product. We also try to reproduce the monomials that the non-separable
    % operator may already produce, for which we reintroduce the additional
    % degree from the integration that was previously subtracted.
    Zd_o_tt = [Zd_joint; maxdeg_xtt+1];
    nd_tt = length(Zd_o_tt);
    % Zd_o should include more monomials, to account for the separability of
    % the kernel. We use as many degrees in x as necessary to achieve this,
    % noting that the x variable vanishes in the separable case.
    nZd_o = ceil(fctr_x_Zx*nZd_x);
    nd_x = ceil(nZd_o/nd_tt);
    Zd_o_x = (0 : nd_x-1)';
    % Consider all combinations of x and tt
    Zd_o = [repmat(Zd_o_x,[nd_tt,1]), kron(Zd_o_tt,ones(nd_x,1))];
    % Finally, store the degrees
    Zd_o = makesparse(Zd_o);
    Zx_o{2} = Zd_o;     Zx_o{3} = Zd_o;
    
else
    Zx_x{2} = zeros(0,2);      Zx_x{3} = zeros(0,2);
    Zx_o{2} = zeros(0,2);      Zx_o{3} = zeros(0,2);
end
% Separability of the y integrals does not affect which monomials we use in
% Zx
degs_xy.Zx = Zx_x;      degs_xo.Zx = Zx_x;
degs_oy.Zx = Zx_o;      degs_oo.Zx = Zx_o;


% % % =============================================================== % % %
% % % Compute monomials for Zy
Zy_y = {[],[],[]};     Zy_o = {[],[],[]};
ydeg_min = ydeg_min_Zy;
nudeg_min = nudeg_min_Zy;

if ~Zexclude(5)
% % Start with multiplier: P.Ryy{1} = Zy{1}(y)' * Qyy{1,1} * Zy{1}(y)
% % Degrees of Zx{1} are just halve degree of P.Rxx{1}
mindeg_y = floor(min(Pdegs.Ryy{1}(:,3))/2);
maxdeg_y = ceil(max(Pdegs.Ryy{1}(:,3))/2);
Zy_y{1} = makesparse((mindeg_y : maxdeg_y)');
Zy_o{1} = Zy_y{1};
else
    Zy_y{1} = zeros(0,1);
    Zy_o{1} = zeros(0,1);
end

if ~Zexclude(6) || ~Zexclude(7)
    % % Next, work on integrators:
    % P.Ryy{2} = ... + int_y^d Zy{2}(mu,y)' * Qyy{2,2} * Zy{2}(mu,nu) d mu
    %           +... + int_c^nu Zy{3}(mu,y)' * Qyy{3,3} * Zy{3}(mu,nu) d mu
    % % We note that each integral produces two terms: a limit that varies in y
    % % or nu, and a constant limit. We focus first on the varying limits.

    % For the varying limits, the joint degrees in y*nu of the result are given
    % by the sum of the joint degrees in y*nu of the monomials Zy{i}(y,nu),
    % plus 1 as a result of the integration. Extract first these joint degrees:
    Pdmat = [Pdegs.Ryy{2}(:,[3,4]);Pdegs.Ryy{3}(:,[3,4])];    % All degrees in [y, nu]
    joint_degs = sum(Pdmat,2)-1;            % Joint degrees in y*nu, subtracting contribution from integral

    % Next, establish the necessary joint degrees in y*nu of the monomials
    % Zy{i}, so that their product returns all the necessary joint degrees:
    mindeg_ynu = max(floor(min(joint_degs)/2),0);
    maxdeg_ynu = floor(max(joint_degs)/2);
    Zd_joint = (mindeg_ynu : maxdeg_ynu)';

    % Finally, for each joint degree, establish all combinations of degrees in y
    % and nu that add up to this joint degree:
    nZd_y = sum(max(Zd_joint+1-ydeg_min-nudeg_min,1));  % Number of possible monomials in (x,tt) with joint degrees in (x*tt) as in Zd_joint
    Zd_y = zeros(nZd_y,2);
    indx_n = 0;
    for k=1:size(Zd_joint,1)
        Zdk_ynu = Zd_joint(k);                                  % Required joint degree in y*nu
        i1y = min(ydeg_min,Zdk_ynu);
        i1nu = min(nudeg_min,Zdk_ynu-i1y);
        Zddk = [(i1y:Zdk_ynu-i1nu)',(Zdk_ynu-i1y:-1:i1nu)'];    % Degrees in y,nu that add up to Zdk_ynu

        % Add to the list of monomials
        indx = indx_n;
        indx_n = indx + size(Zddk,1);
        Zd_y(indx+1:indx_n,:) = Zddk;
    end
    Zd_y = makesparse(Zd_y);
    Zy_y{2} = Zd_y;     Zy_y{3} = Zd_y;

    % % Next, account for the constant limits of the integrals
    % The degrees of nu should produce the desired joint degrees
    Zd_o_nu = [Zd_joint; maxdeg_ynu+1];
    nd_nu = length(Zd_o_nu);
    % Zd_o should include more monomials, to account for the separability of
    % the kernel. We use as many degrees in y as necessary to achieve this,
    % noting that the x variable vanishes in the separable case.
    nZd_o = ceil(fctr_y_Zy*nZd_y);
    nd_y = ceil(nZd_o/nd_nu);
    Zd_o_y = (0 : nd_y-1)';
    % Consider all combinations of y and nu
    Zd_o = [repmat(Zd_o_y,[nd_nu,1]), kron(Zd_o_nu,ones(nd_y,1))];
    % Finally, store the degrees
    Zd_o = makesparse(Zd_o);
    Zy_o{2} = Zd_o;     Zy_o{3} = Zd_o;

else
    Zy_y{2} = zeros(0,2);      Zy_y{3} = zeros(0,2);
    Zy_o{2} = zeros(0,2);      Zy_o{3} = zeros(0,2);
end
% Separability of the x integrals does not affect which monomials we use in
% Zy
degs_xy.Zy = Zy_y;      degs_xo.Zy = Zy_o;
degs_oy.Zy = Zy_y;      degs_oo.Zy = Zy_o;


% % % =============================================================== % % %
% % % Compute monomials Z2
Z2_xy = {[],[],[];[],[],[];[],[],[]};   Z2_xo = {[],[],[];[],[],[];[],[],[]};
Z2_oy = {[],[],[];[],[],[];[],[],[]};   Z2_oo = {[],[],[];[],[],[];[],[],[]};
xdeg_min = xdeg_min_Z2;
ttdeg_min = ttdeg_min_Z2;
ydeg_min = ydeg_min_Z2;
nudeg_min = nudeg_min_Z2;

if ~Zexclude(8)
    % % Start with multiplier: P.R22{1,1} = Z2{1,1}(x,y)' * Q22{1,1} * Z2{1,1}(x,y)
    % % Degrees of Z2{1,1} are established using Newton polytope
    Z2_xy{1,1} = getSOSmonomials(Pdegs.R22{1,1}(:,[1,3]),1);
    Z2_xo{1,1} = Z2_xy{1,1};
    Z2_oy{1,1} = Z2_xy{1,1};
    Z2_oo{1,1} = Z2_xy{1,1};
else
    Z2_xy{1,1} = zeros(0,2);    Z2_xo{1,1} = zeros(0,2);
    Z2_oy{1,1} = zeros(0,2);    Z2_oo{1,1} = zeros(0,2);
end

% % % --------------------------------------------------------------- % % %
if ~Zexclude(9) || ~Zexclude(10)
    % % Next, integrators along the x-direction:
    % P.R22{2,1} = ... + int_x^b Z2{2,1}(eta,x,y)' * Q22{2,2} * Z2{2,1}(eta,tt,y) d eta
    %             +... + int_a^tt Z2{3,1}(eta,x,y)' * Q22{3,3} * Z2{3,1}(eta,tt,y) d eta
    
    % For the varying limits, the joint degrees in x*tt of the result are given
    % by the sum of the joint degrees in x*tt of the monomials Z2{i,1}(x,tt,y),
    % plus 1 as a result of the integration. Extract first the joint degrees,
    % and establish the necessary degrees in [x*tt,y] of the monomials
    % Z2{i,1}, so that their product returns all the necessary joint degrees:
    Pdmat = [Pdegs.R22{2,1}(:,[1,2,3]);Pdegs.R22{3,1}(:,[1,2,3])];  % All degrees in [x, tt, y]
    joint_degs = [sum(Pdmat,[1,2])-1,Pdmat(:,3)];                   % Degrees in [x*tt,y], subtracting contribution from integral
    joint_degs = joint_degs(all(joint_degs>=0,2),:);
    joint_degs = unique(joint_degs,'rows');
    if isempty(joint_degs)
        joint_degs = zeros(1,size(joint_degs,2));
    end
    Zd_joint = getSOSmonomials(joint_degs,1);
    
    % Then, for each joint degree, establish all combinations of degrees in x
    % and tt that add up to this joint degree:
    nZd_x = sum(max(Zd_joint(:,1)+1-xdeg_min-ttdeg_min,1));  % Number of possible monomials in (x,tt,y) with joint degrees in (x*tt,y) as in Zd_joint
    Zd_x = zeros(nZd_x,3);
    indx_n = 0;
    for k=1:size(Zd_joint,1)
        Zdk_xtt = Zd_joint(k,1);                                    % Required joint degree in x*tt
        i1x = min(xdeg_min,Zdk_xtt);
        i1tt = min(ttdeg_min,Zdk_xtt-i1x);
        Zddk_xtt = [(i1x:Zdk_xtt-i1tt)',(Zdk_xtt-i1x:-1:i1tt)'];    % Degrees in x,tt that add up to Zdk_xtt
        Zddk = [Zddk_xtt, Zd_joint(k,2)*ones(size(Zddk_xtt,1),1)];  % Add degrees in y
        
        % Add to the list of monomials
        indx = indx_n;
        indx_n = indx + size(Zddk,1);
        Zd_x(indx+1:indx_n,:) = Zddk;
    end
    Zd_x = makesparse(Zd_x);
    Z2_xy{2,1} = Zd_x;      Z2_xy{3,1} = Zd_x;
    Z2_xo{2,1} = Zd_x;      Z2_xo{3,1} = Zd_x;
    
    
    % % Next, introduce additional monomials to account for the constant limits
    % % of the integrals.
    % The joint degree in x*tt is completely determined by the degree in tt
    joint_degs_o = joint_degs(joint_degs(:,1)==max(joint_degs(:,1)),:) + [1,0];
    joint_degs_o = [joint_degs; joint_degs_o];
    Zd_o_tty = getSOSmonomials(joint_degs_o,1);
    nd_tty = size(Zd_o_tty,1);
    % Allow all degrees in x necessary to obtain sufficient monomials
    nZd_o = ceil(fctr_x_Zx*nZd_x);
    nd_x = ceil(nZd_o/nd_tty);
    Zd_o_x = (0 : nd_x-1)';
    % Consider all appropriate combinations of [x,tt,y]
    Zd_o = [repmat(Zd_o_x,[nd_tty,1]), kron(Zd_o_tty,ones(nd_x,1))];
    % Finally, store the degrees
    Zd_o = makesparse(Zd_o);
    Z2_oy{2,1} = Zd_o;      Z2_oy{3,1} = Zd_o;
    Z2_oo{2,1} = Zd_o;      Z2_oo{3,1} = Zd_o;
    
else
    Z2_xy{2,1} = zeros(0,3);   Z2_xy{3,1} = zeros(0,3);
    Z2_oy{2,1} = zeros(0,3);   Z2_oy{3,1} = zeros(0,3);
    Z2_xo{2,1} = zeros(0,3);   Z2_xo{3,1} = zeros(0,3);
    Z2_oo{2,1} = zeros(0,3);   Z2_oo{3,1} = zeros(0,3);
end

% % % --------------------------------------------------------------- % % %
if ~Zexclude(11) || ~Zexclude(12)
    % % Next, integrators along the y-direction:
    % P.R22{1,2} = ... + int_y^d Z2{1,2}(x,mu,y)' * Q22{4,4} * Z2{1,2}(x,mu,nu) d mu
    %             +... + int_c^nu Z2{1,3}(x,mu,y)' * Q22{5,5} * Z2{1,3}(x,mu,nu) d mu
    
    % % Next, integrators along the x-direction:
    % P.R22{2,1} = ... + int_x^b Z2{2,1}(eta,x,y)' * Q22{2,2} * Z2{2,1}(eta,tt,y) d eta
    %             +... + int_a^tt Z2{3,1}(eta,x,y)' * Q22{3,3} * Z2{3,1}(eta,tt,y) d eta
    
    % For the varying limits, the joint degrees in y*nu of the result are given
    % by the sum of the joint degrees in y*nu of the monomials Z2{1,j}(x,y,nu),
    % plus 1 as a result of the integration. Extract first the joint degrees,
    % and establish the necessary degrees in [x,y*nu] of the monomials
    % Z2{1,j}, so that their product returns all the necessary joint degrees:
    Pdmat = [Pdegs.R22{1,2}(:,[1,3,4]);Pdegs.R22{1,3}(:,[1,3,4])];     % All degrees in [x, y, nu]
    joint_degs = [Pdmat(:,1),sum(Pdmat(:,[2,3]),2)-1];                  % Degrees in [x,y*nu], subtracting contribution from integral
    joint_degs = joint_degs(all(joint_degs>=0,2),:);
    joint_degs = unique(joint_degs,'rows');
    if isempty(joint_degs)
        joint_degs = zeros(1,size(joint_degs,2));
    end
    Zd_joint = getSOSmonomials(joint_degs,1);
    
    % Then, for each joint degree, establish all combinations of degrees in y
    % and nu that add up to this joint degree:
    nZd_y = sum(max(Zd_joint(:,2)+1-ydeg_min-nudeg_min,1));  % Number of possible monomials in (x,y,nu) with joint degrees in (x,y*nu) as in Zd_joint
    Zd_y = zeros(nZd_y,3);
    indx_n = 0;
    for k=1:size(Zd_joint,1)
        Zdk_ynu = Zd_joint(k,2);                                    % Required joint degree in u*nu
        i1y = min(ydeg_min,Zdk_ynu);
        i1nu = min(nudeg_min,Zdk_ynu-i1y);
        Zddk_ynu = [(i1y:Zdk_ynu-i1nu)',(Zdk_ynu-i1y:-1:i1nu)'];    % Degrees in y,nu that add up to Zdk_ynu
        Zddk = [Zd_joint(k,1)*ones(size(Zddk_ynu,1),1), Zddk_ynu];  % Add degrees in x
        
        % Add to the list of monomials
        indx = indx_n;
        indx_n = indx + size(Zddk,1);
        Zd_y(indx+1:indx_n,:) = Zddk;
    end
    Zd_y = makesparse(Zd_y);
    Z2_xy{1,2} = Zd_y;      Z2_xy{1,3} = Zd_y;
    Z2_oy{1,2} = Zd_y;      Z2_oy{1,3} = Zd_y;
    
    
    % % Next, introduce additional monomials to account for the constant limits
    % % of the integrals.
    % The joint degree in y*nu is completely determined by the degree in nu
    joint_degs_o = joint_degs(joint_degs(:,2)==max(joint_degs(:,2)),:) + [0,1];
    joint_degs_o = [joint_degs; joint_degs_o];
    Zd_o_xnu = getSOSmonomials(joint_degs_o,1);
    nd_xnu = size(Zd_o_xnu,1);
    % Allow all degrees in y necessary to obtain sufficient monomials
    nZd_o = ceil(fctr_y_Z2*nZd_y);
    nd_y = ceil(nZd_o/nd_xnu);
    Zd_o_y = (0 : nd_y-1)';
    % Consider all combinations of [x,y,nu]
    Zd_o = [repmat(Zd_o_xnu(:,1),[nd_y,1]), kron(Zd_o_y,ones(nd_xnu,1)), repmat(Zd_o_xnu(:,2),[nd_y,1])];
    % Finally, store the degrees
    Zd_o = makesparse(Zd_o);
    Z2_xo{1,2} = Zd_o;      Z2_xo{1,3} = Zd_o;
    Z2_oo{1,2} = Zd_o;      Z2_oo{1,3} = Zd_o;
    
else
    Z2_xy{1,2} = zeros(0,3);   Z2_xy{1,3} = zeros(0,3);
    Z2_xo{1,2} = zeros(0,3);   Z2_xo{1,3} = zeros(0,3);
    Z2_oy{1,2} = zeros(0,3);   Z2_oy{1,3} = zeros(0,3);
    Z2_oo{1,2} = zeros(0,3);   Z2_oo{1,3} = zeros(0,3);
end

% % % --------------------------------------------------------------- % % %
if any(~Zexclude(13:16))
    % % Finally, integrators along both directions:
    % P.R22{1,2} = ... + int_x^b int_y^d Z2{2,2}(eta,x,mu,y)' * Q22{6,6} * Z2{2,2}(eta,tt,mu,nu) d mu d eta
    %             +... + int_a^tt int_y^d Z2{3,2}(eta,x,mu,y)' * Q22{7,7} * Z2{3,2}(eta,tt,mu,nu) d mu d eta
    %             +... + int_x^b int_c^nu Z2{2,3}(eta,x,mu,y)' * Q22{8,8} * Z2{2,3}(eta,tt,mu,nu) d mu d eta
    %             +... + int_a^tt int_c^nu Z2{3,3}(eta,x,mu,y)' * Q22{9,9} * Z2{3,3}(eta,tt,mu,nu) d mu d eta
    
    % For the varying limits, the joint degrees in x*tt and y*nu of the result
    % are given by the sum of the joint degrees in x*tt and y*nu of the
    % monomials Z2{i,j}(x,tt,y,nu), plus 1 as a result of the integration.
    % Extract first the joint degrees, and establish the necessary degrees in
    % [x*tt,y*nu] of the monomials Z2{i,j}, so that their product returns all
    % the necessary joint degrees:
    Pdmat = [Pdegs.R22{2,2}; Pdegs.R22{3,2}; Pdegs.R22{2,3}; Pdegs.R22{3,3}];
    % Sum degrees of x and tt, and of y and nu
    joint_degs = unique(Pdmat * kron(speye(2),ones(2,1)),'rows');
    joint_degs = joint_degs-[1,1];
    joint_degs = joint_degs(all(joint_degs>=0,2),:);
    if isempty(joint_degs)
        joint_degs = zeros(1,size(joint_degs,2));
    end
    Zd_xy = getSOSmonomials(joint_degs,1);
    nZd_xy = size(Zd_xy,1);
    
    % Then, for each joint degree, establish all combinations of degrees in x
    % and tt and y and nu that add up to this joint degree:
    nZdd_Z = prod(max(Zd_xy+1-[xdeg_min+ttdeg_min,ydeg_min+nudeg_min],1),2); % Number of monomials that produce each desired joint degree
    nZdd_xy = sum(nZdd_Z,1); % Number of monomials in (x,tt,y,nu) with degrees in (x*tt) and (y*nu) as in Zd_xy
    Zdd_xy = zeros(nZdd_xy,4);
    indx_n = 0;
    for k=1:nZd_xy
        Zdk_xtt = Zd_xy(k,1);                          % Required degree in x*tt
        i1x = min(xdeg_min,Zdk_xtt);
        i1tt = min(ttdeg_min,Zdk_xtt-i1x);
        Zddk_xtt = [(i1x:Zdk_xtt-i1tt)',(Zdk_xtt-i1x:-1:i1tt)'];  % Degrees in x,tt that add up to Zdk_xtt
        
        Zdk_ynu = Zd_xy(k,2);          % Required degree in x*tt
        i1y = min(ydeg_min,Zdk_ynu);
        i1nu = min(nudeg_min,Zdk_ynu-i1y);
        Zddk_ynu = [(i1y:Zdk_ynu-i1nu)',(Zdk_ynu-i1y:-1:i1nu)'];  % Degrees in y,nu that add up to Zdk_ynu
        
        % Consider all possible pairings of the degrees in (x,tt) and (y,nu)
        Zddk = [repmat(Zddk_xtt,[Zdk_ynu+1-i1y-i1nu,1]), kron(Zddk_ynu,ones(Zdk_xtt+1-i1x-i1tt,1))];
        
        % Add to the list of monomials
        indx = indx_n;
        indx_n = indx + size(Zddk,1);
        Zdd_xy(indx+1:indx_n,:) = Zddk;
    end
    Zdd_xy = makesparse(Zdd_xy);
    Z2_xy{2,2} = Zdd_xy;        Z2_xy{2,3} = Zdd_xy;
    Z2_xy{3,2} = Zdd_xy;        Z2_xy{3,3} = Zdd_xy;
    
    % The monomials Z2_xy should be sufficient to produce all the monomials
    % that appear in the operator using integration up to x and y (or tt and
    % nu). However, each of the products of these monomials in Zd_xy will also
    % be evaluted at the constant limit of the integrals, so we need to include
    % more freedom in the program to account for this.
    
    % % First, introduce monomials for separable kernel in x-integral. We choose
    % % these monomials such that we can reproduce all the monomials produced
    % % by the non-separable kernels, as well as any other monomial that
    % % appears in the system.
    % Establish a set of joint degrees in [x*tt,y*nu]
    joint_degs_oy = joint_degs(joint_degs(:,1)==max(joint_degs(:,1)),:) + [1,0];    % Reintroduce previously substracted degree
    joint_degs_oy = [joint_degs;joint_degs_oy];
    % Establish degrees for the monomials in Z2_oy. Note: Using a separable
    % kernel in the x-integral, the joint degree in x*tt of the product will
    % be completely determined by the degree in just tt of Z2_oy,
    % int_a^b int_y^d Z_oy{i}(eta,x,mu,y)' * Qij * Z_oy{j}(eta,tt,mu,nu) d mu d eta
    %                          ^  ^
    %                          |  Degree of x in constant limit is degree of tt in Z_oy{i}
    %                          Degree of x in Z_oy{i} will not matter in constant limit
    % The variable x vanishes when performing the integration.
    Zd_oy = getSOSmonomials(joint_degs_oy,1);   % Get degrees in [tt,y*nu] for the monomials Z2_oy
    nZd_oy = size(Zd_oy,1);
    % Then, for each joint degree, establish all combinations of degrees in y
    % and nu that add up to this joint degree:
    nZdd_oy = sum(max(Zd_oy(:,2)+1-ydeg_min-nudeg_min,1),1); % Number of monomials in (x,tt,y,nu) with degrees in tt and (y*nu) as in Zd
    Zdd_oy = zeros(nZdd_oy,4);
    indx_n = 0;
    for k=1:nZd_oy
        Zdk_ynu = Zd_oy(k,2);          % Required degree in y*nu
        i1y = min(ydeg_min,Zdk_ynu);
        i1nu = min(nudeg_min,Zdk_ynu-i1y);
        Zddk_ynu = [(i1y:Zdk_ynu-i1nu)',(Zdk_ynu-i1y:-1:i1nu)'];  % Degrees in y,nu that add up to Zdk_ynu
        Zddk = [repmat(Zd_oy(k,1),[Zdk_ynu+1-i1y-i1nu,1]), Zddk_ynu];
        
        % Add to the list of monomials
        indx = indx_n;
        indx_n = indx + size(Zddk,1);
        Zdd_oy(indx+1:indx_n,[2:4]) = Zddk;
    end
    % Then, we have degrees for tt, y and nu. Finally, introduce degrees in x
    % until we reach a desired number of monomials (sufficient d.o.f.):
    nd_x = ceil(fctr_x_Z2*nZdd_xy/nZdd_oy);
    Zdd_oy = repmat(Zdd_oy,[nd_x,1]);
    Zdd_oy(:,1) = kron((0:nd_x-1)',ones(nZdd_oy,1));
    Zdd_oy = makesparse(Zdd_oy);
    
    Z2_oy{2,2} = Zdd_oy;        Z2_oy{2,3} = Zdd_oy;
    Z2_oy{3,2} = Zdd_oy;        Z2_oy{3,3} = Zdd_oy;
    
    % % Next, introduce monomials for separable kernel in y-integral.
    % Establish a set of joint degrees in [x*tt,y*nu]
    joint_degs_xo = joint_degs(joint_degs(:,2)==max(joint_degs(:,2)),:) + [0,1];    % Reintroduce previously substracted degree
    joint_degs_xo = [joint_degs;joint_degs_xo];
    % Note: Using a separable kernel in the y-integral, the joint degree in
    % y*nu of the product Z2_xo{i}*Z2_xo{j} will be completely determined by
    % the degree in just tt of Z2_xo, the y variable of Z2_xo vanishes.
    Zd_xo = getSOSmonomials(joint_degs_xo,1);   % Get degrees in [x*tt,nu] for the monomials Z2_xo
    nZd_xo = size(Zd_xo,1);
    % Then, for each joint degree, establish all combinations of degrees in x
    % and tt that add up to this joint degree:
    nZdd_xo = sum(max(Zd_xo(:,1)+1-xdeg_min-ttdeg_min,1),1); % Number of monomials in (x,tt,y,nu) with degrees in tt and (y*nu) as in Zd
    Zdd_xo = zeros(nZdd_xo,4);
    indx_n = 0;
    for k=1:nZd_xo
        Zdk_xtt = Zd_xo(k,1);          % Required degree in x*tt
        i1x = min(xdeg_min,Zdk_xtt);
        i1tt = min(ttdeg_min,Zdk_xtt-i1x);
        Zddk_xtt = [(i1x:Zdk_xtt-i1tt)',(Zdk_xtt-i1x:-1:i1tt)'];  % Degrees in x,tt that add up to Zdk_xtt
        Zddk = [Zddk_xtt,repmat(Zd_xo(k,2),[Zdk_xtt+1-i1x-i1tt,1])];
        
        % Add to the list of monomials
        indx = indx_n;
        indx_n = indx + size(Zddk,1);
        Zdd_xo(indx+1:indx_n,[1,2,4]) = Zddk;
    end
    % Then, we have degrees for x, tt, and nu. Finally, introduce degrees in y
    % until we reach a desired number of monomials:
    nd_y = ceil(fctr_y_Z2*nZdd_xy/nZdd_xo);
    Zdd_xo = repmat(Zdd_xo,[nd_y,1]);
    Zdd_xo(:,3) = kron((0:nd_y-1)',ones(nZdd_xo,1));
    Zdd_xo = makesparse(Zdd_xo);
    
    Z2_xo{2,2} = Zdd_xo;        Z2_xo{2,3} = Zdd_xo;
    Z2_xo{3,2} = Zdd_xo;        Z2_xo{3,3} = Zdd_xo;
    
    % % Finally, constant limits in both the x- and y-integrals
    % Establish a set of joint degrees in [x*tt,y*nu]
    joint_degs_oo = joint_degs((joint_degs(:,1)==joint_degs(:,1) | joint_degs(:,2)==max(joint_degs(:,2))),:) + [1,1];    % Reintroduce previously substracted degree
    joint_degs_oo = [joint_degs;joint_degs_oo];
    % Note: Separable kernels in both the x- and y-integrals means that x and y
    % vanish. Thus, the joint degree in x*tt and y*nu of the product is
    % entirely determined by the degree in tt and nu of the monomials
    Zd_nutt = getSOSmonomials(joint_degs_oo,1);   % Get degrees in [tt,nu] for the monomials Z2_oo
    nd_nutt = size(Zd_nutt,1);
    fctr = fctr_x_Z2*fctr_y_Z2;
    nZdd_oo = ceil(fctr * nZdd_xy);
    nZ_xy = ceil(nZdd_oo/(nd_nutt));    % How may monomomials in x and y do we need?
    nd_xy = ceil(sqrt(nZ_xy));
    Zdd_x = (0 : nd_xy-1)';
    Zdd_y = (0 : nd_xy-1)';
    Zdd_oo = [repmat(Zdd_x,[nd_xy*nd_nutt,1]), kron(Zd_nutt(:,1),ones(nd_xy^2,1)),...
        repmat(kron(Zdd_y,ones(nd_nutt,1)),[nd_xy,1]), kron(Zd_nutt(:,2),ones(nd_xy^2,1))];
    Zdd_oo = makesparse(Zdd_oo);
    
    Z2_oo{2,2} = Zdd_oo;        Z2_oo{2,3} = Zdd_oo;
    Z2_oo{3,2} = Zdd_oo;        Z2_oo{3,3} = Zdd_oo;
    
else
    Z2_xy{2,2} = zeros(0,4);   Z2_xy{2,3} = zeros(0,4);
    Z2_xy{3,2} = zeros(0,4);   Z2_xy{3,3} = zeros(0,4);
    
    Z2_xo{2,2} = zeros(0,4);   Z2_xo{2,3} = zeros(0,4);
    Z2_xo{3,2} = zeros(0,4);   Z2_xo{3,3} = zeros(0,4);
    
    Z2_oy{2,2} = zeros(0,4);   Z2_oy{2,3} = zeros(0,4);
    Z2_oy{3,2} = zeros(0,4);   Z2_oy{3,3} = zeros(0,4);
    
    Z2_oo{2,2} = zeros(0,4);   Z2_oo{2,3} = zeros(0,4);
    Z2_oo{3,2} = zeros(0,4);   Z2_oo{3,3} = zeros(0,4);
end

degs_xy.Z2 = Z2_xy;     degs_xo.Z2 = Z2_xo;
degs_oy.Z2 = Z2_oy;     degs_oo.Z2 = Z2_oo;

end



function Zd = getSOSmonomials(ZZdegs,use_conv)

% For a given set of monomial degrees ZZdegs, establish a set of
% monomial degrees Zd in the same variables such that a linear combination
% of the monomials defined by ZZdegs could be represented by an SOS
% combination Zd'*Q*Zd of the monomials defined by degrees Zd.

if nargin==1
    use_conv=0;
end

% First, initialize a full set of monomials in the same variables as ZZdegs
% between the minimal and maximal joint degrees in all variables
maxdeg = full(max(sum(ZZdegs,2)));    % maximum joint degree in all variables of the monomials
mindeg = full(min(sum(ZZdegs,2)));    % minimum joint degree in all variables nu of the monomials
Zd = monomials(2,floor(mindeg/2):ceil(maxdeg/2));   % consider of joint degrees between mindeg/2 and maxdeg/2

% Next, we discard some unnecessary monomials
maxdegree = sparse(max(ZZdegs,[],1)/2); % row of max half degrees in each variable separately
mindegree = sparse(min(ZZdegs,[],1)/2); % row of min half degrees in each variable separately
Zdummy1 = bsxfun(@minus,maxdegree,Zd);  % maxdegree monomial minus each monomial
Zdummy2 = bsxfun(@minus,Zd,mindegree);  % each monomial minus mindegree monomial
%[I,~] = find([Zdummy1 Zdummy2]<0);      % rows of monomials that exceed the necessary max degree or are less than the min degree
%IND = setdiff(1:size(Zd,1),I,'stable'); % rows not listed in I
IND = any(Zdummy1>=0,2) | any(Zdummy2>=0,2);          % rows of monomials with max and min degrees within the allowed bounds
Zd = Zd(IND,:);

% Finally, use the Newton polytope to further reduce the set of monomials
% if desired
if use_conv==1
    Zd = inconvhull(full(Zd),full(ZZdegs./2));    % Retain only monomials in convex hull of degs_tot/2 (Newton polytope)
% elseif use_cov==2
%--------------------------------------------------------------------------
% % % Multipartite reduction, not implemented yet...
%     Z2 = ZZdegs./2;           % lots of fractional degrees
%     info2 = sos.expr.multipart{i};  % the vectors of independent variables
%     sizeinfo2m = length(info2);
%     vecindex = [];
%     for indm = 1:sizeinfo2m % for each set of independent variables (first true ind, then matrix)
%         sizeinfo2n(indm) = length(info2{indm}); % number of variables in cell
%         for indn = 1:sizeinfo2n(indm) % scroll through the matrix variables,
%             
%             % PJS 9/12/13: Update code to handle polynomial objects
%             var = info2{indm}(indn);
%             cvartable = char(sos.varmat.vartable);
%             
%             if ispvar(var)
%                 % Convert to string representation
%                 var = var.varname;
%             end
%             varcheckindex = find(strcmp(var,sos.vartable));
%             if ~isempty(varcheckindex)
%                 vecindex{indm}(indn) = varcheckindex;
%             else
%                 vecindex{indm}(indn) = length(info2{1}) + find(strcmp(var,cvartable));
%             end
%             
%             % PJS 9/12/13: Original Code to handle polynomial objects
%             %vecindex{indm}(indn) = find(strcmp(info2{indm}(indn).varname,sos.vartable));
%             
%         end
%     end
%     Zmp = sparsemultipart(full(Zd),full(Z2),vecindex);
%     Zmp = makesparse(Zmp);
%     if ~isempty(Zmp)        % Fix in case result is empty, but might not be appropriate...
%         Zd = Zmp;
%     end
%--------------------------------------------------------------------------
end

Zd = makesparse(Zd);

end