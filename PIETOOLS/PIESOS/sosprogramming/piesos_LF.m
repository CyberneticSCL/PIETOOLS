function [prog,V,dV,Pmat,Zop] = piesos_LF(prog,PIE,xvars,degmat,opts)
% [PROG,V,DV,PMAT,ZOP] = PIESOS_LF(PROG,PIE,XVARS,DEGMAT,OPTS) declares a
% distributed polynomial Lyapunov functional V(x), in terms of the PIE
% state x, and computes its Lie derivative dV(x) along the PIE
%   d/dt (T*x(t)) = f(x(t)).
% With u=T*x the PDE state, the functional is of one of the forms
%   (PDE states)  V(x) = <Zop*Z(u), Pmat*Zop*Z(u)>,   u=T*x;
%   (PIE states)  V(x) = <Zop*Z(x), Pmat*Zop*Z(u)>,   u=T*x;
% for distributed monomials Z, a basis of tensor-PI operators Zop, and a
% matrix Pmat(s) of polynomials in the primary variables s. For PDE states,
% Pmat is symmetric (or positive semidefinite, with opts.pos), and
%   dV(x) = 2*<Zop*Z(u), Pmat*Zop*(d/dt Z(u))>,   d/dt u = f(x).
% For PIE states, the symmetry condition
%   <Zop*Z(y), Pmat*Zop*Z(T*x)> = <Zop*Z(x), Pmat*Zop*Z(T*y)>
% is enforced for independent y and x, so that dV has the same form.
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program (see `piesos_program`);
% - PIE:    'struct' representing the PIE for which a Lyapunov functional
%           is sought, as output by `convert_PIETOOLS_PDE_nonlinear`, with
%           fields T, f, x, varname_pde, varname_pie, x_tab, vars, and dom;
% - xvars:  n x 1 'polyopvar' array of individual state variables, or
%           cellstr with their names, specifying the state variables on
%           which V depends. These may be PDE state variables, PIE state
%           variables, or both, in which case V is the sum of a functional
%           in the PDE states and one in the PIE states. Can also be
%           specified as 'pde' or 'pie', to use all PDE or all PIE state
%           variables, respectively (PIE.varname_pde, PIE.varname_pie);
% - degmat: m x n array of integers specifying the distributed monomial
%           degrees of the monomials Z in each of the state variables. If
%           specified as m x 1 array, the degrees are taken to be the
%           cumulative degrees of the monomials in all of the variables
%           (then Z may contain more than m rows);
% - opts:   (optional) struct with fields:
%           - 'degs': struct or array specifying the degrees of the
%               monomials to use in the parameterization of the operators
%               Zop. See `piesos_quadvar` for more details;
%           - 'exclude': 3x1 boolean array indicating whether to exclude
%               the multiplier term (exclude(1)=1), lower integral
%               (exclude(2)=1), or upper integral (exclude(3)=1) from the
%               definition of basis operator Zop;
%           - 'psatz': 1xp array of integers specifying whether to use a
%               psatz multiplier to enforce positivity only on the interval
%               over which the operator integrates. If psatz=0, no such
%               multiplier is used. If psatz=1, the matrix Pmat will be
%               multiplied by (s-a)*(b-s), where [a,b]=Z.dom. If
%               psatz=[0,1], we have Pmat = Pmat1 + (s-a)*(b-s)*Pmat2;
%           - 'sep': scalar boolean specifying whether to enforce
%               separability of the operator;
%           - 'pos': scalar boolean, set to true to declare V as a sum of
%               squares (only for PDE state variables). Otherwise, Pmat is
%               symmetric but indefinite;
%           - 'derivative': scalar boolean, set to false to skip computing
%               the derivative dV, which is then returned as [] (default
%               true). It is also skipped if fewer than 3 outputs are
%               requested. Outputs ignored with ~, as in [prog,V,~,Pmat],
%               count as requested, as nargout includes them (MATLAB has
%               no isargout): use this option instead;
%
% OUTPUTS
% - prog:   'struct' specifying the PIESOS program with the decision
%           variables defining Pmat added, and, for PIE state variables,
%           the symmetry condition on Pmat, and, if opts.pos==true, the
%           condition Pmat>=0;
% - V:      scalar 'polyopvar' object representing a distributed polynomial
%           functional in the PIE state variables;
% - dV:     scalar 'polyopvar' object representing the Lie derivative of V
%           along solutions to the PIE, or [] if it was not computed (see
%           opts.derivative). prog, V, Pmat and Zop do not depend on this;
% - Pmat:   m x m cell, with Pmat{i,j} a 'dpvar' object in the primary
%           variables, defining the matrix Pmat(s) (a block for each
%           monomial in Z, see `piesos_quadvar`);
% - Zop:    m x n 'tensopmat' object representing the monomial basis
%           operator acting on the distributed monomial basis defined by
%           Z. If xvars has PDE and PIE state variables, Pmat and Zop are
%           2 x 1 cells with the output for the PDE and for the PIE
%           states;
%
% NOTES
% - The Lie derivative composes the functional with the operators defining
%   the PIE, also for nonlinear terms in f, the pointwise products of
%   operators acting on the state, which give nested tensor-PI operators.
% - The derivative is computed in the linear representation of V, as a
%   polynomial in x with a functional as coefficient. `piesos_LF_alt`
%   computes it in the quadratic representation, which does not support
%   nonlinear terms in f, unless Zop only multiplies.
% - Based on the function `V3_Liediff` by CR (see Examples/LocalStability).
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_LF
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
% DJ, 10/15/2026: Initial coding, based on V3_Liediff by CR;
% DJ, 10/23/2026: Fix for vector-valued states;

% Check the options
if nargin<=4
    opts = struct();
end
if ~isa(opts,'struct')
    error("Options must be specified as a 'struct'.")
end
% Skip the derivative if it is not requested (also for the recursive calls below)
if isfield(opts,'derivative') && ~(isscalar(opts.derivative) && (islogical(opts.derivative) || isnumeric(opts.derivative)))
    error("Option 'derivative' must be specified as a scalar boolean.")
end
if nargout<3
    opts.derivative = false;
end
get_dV = ~isfield(opts,'derivative') || logical(opts.derivative);

% Extract the PDE and PIE state variables
varname_pde = PIE.varname_pde;
varname_pie = PIE.varname_pie;

% Check that the state variables are properly declared
if (ischar(xvars) || (isstring(xvars) && isscalar(xvars))) && ~isempty(xvars)
    % Use all PDE or all PIE state variables
    switch lower(char(xvars))
        case 'pde'
            xvars = varname_pde;
        case 'pie'
            xvars = varname_pie;
        otherwise
            error("State variables specified as char must be 'pde' or 'pie'.")
    end
end
if isempty(xvars)
    error("At least one state variable must be specified.")
elseif iscellstr(xvars)
    % Only names of state variables are provided
    varnames = xvars(:);
elseif isa(xvars,'polyopvar')
    % State variables are already provided as polyopvar objects
    [is_pvar_xvars,varnames,dims_x,pvars_x,doms_x] = ispvar(xvars);
    if ~is_pvar_xvars
        error("State variables must be specified as 'cellstr' or 'polyopvar' array of individual variables.")
    end
else
    error("State variables must be specified as 'cellstr' or 'polyopvar' array of individual variables.")
end
% Make sure the state variables appear in the PDE or PIE
[~,idcs_pde,xidcs_pde] = intersect(varname_pde,varnames);
[~,idcs_pie,xidcs_pie] = intersect(varname_pie,varnames);
if numel(xidcs_pde)+numel(xidcs_pie)~=numel(unique(varnames))
    error("All variables must appear as state in either the PIE or PDE representation.")
end

if ~isempty(idcs_pde) && ~isempty(idcs_pie)
    % Take V as the sum of two functionals
    if size(degmat,2)==numel(varnames) && numel(varnames)>1
        % Monomials in the variables of each functional
        degmat_pde = degmat(:,xidcs_pde);
        degmat_pde = degmat_pde(any(degmat_pde,2),:);
        degmat_pie = degmat(:,xidcs_pie);
        degmat_pie = degmat_pie(any(degmat_pie,2),:);
    else
        degmat_pde = degmat;
        degmat_pie = degmat;
    end
    if isempty(degmat_pde) || isempty(degmat_pie)
        % All monomials are in the variables of one of the functionals
        if isempty(degmat_pde)
            [prog,V,dV,Pmat,Zop] = piesos_LF(prog,PIE,varnames(xidcs_pie),degmat_pie,opts);
        else
            [prog,V,dV,Pmat,Zop] = piesos_LF(prog,PIE,varnames(xidcs_pde),degmat_pde,opts);
        end
        return
    end
    [prog,V1,dV1,Pmat1,Zop1] = piesos_LF(prog,PIE,varnames(xidcs_pde),degmat_pde,opts);
    [prog,V2,dV2,Pmat2,Zop2] = piesos_LF(prog,PIE,varnames(xidcs_pie),degmat_pie,opts);
    V = V1+V2;
    dV = dV1+dV2;                       % [] if the derivative is not computed (opts is passed on)
    Pmat = {Pmat1;Pmat2};
    Zop = {Zop1;Zop2};
    return
elseif ~isempty(idcs_pde)
    use_pde = true;
    idcs_PIE = idcs_pde;
    xidcs = xidcs_pde;
else
    use_pde = false;
    idcs_PIE = idcs_pie;
    xidcs = xidcs_pie;
end
if size(degmat,2)==numel(varnames) && numel(varnames)>1
    degmat = degmat(:,xidcs);
end

% Determine the function space in which the states exist
pvars = pvar2varname(PIE.vars(:,1))';                       % spatial variables
xsize = PIE.x_tab(idcs_PIE,2);                              % state size
xdepmat = logical(PIE.x_tab(idcs_PIE,2+(1:numel(pvars))));  % dependence on spatial variables
dom = PIE.dom;                                              % spatial domain
if isa(xvars,'polyopvar')
    % Check that the provided state variables exist in the same domain as
    % their counterparts in the PIE
    [~,idcs_x] = ismember(varnames(xidcs),varnames);
    for i=1:numel(xidcs)
        deps_i = xdepmat(i,:);
        if ~isequal(dims_x(idcs_x(i)),xsize(i)) ||...
              ~isequal(pvars_x{idcs_x(i)},pvars(deps_i)) ||...
                ~isequal(doms_x{idcs_x(i)},dom(deps_i,:))
            error("Function spaces of state variables in PIE do not match those of provided variables.")
        end
    end
end

% Extract the map from PIE to PDE state, and the right-hand side of the PIE
Top = PIE.T;
fx = PIE.f;
if ~isequal(idcs_PIE(:)',1:numel(varname_pie))
    % Indexing is not supported for operators of multiple components, so
    % only index if V depends on a subset of the states
    % A vector-valued state occupies several rows of T and f
    xrow0 = cumsum([0;PIE.x_tab(:,2)]);
    rows_PIE = cell2mat(arrayfun(@(i) xrow0(i)+1:xrow0(i+1),idcs_PIE(:)','UniformOutput',false));
    Top = Top(rows_PIE,:);
    fx = fx(rows_PIE,:);                    % d/dt (Tx) = f(x)
end
Tx = Top*PIE.x;                             % PDE state in terms of PIE state
names_pde = varname_pde(idcs_PIE);          % PDE states, in the order of the rows of Tx, fx
names_x = varname_pie(idcs_PIE);            % associated PIE states

% Create a basis of monomials of degree degmat in the PDE state variables,
% and the same monomials in the PIE state variables
xu = polyopvar(names_pde,xsize,pvars,dom,xdepmat);
Zu = dmonomials(xu,degmat);
if any(sum(Zu.Z>0,2)==0)
    error("Constant (degree zero) monomials are not supported.")
end
if use_pde
    % The left monomials are composed with T below: use different names
    names_L = strcat(names_x,'_L');
else
    names_L = names_x;
end
[~,ord] = ismember(Zu.varname,names_pde);
Zx = subs(Zu,Zu.varname,names_L(ord));

% Declare the functional
%   Vxu(x,u) = <Zop Z(x), Pmat Zop Z(u)>
% in terms of the monomials Z(x) in PIE states x and Z(u) in PDE states u.
% (piesos_quadvar acts on the distributed monomials Z, not on the operators
% defining them, so T is composed with the left monomials below.)
if use_pde
    % Symmetric Pmat is required to compute the derivative of V(u,u)
    if ~(isfield(opts,'pos') && opts.pos)
        opts.sym = true;
    end
elseif isfield(opts,'pos') && opts.pos
    % This would require declaring another functional W(x)>=0, and
    % enforcing V(x)=W(x)
    error("Enforcing positivity of mixed PDE-PIE state LF is currently not supported.")
end
[prog,Vxu,Pmat,ZopL,ZopR] = piesos_quadvar(prog,Zx,Zu,opts);                % <-- is there any advantage using piesos_polyvar instead for non-sym case?
if ~isequal(ZopL,ZopR)
    error("The same operators must be used on the left and right of the Lyapunov functional.")
end
Zop = ZopL;

% Compute V(x) and its derivative
%   dV(x) = 2*<Zop Z(.), Pmat Zop d/dt Z(u)>
% where d/dt u = f(x), so that the derivative acts only on the right
% monomials, Z(u), and the symmetry of V accounts for the left monomials.
% Take the derivative of the (linear) representation V = K*Z
if use_pde
    % Compose the left monomials with T: x_L -> T*x
    Wxu = subs(Vxu,names_L,Tx);
else
    Wxu = Vxu;
end
% Express V in terms of solely the PIE state x
V = subs(Wxu,names_pde,Tx);
% Compute the derivative: in each monomial, replace each factor u in turn
% by f(x), and the other factors u by T*x, based on the product rule
% The derivative is computed only if requested (get_dV); the unused states of V are removed after
if get_dV
    dV = 0;
    for k=find(~cellfun('isempty',Wxu.C.ops(1,:)))
        Wk = Wxu{1,k};
        zk = full(Wk.Z(1,:));
        [is_u,rows] = ismember(Wk.varname(zk(zk>0)),names_pde);     % factors u, and their rows in T, f
        fctrs = find(is_u(:))';                                     % row, for the loops
        for j=fctrs
            % Rows of f of the state of factor j
            xrowsel = cumsum([0;xsize(:)]);
            if size(fx,1)==xsize(rows(j))
                fj = fx;                % indexing not supported for multiple components
            else
                fj = fx(xrowsel(rows(j))+1:xrowsel(rows(j)+1),:);
            end
            dV = dV + subs(subs(Wk,j,fj),names_pde,Tx);
        end
    end
    dV = 2*dV;
    dV = remove_unused_states(dV);
else
    dV = [];
end
V = remove_unused_states(V);

if ~use_pde
    % Enforce the symmetry condition V(xL,T*xR) = V(xR,T*xL) for independent
    % copies xL and xR of the PIE state
    TxL = subs(Tx,Tx.varname,strcat(Tx.varname,'_L'));
    TxR = subs(Tx,Tx.varname,strcat(Tx.varname,'_R'));
    WLR = subs(subs(Vxu,names_x,strcat(names_x,'_L')),names_pde,TxR);
    WRL = subs(subs(Vxu,names_x,strcat(names_x,'_R')),names_pde,TxL);
    prog = piesos_eq(prog,WLR-WRL);
    if isfield(opts,'pos') && opts.pos
        % Enforce positivity using piesos_ineq
        prog = piesos_ineq(prog,V);
    end
end

end