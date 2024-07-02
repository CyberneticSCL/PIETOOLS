function Eint = int_simple(E,ivar,x1,x2)
% Eint = int(E,var,L,U) integrates the dpvar or polynomial E in independent 
% variable var from lower limit L to upper limit U
%
% NOTE: This function is distinct from the standard "int" functions in that
% it allows only (scalar) 'pvar' or 'double' expressions for the upper and
% lower limit. More complicated expressions will produce an error!
% 
% INPUTS:
% E: A dpvar or polynomial class object
% ivar: A scalar pvar variable (1x1 polynomial object)
% x1: Either a scalar pvar variable or scalar value
%       - If a pvar variable, it cannot match input "var"
%       - This input provides the lower limit of the definite integral
% x2: Either a scalar pvar variable or scalar value
%       - If a pvar variable, it cannot match input "var"
%       - This input provides the upper limit of the definite integral
% 
% OUTPUTS:
% Eint: integrated dpvar object: Eint = int_x1^x2 (E) d dx
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ, MP, SS - 04/14/2022
% DJ, 05/23/22 - Allow inputs of class 'double'

if nargin<4
    error('At least 4 input arguments are required')
elseif nargin<5
    tol = 1e-15;    % Coefficients below tolerance are set to zero
end

% Conver inputs of class 'double' to 'polynomial'
if isa(E,'double')
    E = polynomial(E);
end

% Correction for zero dimension case
if numel(E.C)==0
    if isa(E,'polynomial')
        Eint = polynomial(zeros(size(E)));
    elseif isa(E,'dpvar')
        Eint = dpvar(zeros(size(E)));
    end
    return
end

if ispvar(ivar)
    ivar = ivar.varname{1};
elseif ~ischar(ivar)
    error('Integration variable must be of type "polynomial" or "char"');
end

% Find index of variable we are integrating with respect to
ivar_idx = strcmp(ivar,E.varname);

% Separately deal with case where the integration variable does not appear
% in the integrand
if ~any(ivar_idx)
    Eint = E*(x2-x1);
    return
end

% Distinguish case where input is dpvar and polynomial
if isa(E,'dpvar')
    % Integrate the monomial terms specified by E.degmat
    Z_deg = E.degmat;
    nz = size(Z_deg,1);
    Z_var = E.varname;
    nr_Z = nz;  nc_Z = 1;
    Z_coef = speye(nz);
elseif isa(E,'polynomial')
    % Otherwise, just perform direct integration
    Z_deg = E.degmat;
    nz = size(Z_deg,1);
    Z_var = E.varname;
    [nr_Z,nc_Z] = size(E);
    Z_coef = E.coefficient;
else
    error('Input must be of type "polynomial" or "dpvar"')
end

% Perform the integration
Zint_coef = Z_coef./(Z_deg(:,ivar_idx)+1);
Zint_deg = Z_deg;
Zint_deg(:,ivar_idx) = Z_deg(:,ivar_idx)+1;
Zint_var = Z_var;

% Substitue the limits, distinguishing 4 cases:
if isa(x1,'double') && isa(x2,'double')
    Zint_coef = Zint_coef.*(x2.^Zint_deg(:,ivar_idx)) - Zint_coef.*(x1.^Zint_deg(:,ivar_idx));
    Zint_deg(:,ivar_idx) = [];
    Zint_var(ivar_idx) = [];
    
elseif isa(x1,'double') && ispvar(x2)
    % Check if the new variable already appears in the expression
    var2 = x2.varname{1};
    if strcmp(var2,ivar)
        error('The limits of the integral cannot depend on the integration variable!')
    end
    var2_idx = strcmp(var2,Zint_var);
    if any(var2_idx)
        % If the variable already appears, combine contributions from
        % old and new variable, and remove the old variable
        Zint_coef = [Zint_coef;
                    -Zint_coef.*(x1.^Zint_deg(:,ivar_idx))];    % Subtract expression at x1 from that at x2
        Zint_var = Zint_var(~ivar_idx);                         % Remove the old variable
        Zint_deg2 = Zint_deg;
        Zint_deg2(:,var2_idx) = Zint_deg(:,var2_idx) + Zint_deg(:,ivar_idx);
        Zint_deg = [Zint_deg2(:,~ivar_idx);
            Zint_deg(:,~ivar_idx)];
    else
        % If the new variable does not appear yet, just replace the
        % old variable
        Zint_coef = [Zint_coef;
                    -Zint_coef.*(x1.^Zint_deg(:,ivar_idx))];
        Zint_var{ivar_idx} = var2;                              % Replace the old variable
        Zint_deg1 = Zint_deg;
        Zint_deg1(:,ivar_idx) = 0;
        Zint_deg = [Zint_deg; Zint_deg1];
    end
    
elseif ispvar(x1) && isa(x2,'double')
    % Check if the new variable already appears in the expression
    var1 = x1.varname{1};
    if strcmp(var1,ivar)
        error('The limits of the integral cannot depend on the integration variable!')
    end
    var1_idx = strcmp(var1,Zint_var);
    if any(var1_idx)
        % If the variable already appears, combine contributions from
        % old and new variable, and remove the old variable
        Zint_coef = [Zint_coef.*(x2.^Zint_deg(:,ivar_idx));
                    -Zint_coef];                            % Subtract expression at x1 from that at x2
        Zint_deg1 = Zint_deg;
        Zint_deg1(:,var1_idx) = Zint_deg(:,var1_idx) + Zint_deg(:,ivar_idx);
        Zint_deg = [Zint_deg(:,~ivar_idx);
                    Zint_deg1(:,~ivar_idx)];
        Zint_var = Zint_var(~ivar_idx);                     % Remove the old variable
    else
        % If the new variable does not appear yet, just replace the
        % old variable
        Zint_coef = [Zint_coef.*(x2.^Zint_deg(:,ivar_idx));
            -Zint_coef];
        Zint_var{ivar_idx} = var1;          % Replace the old variable
        Zint_deg2 = Zint_deg;
        Zint_deg2(:,ivar_idx) = 0;
        Zint_deg = [Zint_deg2; Zint_deg];
    end
    
elseif ispvar(x1) && ispvar(x2)
    % Check if the new variables already appears in the expression
    var1 = x1.varname{1};
    var1_idx = strcmp(var1,Zint_var);
    var2 = x2.varname{1};
    var2_idx = strcmp(var2,Zint_var);
    if strcmp(var1,ivar) || strcmp(var2,ivar)
        error('The limits of the integral cannot depend on the integration variable!')
    end
    if any(var1_idx) && any(var2_idx)
        % If both variables appear, combine contributions from old and
        % new variable, and remove the old variable
        Zint_coef = [Zint_coef; -Zint_coef];    % Subtract expression at x1 from that at x2
        Zint_var = Zint_var(~ivar_idx);         % Remove the old variable
        Zint_deg1 = Zint_deg;
        Zint_deg1(:,var1_idx) = Zint_deg(:,var1_idx) + Zint_deg(:,ivar_idx);
        Zint_deg2 = Zint_deg;
        Zint_deg2(:,var2_idx) = Zint_deg(:,var2_idx) + Zint_deg(:,ivar_idx);
        Zint_deg = [Zint_deg2(:,~ivar_idx); Zint_deg1(:,~ivar_idx)];
    elseif any(var2_idx)
        % Otherwise, combine contributions from old variable and x2,
        % and replace the old variable with x1
        Zint_coef = [Zint_coef; -Zint_coef];
        Zint_var{ivar_idx} = var1;              % Replace the old variable
        Zint_deg2 = Zint_deg;
        Zint_deg2(:,var2_idx) = Zint_deg(:,var2_idx) + Zint_deg(:,ivar_idx);
        Zint_deg2(:,ivar_idx) = 0;
        Zint_deg = [Zint_deg2; Zint_deg];
    elseif any(var1_idx)
        % Or, combine contributions from old variable and x1,
        % and replace the old variable with x2
        Zint_coef = [Zint_coef; -Zint_coef];
        Zint_var{ivar_idx} = var2;              % Replace the old variable
        Zint_deg1 = Zint_deg;
        Zint_deg1(:,var1_idx) = Zint_deg(:,var1_idx) + Zint_deg(:,ivar_idx);
        Zint_deg1(:,ivar_idx) = 0;
        Zint_deg = [Zint_deg; Zint_deg1];
    else
        % If neither of the new variable appears yet, just replace
        % the old variable by x2
        Zint_coef = [Zint_coef; -Zint_coef];
        Zint_var{ivar_idx} = var2;              % Replace the old variable
        Zint_var = [Zint_var;{var1}];           % Add the other variable
        Zint_deg1 = Zint_deg;
        Zint_deg1(:,ivar_idx) = 0;
        Zint_deg1 = [Zint_deg1,Zint_deg(:,ivar_idx)];
        Zint_deg = [[Zint_deg,zeros(nz,1)]; Zint_deg1];
    end
    
else
    error('Upper and lower limits must be of type "double" or "pvar"; try calling the standard function "int" instead')
end

% Construct the new polynomial, and combine the different terms
Eint = combine(polynomial(Zint_coef,Zint_deg,Zint_var,[nr_Z nc_Z]));


% If our input is dpvar, adjust the coefficients to account for the
% decision variables: multiply each row/col block of E.C with the 
% integration matrix Eint.coeff
if isa(E,'dpvar')
    Cn = E.C*kron(speye(E.matdim(2)),Eint.coeff');

    % Set numerically insignificant coefficients to zero
    [i,j,vals] = find(Cn);
    idx = (abs(vals)>=tol);
    Cn = sparse(i(idx),j(idx),vals(idx),size(Cn,1),size(Cn,2));

    Eint = dpvar(Cn, Eint.degmat, Eint.varname, E.dvarname, E.matdim);
end

end