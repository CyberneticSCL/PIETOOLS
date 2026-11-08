function DP = subs(P,invar,inval,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dP = subs(P) evaluates variables invar of (Pu) at values inval
% (P*u)(invar=inval) = DP*[u; u(invar=inval)]
% 
% INPUT 
%   P: opvar object
%   invar: nx1 cell or polynomial
%           - Each element invar{i} or invar(i) should correspond to a
%             variable appearing in the opvar object
%           - In this function, P will be evaluated wrt each of the
%             variables invar(i) and associated values inval(i) separately
%   opts: optional inputs (see OUTPUT for details)
% 
% OUTPUT 
%   DP: opvar operator P so that P*u|_invar=inval = DP*[u;u(invar=inval)]
%           - DP will be such that
%             P*[u0;u1]|_invar=inval = DP*[u0;u1(p);u1]
%           - if opts=='pure', we assume P*u|invar=inval does not depend on
%             u(invar=inval), and DP will be such that
%             P*[u0;u1]|_invar=inval = DP*[u0;u1]
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - subs
%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 10/16/2024

% % Extract the inputs
if nargin>4
    error('At most 4 inputs are accepted')
end
if isempty(invar)
    % No variable to substitute; nothing to do...
    DP = P;
    return
end
if nargin==2 
    if isnumeric(invar)
        inval = invar;
        invar = P.var1;
        opts = '';
    else
        error('A position at which to evaluate the operator must be specified.')
    end
elseif nargin==3
    if isnumeric(invar)
        opts = inval;
        inval = invar;
        invar = P.var1;
    elseif isa(inval,'char')
        error('A position at which to evaluate the operator must be specified.')
    else
        opts = '';
    end
end
% % Error checking
% if ndims(invar)>=3 || ~any(size(invar)==1)
%     error('Input variables must be an nx1 cell of strings, or nx1 polynomial')
% end
if isa(invar,'cell') && ~(length(unique(invar))==length(invar))
    error('Input variable names must be unique')
elseif isa(invar,'polynomial') && ~(length(unique(invar.varname))==length(invar))
    error('Input variable names must be unique')
end
if ~all(size(invar)==size(inval))
    error('The number of input variables must match the number of values at which to evaluate them.')
end


% % Evaluate with respect to each element of invar separately
addvar = invar(:);      addval = inval(:);
% Avoid unnecessary evaluations at e.g. x = x.
rm_vars = isequal(addvar,addval);
addvar = addvar(~rm_vars);
addval = addval(~rm_vars);
if isempty(addvar)
    DP = P;
    return
end
while ~isempty(addvar)
    delvar = addvar(1);         ppp = addval(1);
    addvar = addvar(2:end);     addval = addval(2:end);
    if isa(delvar,'polynomial')
        delvar = delvar.varname;
    end

    % if ~any(ismember(P.var1.varname,delvar))
    %     error('Each input variable must correspond to one of the variables appearing in the opvar object')
    % end
    
    if ispvar(ppp) && all(isequal(delvar,ppp))
        % If we're evaluating at x=x, there's nothing to do...
        DP = P;
        continue
    else
        try ppp = double(ppp);
        catch
            error('Evaluation of the operator at a varying position is not supported')
        end
    end

    % % % Perform the evaluation % % %
    opvar DP;
    DP.I = P.I;     
    DP.var1 = P.var1;   DP.var2 = P.var2;
    dim_new = P.dim;
    if strcmp(P.var1(1).varname,delvar) || strcmp(P.var2(1).varname,delvar)
        % The variable to evalute is the first one
        if strcmp(P.var2(1).varname,delvar)
            warning(['The proposed variable ',delvar,' for substitution is a dummy variable of the opvar.',...
                        ' Substituting the variable ',P.var1(1).varname,' instead.']); 
        end
        sss = P.var1(1);    ttt = P.var2(1);
        if ppp==P.I(1,1)
            eval_up = 0;
        elseif ppp==P.I(1,2)
            eval_up = 1;
        else
            error('Only evaluation at upper and lower boundaries is currently supported')
        end

        % If desired, and possible, do not evaluate the state itself at x=p
        if strcmpi(opts,'exclude') || strcmpi(opts,'pure')
            % Check if the output can be expressed in terms of the original
            % state, not evaluated at x=p
            DR0 = polynomial(subs(P.R.R0,sss,ppp));
            tol = 1e-10;
            if any(any(abs(DR0.C)>tol))
                error(['The desired option ''',opts,''' is not allowed: Delta*P*u is not independent of Delta*u'])
            end
            % Adjust the dimensions, noting that the output does not vary
            % in space
            dim_new(1,1) = P.dim(1,1) + P.dim(2,1);
            dim_new(2,1) = 0;   % The output cannot be a function of x
            DP.dim = dim_new;

            % Build the new map to R^n
            DP.P = [P.P; subs(polynomial(P.Q2),sss,ppp)];
            if eval_up %evaluate at upper boundary
                DP.Q1 = [P.Q1; subs(polynomial(P.R.R1),[sss;ttt],[ppp;sss])];
            else
                DP.Q1 = [P.Q1; subs(polynomial(P.R.R2),[sss;ttt],[ppp;sss])];
            end
        else
            % Adjust the dimensions, noting that our state becomes
            % [u0; u1] --> [u0; u1(p); u1]
            dim_new(1,1) = P.dim(1,1) + P.dim(2,1);
            dim_new(2,1) = 0;   % The output cannot be a function of x
            dim_new(1,2) = P.dim(1,2) +P.dim(2,2);
            dim_new(2,2) = P.dim(2,2);
            DP.dim = dim_new;
    
            % Build the new map to R^n
            DP.P = [P.P, zeros(size(P.P,1),size(P.R.R0,2));
                    polynomial(subs(P.Q2,sss,ppp)), polynomial(subs(P.R.R0,sss,ppp))];
            if eval_up %evaluate at upper boundary
                DP.Q1 = [P.Q1; polynomial(subs(P.R.R1,[sss;ttt],[ppp;sss]))];
            else
                DP.Q1 = [P.Q1; polynomial(subs(P.R.R2,[sss;ttt],[ppp;sss]))];
            end
        end
    else
        % The variable to evalute is some other variable which the operator
        % may or may not depend on
        %warning(['The proposed variable to substitute does not appear in either var1 or var2 of the opvar object.',...
        %            'No substitution is performed'])
        sss = polynomial(delvar);

        DP.dim = dim_new;

        % Build the new map to R^n
        DP.P = subs(polynomial(P.P,sss,ppp));
        DP.Q1 = subs(polynomial(P.Q1),sss,ppp);
        DP.Q2 = subs(polynomial(P.Q1),sss,ppp);
        DP.R.R0 = subs(polynomial(P.R.R0),sss,ppp);
        DP.R.R1 = subs(polynomial(P.R.R2),sss,ppp);
        DP.R.R2 = subs(polynomial(P.R.R1),sss,ppp);
    end
    
    P = DP;

end

end