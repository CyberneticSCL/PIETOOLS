function [DP] = delta_opvar2d(P,invar,inval,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [[dP] = delta_opvar(P) evaluates variables invar of (Pu) at values inval
% (P*u)(invar=inval) = DP*[u; u(invar=inval)]
% 
% INPUT 
%   P: opvar2d object
%   invar: nx1 cell or polynomial
%           - Each element invar{i} or invar(i) should correspond to one of
%             the variables P.var1 appearing in the opvar2d object
%           - In this function, P will be evaluated wrt each of the
%             variables invar(i) and associated values inval(i) separately
%   opts: optional inputs (see OUTPUT for details)
% 
% OUTPUT 
%   DP: opvar2d operator P so that P*u|_invar=inval = DP*[u;u(invar=inval)]
%           - DP will be such that
%             P*[u0;ux;uy;u2] = DP*[u0;ux(p);ux;uy;u2(p,y);u2]
%             if invar==ss1, or
%             P*[u0;ux;uy;u2] = DP*[u0;uy(p);ux;u2(x,p);uy;u2]
%             if invar==ss2
%           - if opts=='pure', we assume P*u|invar=inval does not depend on
%             u(invar=inval), and DP will be such that
%             P*[u0;ux;uy;u2] = DP*[u0;ux;uy;u2]
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - delta_opvar2d
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
% Initial coding DJ - 08/27/2021
% DJ, 06/17/2022 - Support edge boundary inputs (e.g. inval=(x,1)).

% % Error checking
if ~isa(P,'opvar2d')
    error("Input function must be an opvar2d class object");
end

if ndims(invar)>=3 || ~any(size(invar)==1)
    error('Input variables must be an nx1 cell of strings, or nx1 polynomial')
end
if isa(invar,'cell') && ~(length(unique(invar))==length(invar))
    error('Input variable names must be unique')
elseif isa(invar,'polynomial') && ~(length(unique(invar.varname))==length(invar))
    error('Input variable names must be unique')
end
if ~all(size(invar)==size(inval))
    error('The number of input variables must match the number of values at which to evaluate them')
end


% % Evaluate with respect to each element of invar separately
addvar = invar(:);      addval = inval(:);
% Avoid unnecessary evaluations at e.g. x = x.
rm_vars = isequal(addvar,addval);
addvar = addvar(~rm_vars);
addval = addval(~rm_vars);
while ~isempty(addvar)
    delvar = addvar(1);         ppp = addval(1);
    addvar = addvar(2:end);     addval = addval(2:end);
    if isa(delvar,'polynomial')
        delvar = delvar.varname;
    end

    if ~any(ismember(P.var1.varname,delvar))
        error('Each input variable must correspond to one of the variables appearing in the opvar2d object')
    end
    
    if ispvar(ppp) && all(isequal(delvar,ppp))
        DP = P;
        continue
    else
        try ppp = double(ppp);
        catch
            error('Evaluation of the operator at a varying position is not supported')
        end
    end

    % % % Perform the evaluation % % %
    opvar2d DP;
    DP.I = P.I;     
    DP.var1 = P.var1;   DP.var2 = P.var2;
    dim_new = P.dim;
    if strcmp(P.var1(1).varname,delvar) || strcmp(P.var2(1).varname,delvar)
        % The variable to evalute is the first one
        if strcmp(P.var2(1).varname,delvar)
            warning(['The proposed variable ',delvar,' for substitution is a dummy variable of the opvar2d.',...
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

        % Adjust the dimensions, noting that our state becomes
        % [u0; ux; uy; u2] --> [u0; ux(p); ux; uy; u2(p,x); u2]
        dim_new(1,1) = P.dim(1,1) + P.dim(2,1);
        dim_new(2,1) = 0;   % The output cannot be a function in x
        dim_new(3,1) = P.dim(3,1) + P.dim(4,1);
        dim_new(4,1) = 0;   % The output cannot be a function in x,y
        dim_new(1,2) = P.dim(1,2) + P.dim(2,2); % state [ux] becomes [ux(p); ux];
        dim_new(2,2) = P.dim(2,2);
        dim_new(3,2) = P.dim(3,2) + P.dim(4,2); % state [u2] becomes [u2(p,y); u2];
        dim_new(4,2) = P.dim(4,2);
        DP.dim = dim_new;

        % Build the new map to R^n
        R00 = [P.R00; polynomial(subs(P.Rx0,sss,ppp))];
        DP.R00 = [R00, [zeros(size(P.R00,1),size(P.Rxx{1},2)); 
                        polynomial(subs(P.Rxx{1},sss,ppp))]];
        if eval_up %evaluate at upper boundary
            DP.R0x = [P.R0x; polynomial(subs(P.Rxx{2},[sss;ttt],[ppp;sss]))];
        else
            DP.R0x = [P.R0x; polynomial(subs(P.Rxx{3},[sss;ttt],[ppp;sss]))];
        end
        R0y = [P.R0y; polynomial(subs(P.Rxy,sss,ppp))];
        DP.R0y = [R0y, [zeros(size(P.R0y,1),size(P.Rx2{1},2)); 
                        polynomial(subs(P.Rx2{1},sss,ppp))]];
        if eval_up
            DP.R02 = [P.R02; polynomial(subs(P.Rx2{2},[sss;ttt],[ppp;sss]))];
        else
            DP.R02 = [P.R02; polynomial(subs(P.Rx2{3},[sss;ttt],[ppp;sss]))];
        end

        % Build the new map to L2[y]
        Ry0 = [P.Ry0; polynomial(subs(P.R20,sss,ppp))];
        DP.Ry0 = [Ry0, [zeros(size(P.Ry0,1),size(P.R2x{1},2)); 
                        polynomial(subs(P.R2x{1},sss,ppp))]];
        if eval_up
            DP.Ryx = [P.Ryx; polynomial(subs(P.R2x{2},[sss;ttt],[ppp;sss]))];
        else
            DP.Ryx = [P.Ryx; polynomial(subs(P.R2x{3},[sss;ttt],[ppp;sss]))];
        end
        Ryy{1} = [P.Ryy{1}; polynomial(subs(P.R2y{1},sss,ppp))];
        DP.Ryy{1} = [Ryy{1}, [zeros(size(P.Ryy{1},1),size(P.R22{1,1},2)); 
                                polynomial(subs(P.R22{1,1},sss,ppp))]];
        Ryy{2} = [P.Ryy{2}; polynomial(subs(P.R2y{2},sss,ppp))];
        DP.Ryy{2} = [Ryy{2}, [zeros(size(P.Ryy{2},1),size(P.R22{1,2},2)); 
                                polynomial(subs(P.R22{1,2},sss,ppp))]];
        Ryy{3} = [P.Ryy{3}; polynomial(subs(P.R2y{3},sss,ppp))];
        DP.Ryy{3} = [Ryy{3}, [zeros(size(P.Ryy{3},1),size(P.R22{1,3},2)); 
                                polynomial(subs(P.R22{1,3},sss,ppp))]];
        if eval_up
            DP.Ry2{1} = [P.Ry2{1}; polynomial(subs(P.R22{2,1},[sss;ttt],[ppp;sss]))];
            DP.Ry2{2} = [P.Ry2{2}; polynomial(subs(P.R22{2,2},[sss;ttt],[ppp;sss]))];
            DP.Ry2{3} = [P.Ry2{3}; polynomial(subs(P.R22{2,3},[sss;ttt],[ppp;sss]))];
        else
            DP.Ry2{1} = [P.Ry2{1}; polynomial(subs(P.R22{3,1},[sss;ttt],[ppp;sss]))];
            DP.Ry2{2} = [P.Ry2{2}; polynomial(subs(P.R22{3,2},[sss;ttt],[ppp;sss]))];
            DP.Ry2{3} = [P.Ry2{3}; polynomial(subs(P.R22{3,3},[sss;ttt],[ppp;sss]))];
        end

        % If desired, and possible, get rid of ux(p) and u2(p,y) contributions
        if nargin==4 && (strcmpi(opts,'exclude') || strcmpi(opts,'pure'))
            exclude_c = [P.dim(1,2)+1:P.dim(1,2)+P.dim(2,2), P.dim(1,2)+2*P.dim(2,2)+P.dim(3,2)+1: P.dim(1,2)+2*P.dim(2,2)+P.dim(3,2)+P.dim(4,2)];
            if ~isempty(exclude_c) && DP(:,exclude_c)==0
                include_c = 1:size(DP,2);
                include_c(exclude_c) = [];
                DP = DP(:,include_c);
            elseif ~isempty(exclude_c)
                error(['The desired option ''',opts,''' is not allowed: Delta*P*u is not independent of Delta*u'])
            end
        end

    elseif strcmp(P.var1(2).varname,delvar) || strcmp(P.var2(2).varname,delvar)
        % The variable to evalute is the second one
        if strcmp(P.var2(1).varname,delvar)
            warning(['The proposed variable ',delvar,' for substitution is a dummy variable of the opvar2d.',...
                        ' Substituting the variable ',P.var1(2).varname,' instead.']); 
        end
        sss = P.var1(2);    ttt = P.var2(2);
        if ppp==P.I(2,1)
            eval_up = 0;
        elseif ppp==P.I(2,2)
            eval_up = 1;
        else
            error('Only evaluation at upper and lower boundaries is currently supported')
        end

        % Adjust the dimensions, noting that our state becomes
        % [u0; ux; uy; u2] --> [u0; uy(p); ux; u2(x,p); uy; u2]
        dim_new(1,1) = P.dim(1,1) + P.dim(3,1);
        dim_new(2,1) = P.dim(2,1) + P.dim(4,1);
        dim_new(3,1) = 0;   % The output cannot be a function in y
        dim_new(4,1) = 0;   % The output cannot be a function in x,y
        dim_new(1,2) = P.dim(1,2) + P.dim(3,2); % state [uy] becomes [uy(p); uy];
        dim_new(2,2) = P.dim(2,2) + P.dim(4,2); % state [u2] becomes [u2(x,p); u2];
        dim_new(3,2) = P.dim(3,2);
        dim_new(4,2) = P.dim(4,2);
        DP.dim = dim_new;

        % Build the new map to R^n
        R00 = [P.R00; polynomial(subs(P.Ry0,sss,ppp))];
        DP.R00 = [R00, [zeros(size(P.R00,1),size(P.Ryy{1},2)); 
                        polynomial(subs(P.Ryy{1},sss,ppp))]];
        R0x = [P.R0x; polynomial(subs(P.Ryx,sss,ppp))];
        DP.R0x = [R0x, [zeros(size(P.R0x,1),size(P.Ry2{1},2)); 
                        polynomial(subs(P.Ry2{1},sss,ppp))]];
        if eval_up %evaluate at upper boundary
            DP.R0y = [P.R0y; polynomial(subs(P.Ryy{2},[sss;ttt],[ppp;sss]))];
        else
            DP.R0y = [P.R0y; polynomial(subs(P.Ryy{3},[sss;ttt],[ppp;sss]))];
        end
        if eval_up
            DP.R02 = [P.R02; polynomial(subs(P.Ry2{2},[sss;ttt],[ppp;sss]))];
        else
            DP.R02 = [P.R02; polynomial(subs(P.Ry2{3},[sss;ttt],[ppp;sss]))];
        end

        % Build the new map to L2[x]
        Rx0 = [P.Rx0; polynomial(subs(P.R20,sss,ppp))];
        DP.Rx0 = [Rx0, [zeros(size(P.Rx0,1),size(P.R2y{1},2)); 
                        polynomial(subs(P.R2y{1},sss,ppp))]];
        Rxx{1} = [P.Rxx{1}; polynomial(subs(P.R2x{1},sss,ppp))];
        DP.Rxx{1} = [Rxx{1}, [zeros(size(P.Rxx{1},1),size(P.R22{1,1},2)); 
                                polynomial(subs(P.R22{1,1},sss,ppp))]];
        Rxx{2} = [P.Rxx{2}; polynomial(subs(P.R2x{2},sss,ppp))];
        DP.Rxx{2} = [Rxx{2}, [zeros(size(P.Rxx{2},1),size(P.R22{2,1},2)); 
                                polynomial(subs(P.R22{2,1},sss,ppp))]];
        Rxx{3} = [P.Rxx{3}; polynomial(subs(P.R2x{3},sss,ppp))];
        DP.Rxx{3} = [Rxx{3}, [zeros(size(P.Rxx{3},1),size(P.R22{3,1},2)); 
                                polynomial(subs(P.R22{3,1},sss,ppp))]];
        if eval_up
            DP.Rxy = [P.Rxy; polynomial(subs(P.R2y{2},[sss;ttt],[ppp;sss]))];
        else
            DP.Rxy = [P.Rxy; polynomial(subs(P.R2y{3},[sss;ttt],[ppp;sss]))];
        end
        if eval_up
            DP.Rx2{1} = [P.Rx2{1}; polynomial(subs(P.R22{1,2},[sss;ttt],[ppp;sss]))];
            DP.Rx2{2} = [P.Rx2{2}; polynomial(subs(P.R22{2,2},[sss;ttt],[ppp;sss]))];
            DP.Rx2{3} = [P.Rx2{3}; polynomial(subs(P.R22{3,2},[sss;ttt],[ppp;sss]))];
        else
            DP.Rx2{1} = [P.Rx2{1}; polynomial(subs(P.R22{1,3},[sss;ttt],[ppp;sss]))];
            DP.Rx2{2} = [P.Rx2{2}; polynomial(subs(P.R22{2,3},[sss;ttt],[ppp;sss]))];
            DP.Rx2{3} = [P.Rx2{3}; polynomial(subs(P.R22{3,3},[sss;ttt],[ppp;sss]))];
        end

        % If desired, and possible, get rid of ux(p) and u2(p,y) contributions
        if nargin==4 && (strcmpi(opts,'exclude') || strcmpi(opts,'pure'))
            exclude_c = [P.dim(1,2)+1:P.dim(1,2)+P.dim(3,2), P.dim(1,2)+2*P.dim(3,2)+P.dim(2,2)+1: P.dim(1,2)+2*P.dim(3,2)+P.dim(2,2)+P.dim(4,2)];
            if ~isempty(exclude_c) && DP(:,exclude_c)==0
                include_c = 1:size(DP,2);
                include_c(exclude_c) = [];
                DP = DP(:,include_c);
            elseif ~isempty(exclude_c)
                error(['The desired option ''',opts,''' is not allowed: Delta*P*u is not independent of Delta*u'])
            end
        end
        
    else
        warning(['The proposed variable to substitute does not appear in either var1 or var2 of the opvar2d object.',...
                    'No substitution is performed'])
    end
    
    P = DP;

end

% % Continue evaluation with respect to remaining variables
% if ~isempty(addvar) && nargin==4
%     DP = delta_opvar2d(DP,addvar,addval,opts);
% elseif ~isempty(addvar)
%     DP = diff_opvar2d(DP,addvar,addval);
% end

end