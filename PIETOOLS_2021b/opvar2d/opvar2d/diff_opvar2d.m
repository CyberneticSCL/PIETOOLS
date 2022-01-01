function [dP] = diff_opvar2d(P,invar,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [[dP] = diff_opvar(P) calculates spatial derivative of (Pu)(s) with respect to
% variable var so that d(P*u) = dP*[u;du]
% 
% INPUT 
%   P: opvar2d object
%   invar: nx1 cell or polynomial
%           - Each element invar{i} or invar(i) should correspond to one of
%             the variables P.var1 appearing in the opvar2d object
%           - In this function, P will be differentiated wrt to invar(1)
%             first, the result of which will be differentiated wrt
%             invar(2), etc. until no variables are left
%           - If invar(i) = x^p for some power p, differentiation wrt x
%             will be performed p times
%           - If invar(i) = x*y, invar(i) will be split into two separate
%             variables [x;y]
%   opts: optional inputs (see OUTPUT for details)
% 
% OUTPUT 
%   dP: derivative of P with respect to variable invar
%           - dP will be such that
%             (d/d invar)*P*[u0;ux;uy;u2] = dP*[u0;ux;ux_x;uy;u2;u2_x]
%             if invar==ss1, or
%             (d/d invar)*P*[u0;ux;uy;u2] = dP*[u0;ux;uy;uy_y;u2;u2_y]
%             if invar==ss2
%           - if opts=='pure', we assume (d/d invar)*P*u does not depend on
%             any derivatives of u, and dP will be such that
%             (d/d invar)*P*[u0;ux;uy;u2] = dP*[u0;ux;uy;u2]
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - diff_opvar2d
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

% % Error checking
if ~isa(P,'opvar2d')
    error("Input diff function must be an opvar2d class object");
end
if isempty(invar)
    dP = P;
    return
end

if ndims(invar)>=3 || ~any(size(invar)==1)
    error('Differentiation variables must be input as an nx1 cell of strings, or nx1 polynomial')
end
   
% % Differentiate with respect to each element of invar separately
difvar = invar(1);
addvar = invar(2:end);
% For polynomial variable x^p*y^q, decompose into [x;x;...;x;y;...;y]
if isa(difvar,'polynomial')
    try
        double(difvar);
        isdble = 1;
    catch
        isdble = 0;
    end
    if isdble
        % If difvar is double, we're not differentiating
        [dP] = diff_opvar2d(P,addvar,opts);
        return
    else
        if any(difvar.matdim~=1) || any(difvar.coef~=1)
            error('Each element of polynomial variable array must correspond to a single variable with coefficient 1')
        end
        if length(difvar.varname)>1
            addvar = [difvar;addvar];
            difdeg = difvar.degmat;
            adddeg = difvar.degmat;
            difdeg(2:end) = 0;
            adddeg(1) = 0;
            difvar.degmat = difdeg;
            addvar.degmat = adddeg;
        end
        if difvar.degmat>1
            repval = difvar.degmat-1;
            difvar.degmat = 1;
            addvar = [repmat(difvar,[repval,1]);addvar];
        end
        difvar = difvar.varname;
    end
end

if ~any(ismember(P.var1.varname,difvar))
    error('Each differentiation variable must be one of the variables appearing in the opvar2d object')
end

% % % Perform the differentiation % % %
opvar2d dP;
dP.I = P.I;
dim_new = P.dim;
if strcmp(P.var1(1).varname,difvar) || strcmp(P.var2(1).varname,difvar)
    % The differentiation variable is the first variable
    sss = P.var1(1);    ttt = P.var2(1);
    if nargin==3 && (strcmpi(opts,'exclude') || strcmpi(opts,'pure'))
        % In this case, our state [u0; ux; uy; u2] doesn't change
        if any(~isequal(dP.Rxx{1},0)) || any(~isequal(dP.Rx2{1},0)) || any(~isequal(dP.R2x{1},0)) || ...
           any(~isequal(dP.R22{1,1},0)) || any(~isequal(dP.R22{1,2},0)) || any(~isequal(dP.R22{1,3},0))
            error('"Pure differentiation" of PI operator is only possible when there are no multiplier terms along the dimension of differentiation')
        else
            dP.dim = dim_new;
            
            % Take the derivatives for the map to L2[x]
            if dim_new(2,1)~=0
            if dim_new(1,2)~=0
            dP.Rx0 = diff(P.Rx0,sss);
            end
            if dim_new(2,2)~=0
            dP.Rxx{1} = subs(P.Rxx{2}-P.Rxx{3},ttt,sss);
            dP.Rxx{2} = diff(P.Rxx{2},sss);
            dP.Rxx{3} = diff(P.Rxx{3},sss);
            end
            if dim_new(3,2)~=0
            dP.Rxy = diff(P.Rxy,sss);
            end
            if dim_new(4,2)~=0
            dP.Rx2{1} = subs(P.Rx2{2}-P.Rx2{3},ttt,sss);
            dP.Rx2{2} = diff(P.Rx2{2},sss);
            dP.Rx2{3} = diff(P.Rx2{3},sss);
            end
            end
            % Take the derivatives for the map to L2[x,y]
            if dim_new(4,1)~=0
            if dim_new(1,2)~=0
            dP.R20 = diff(P.R20,sss);
            end
            if dim_new(2,2)~=0
            dP.R2x{1} = subs(P.R2x{2}-P.R2x{3},ttt,sss);
            dP.R2x{2} = diff(P.R2x{2},sss);
            dP.R2x{3} = diff(P.R2x{3},sss);
            end
            if dim_new(3,2)~=0
            dP.R2y{1} = diff(P.R2y{1},sss);
            dP.R2y{2} = diff(P.R2y{2},sss);
            dP.R2y{3} = diff(P.R2y{3},sss);
            end
            if dim_new(4,2)~=0
            dP.R22{1,1} = subs(P.R22{2,1}-P.R22{3,1},ttt,sss);
            dP.R22{1,2} = subs(P.R22{2,2}-P.R22{3,2},ttt,sss);
            dP.R22{1,3} = subs(P.R22{2,3}-P.R22{3,3},ttt,sss);
            dP.R22{2,1} = diff(P.R22{2,1},sss);
            dP.R22{2,2} = diff(P.R22{2,2},sss);
            dP.R22{2,3} = diff(P.R22{2,3},sss);
            dP.R22{3,1} = diff(P.R22{3,1},sss);
            dP.R22{3,2} = diff(P.R22{3,2},sss);
            dP.R22{3,3} = diff(P.R22{3,3},sss);
            end
            end
        end
    else
        % In this case, our state extends to [u0; ux; ux_x; uy; u2; u2_x];
        dim_new(2,2) = 2*dim_new(2,2);  % state [ux] becomes [ux; ux_x];
        dim_new(4,2) = 2*dim_new(4,2);  % state [u2] becomes [u2; u2_x];
        dP.dim = dim_new;
        
        % Take the derivatives for the map to L2[x]
        if dim_new(2,1)~=0
        if dim_new(1,2)~=0
        dP.Rx0 = diff(P.Rx0,sss);
        end
        if dim_new(2,2)~=0
        dP.Rxx{1} = [diff(P.Rxx{1},sss) + subs(P.Rxx{2}-P.Rxx{3},ttt,sss), P.Rxx{1}];
        dP.Rxx{2} = [diff(P.Rxx{2},sss) zeros(size(P.Rxx{2}))];
        dP.Rxx{3} = [diff(P.Rxx{3},sss) zeros(size(P.Rxx{3}))];
        end
        if dim_new(3,2)~=0
        dP.Rxy = diff(P.Rxy,sss);
        end
        if dim_new(4,2)~=0
        dP.Rx2{1} = [diff(P.Rx2{1},sss) + subs(P.Rx2{2}-P.Rx2{3},ttt,sss), P.Rx2{1}];
        dP.Rx2{2} = [diff(P.Rx2{2},sss) zeros(size(P.Rx2{2}))];
        dP.Rx2{3} = [diff(P.Rx2{3},sss) zeros(size(P.Rx2{3}))];
        end
        end
        % Take the derivatives for the map to L2[x,y]
        if dim_new(4,1)~=0
        if dim_new(1,2)~=0
        dP.R20 = diff(P.R20,sss);
        end
        if dim_new(2,2)~=0
        dP.R2x{1} = [diff(P.R2x{1},sss) + subs(P.R2x{2}-P.R2x{3},ttt,sss), P.R2x{1}];
        dP.R2x{2} = [diff(P.R2x{2},sss) zeros(size(P.R2x{2}))];
        dP.R2x{3} = [diff(P.R2x{3},sss) zeros(size(P.R2x{3}))];
        end
        if dim_new(3,2)~=0
        dP.R2y{1} = diff(P.R2y{1},sss);
        dP.R2y{2} = diff(P.R2y{2},sss);
        dP.R2y{3} = diff(P.R2y{3},sss);
        end
        if dim_new(4,2)~=0
        dP.R22{1,1} = [diff(P.R22{1,1},sss) + subs(P.R22{2,1}-P.R22{3,1},ttt,sss), P.R22{1,1}];
        dP.R22{1,2} = [diff(P.R22{1,2},sss) + subs(P.R22{2,2}-P.R22{3,2},ttt,sss), P.R22{1,2}];
        dP.R22{1,3} = [diff(P.R22{1,3},sss) + subs(P.R22{2,3}-P.R22{3,3},ttt,sss), P.R22{1,3}];
        dP.R22{2,1} = [diff(P.R22{2,1},sss) zeros(size(P.R22{2,1}))];
        dP.R22{2,2} = [diff(P.R22{2,2},sss) zeros(size(P.R22{2,2}))];
        dP.R22{2,3} = [diff(P.R22{2,3},sss) zeros(size(P.R22{2,3}))];
        dP.R22{3,1} = [diff(P.R22{3,1},sss) zeros(size(P.R22{3,1}))];
        dP.R22{3,2} = [diff(P.R22{3,2},sss) zeros(size(P.R22{3,2}))];
        dP.R22{3,3} = [diff(P.R22{3,3},sss) zeros(size(P.R22{3,3}))];
        end
        end
    end
    
else
    
    % The differentiation variable is the second variable
    sss = P.var1(2);    ttt = P.var2(2);
    if nargin==3 && (strcmpi(opts,'exclude') || strcmpi(opts,'pure'))
        % In this case, our state [u0; ux; uy; u2] doesn't change
        if any(~isequal(dP.Ryy{1},0)) || any(~isequal(dP.Ry2{1},0)) || any(~isequal(dP.R2y{1},0)) || ...
           any(~isequal(dP.R22{1,1},0)) || any(~isequal(dP.R22{2,1},0)) || any(~isequal(dP.R22{3,1},0))
            error('"Pure differentiation" of PI operator is only possible when there are no multiplier terms along the dimension of differentiation')
        else
            dP.dim = dim_new;
            
            % Take the derivatives for the map to L2[y]
            dP.Ry0 = diff(P.Ry0,sss);
            
            dP.Ryy{1} = subs(P.Ryy{2}-P.Ryy{3},ttt,sss);
            dP.Ryy{2} = diff(P.Ryy{2},sss);
            dP.Ryy{3} = diff(P.Ryy{3},sss);
            
            dP.Ryx = diff(P.Ryx,sss);
            
            dP.Ry2{1} = subs(P.Ry2{2}-P.Ry2{3},ttt,sss);
            dP.Ry2{2} = diff(P.Ry2{2},sss);
            dP.Ry2{3} = diff(P.Ry2{3},sss);
            
            % Take the derivatives for the map to L2[x,y]
            dP.R20 = diff(P.R20,sss);
            
            dP.R2x{1} = diff(P.R2x{1},sss);
            dP.R2x{2} = diff(P.R2x{2},sss);
            dP.R2x{3} = diff(P.R2x{3},sss);
            
            dP.R2y{1} = subs(P.R2y{2}-P.R2y{3},ttt,sss);
            dP.R2y{2} = diff(P.R2y{2},sss);
            dP.R2y{3} = diff(P.R2y{3},sss);
            
            dP.R22{1,1} = subs(P.R22{1,2}-P.R22{1,3},ttt,sss);
            dP.R22{2,1} = subs(P.R22{2,2}-P.R22{2,3},ttt,sss);
            dP.R22{3,1} = subs(P.R22{3,2}-P.R22{3,3},ttt,sss);
            dP.R22{1,2} = diff(P.R22{1,2},sss);
            dP.R22{2,2} = diff(P.R22{2,2},sss);
            dP.R22{3,2} = diff(P.R22{3,2},sss);
            dP.R22{1,3} = diff(P.R22{1,3},sss);
            dP.R22{2,3} = diff(P.R22{2,3},sss);
            dP.R22{3,3} = diff(P.R22{3,3},sss);
        end
    else
        % In this case, our state extends to [u0; ux; uy; uy_y; u2; u2_y];
        dim_new(3,2) = 2*dim_new(3,2);  % state [uy] becomes [uy; uy_y];
        dim_new(4,2) = 2*dim_new(4,2);  % state [u2] becomes [u2; u2_y];
        dP.dim = dim_new;
        
        % Take the derivatives for the map to L2[y]
        dP.Ry0 = diff(P.Ry0,sss);
        
        dP.Ryy{1} = [diff(P.Ryy{1},sss) + subs(P.Ryy{2}-P.Ryy{3},ttt,sss), P.Ryy{1}];
        dP.Ryy{2} = [diff(P.Ryy{2},sss) zeros(size(P.Ryy{2}))];
        dP.Ryy{3} = [diff(P.Ryy{3},sss) zeros(size(P.Ryy{3}))];
        
        dP.Ryx = diff(P.Ryx,sss);
        
        dP.Ry2{1} = [diff(P.Ry2{1},sss) + subs(P.Ry2{2}-P.Ry2{3},ttt,sss), P.Ry2{1}];
        dP.Ry2{2} = [diff(P.Ry2{2},sss) zeros(size(P.Ry2{2}))];
        dP.Ry2{3} = [diff(P.Ry2{3},sss) zeros(size(P.Ry2{3}))];
        
        % Take the derivatives for the map to L2[x,y]
        dP.R20 = diff(P.R20,sss);
        
        dP.R2x{1} = diff(P.R2x{1},sss);
        dP.R2x{2} = diff(P.R2x{2},sss);
        dP.R2x{3} = diff(P.R2x{3},sss);
        
        dP.R2y{1} = [diff(P.R2y{1},sss) + subs(P.R2y{2}-P.R2y{3},ttt,sss), P.R2y{1}];
        dP.R2y{2} = [diff(P.R2y{2},sss) zeros(size(P.R2y{2}))];
        dP.R2y{3} = [diff(P.R2y{3},sss) zeros(size(P.R2y{3}))];
        
        dP.R22{1,1} = [diff(P.R22{1,1},sss) + subs(P.R22{1,2}-P.R22{1,3},ttt,sss), P.R22{1,1}];
        dP.R22{2,1} = [diff(P.R22{2,1},sss) + subs(P.R22{2,2}-P.R22{2,3},ttt,sss), P.R22{2,1}];
        dP.R22{3,1} = [diff(P.R22{3,1},sss) + subs(P.R22{3,2}-P.R22{3,3},ttt,sss), P.R22{3,1}];
        dP.R22{1,2} = [diff(P.R22{1,2},sss) zeros(size(P.R22{1,2}))];
        dP.R22{2,2} = [diff(P.R22{2,2},sss) zeros(size(P.R22{2,2}))];
        dP.R22{3,2} = [diff(P.R22{3,2},sss) zeros(size(P.R22{3,2}))];
        dP.R22{1,3} = [diff(P.R22{1,3},sss) zeros(size(P.R22{1,3}))];
        dP.R22{2,3} = [diff(P.R22{2,3},sss) zeros(size(P.R22{2,3}))];
        dP.R22{3,3} = [diff(P.R22{3,3},sss) zeros(size(P.R22{3,3}))];
    
    end
    
end

% Get rid of terms which are practically zero
dP = zremove(dP);

% Continue differentiation with respect to remaining variables
if ~isempty(addvar) && nargin==3
    dP = diff_opvar2d(dP,addvar,opts);
elseif ~isempty(addvar)
    dP = diff_opvar2d(dP,addvar);
end

end