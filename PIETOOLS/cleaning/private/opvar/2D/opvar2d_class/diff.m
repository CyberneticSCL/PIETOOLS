function [dP] = diff(P,invar,deg,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [[dP] = diff(P) calculates spatial derivative of (Pu)(s) with respect to
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
%   deg:    nx1 array of integers, defining to what order P should be
%           differentiated with respect to each variable
%   opts:   Optional inputs (see OUTPUT for details)
% 
% OUTPUT 
%   dP: composition of P with a spatial derivative wrt variable invar:
%           - dP will be such that
%             (d/d invar)*P*[u0;ux;uy;u2] = dP*[u0;ux;ux_x;uy;u2;u2_x]
%             if invar==s1, or
%             (d/d invar)*P*[u0;ux;uy;u2] = dP*[u0;ux;uy;uy_y;u2;u2_y]
%             if invar==s2
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
% Initial coding DJ - 08/27/2021
% DJ, 01/24/2026: Correction in check if opts='pure' is supported.

% % Error checking
if ~isa(P,'opvar2d')
    error("Input diff function must be an opvar2d class object");
end
if isempty(invar)
    dP = P;
    return
else
    P = poly_opvar2d(P);   % Convert all parameters to polynomial for differentiation
end
if ndims(invar)>=3 || ~any(size(invar)==1)
    error('Differentiation variables must be input as an nx1 cell of strings, or nx1 polynomial')
end
if nargin>4
    error('At most 4 inputs are accepted')
elseif nargin==4
    invar = invar.^deg;
elseif nargin==3 && isa(deg,'char')
    opts = deg;
elseif nargin==3 && all(size(deg)==size(invar))
    opts = '';
    invar = invar.^deg;
elseif nargin==3
    error('An order of derivative must be specified for each variable')
end
   

% % Differentiate with respect to each element of invar separately
addvar = invar(:);
while ~isempty(addvar)
    difvar = addvar(1);     % Differentiate wrt diffvar now
    addvar = addvar((2:end)'); % Differentiate wrt remaining variables later
    
    % For polynomial variable x^p*y^q, decompose into [x;x;...;x;y;...;y]
    if isa(difvar,'polynomial')
        % First, see if the variable is a constant
        try
            double(difvar);
            isdble = 1;
        catch
            isdble = 0;
        end
        % If the variable is constant, we will not differentiate
        if ~isdble
            if any(difvar.coef~=1) || size(difvar.degmat,1)>1
                error('Each element of polynomial variable array must correspond to a single variable with coefficient 1')
            end
            if length(difvar.varname)>1
                % If difvar contains multiple variables, extract only first variable for differentiation now
                newvar = difvar;    % Variables to differentiate wrt later
                
                difdeg = difvar.degmat;
                difdeg(2:end) = 0;      % For now, differentiate only wrt first variable
                difvar.degmat = difdeg;
                
                newdeg = newvar.degmat;
                newdeg(1) = 0;          % Differentiate wrt remaining later
                newvar.degmat = newdeg;
                addvar = [newvar;addvar];
            end
            if difvar.degmat>1
                % If we want to differentiate wrt x^p, differentiate p times with respect to x
                repval = difvar.degmat-1;
                difvar.degmat = 1;
                addvar = [repmat(difvar,[repval,1]);addvar];
            end
            difvar = difvar.varname;
        end
    end
    
    % Next, differentiate with respect to difvar, unless difvar is a
    % constant
    if ~isdble        
        if ~any(ismember(P.var1.varname,difvar))
            error('Each differentiation variable must be one of the variables appearing in the opvar2d object')
        end
        
        % % % Perform the differentiation % % %
        opvar2d dP;
        dP.I = P.I;  
        dP.var1 = P.var1;   dP.var2 = P.var2;
        dim_new = P.dim;
        if strcmp(P.var1(1).varname,difvar) || strcmp(P.var2(1).varname,difvar)
            % The differentiation variable is the first variable
            sss = P.var1(1);    ttt = P.var2(1);
            if (strcmpi(opts,'exclude') || strcmpi(opts,'pure'))
                % In this case, our state [u0; ux; uy; u2] doesn't change
                P = clean_opvar(P,1e-12);
                if any(any(~isequal(P.Rxx{1},zeros(size(P.Rxx{1}))))) ||...
                     any(any(~isequal(P.Rx2{1},zeros(size(P.Rx2{1}))))) ||...
                       any(any(~isequal(P.R2x{1},zeros(size(P.R2x{1}))))) || ...
                         any(any(~isequal(P.R22{1,1},zeros(size(P.R22{1,1}))))) ||...
                           any(any(~isequal(P.R22{1,2},zeros(size(P.R22{1,2}))))) ||...
                             any(any(~isequal(P.R22{1,3},zeros(size(P.R22{1,3})))))
                    error('"Pure differentiation" of PI operator is only possible when there are no multiplier terms along the dimension of differentiation')
                else
                    dP.dim = dim_new;
                    
                    % Take the derivatives for the map to L2[x]
                    if dim_new(2,1)~=0
                        if dim_new(1,2)~=0
                            dP.Rx0 = diff(P.Rx0,sss);
                        end
                        if dim_new(2,2)~=0
                            dP.Rxx{1} = polynomial(subs(P.Rxx{2}-P.Rxx{3},ttt,sss));
                            dP.Rxx{2} = diff(P.Rxx{2},sss);
                            dP.Rxx{3} = diff(P.Rxx{3},sss);
                        end
                        if dim_new(3,2)~=0
                            dP.Rxy = diff(P.Rxy,sss);
                        end
                        if dim_new(4,2)~=0
                            dP.Rx2{1} = polynomial(subs(P.Rx2{2}-P.Rx2{3},ttt,sss));
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
                            dP.R2x{1} = polynomial(subs(P.R2x{2}-P.R2x{3},ttt,sss));
                            dP.R2x{2} = diff(P.R2x{2},sss);
                            dP.R2x{3} = diff(P.R2x{3},sss);
                        end
                        if dim_new(3,2)~=0
                            dP.R2y{1} = diff(P.R2y{1},sss);
                            dP.R2y{2} = diff(P.R2y{2},sss);
                            dP.R2y{3} = diff(P.R2y{3},sss);
                        end
                        if dim_new(4,2)~=0
                            dP.R22{1,1} = polynomial(subs(P.R22{2,1}-P.R22{3,1},ttt,sss));
                            dP.R22{1,2} = polynomial(subs(P.R22{2,2}-P.R22{3,2},ttt,sss));
                            dP.R22{1,3} = polynomial(subs(P.R22{2,3}-P.R22{3,3},ttt,sss));
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
                        dP.Rxx{1} = [diff(P.Rxx{1},sss) + polynomial(subs(P.Rxx{2}-P.Rxx{3},ttt,sss)), P.Rxx{1}];
                        dP.Rxx{2} = [diff(P.Rxx{2},sss) zeros(size(P.Rxx{2}))];
                        dP.Rxx{3} = [diff(P.Rxx{3},sss) zeros(size(P.Rxx{3}))];
                    end
                    if dim_new(3,2)~=0
                        dP.Rxy = diff(P.Rxy,sss);
                    end
                    if dim_new(4,2)~=0
                        dP.Rx2{1} = [diff(P.Rx2{1},sss) + polynomial(subs(P.Rx2{2}-P.Rx2{3},ttt,sss)), P.Rx2{1}];
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
                        dP.R2x{1} = [diff(P.R2x{1},sss) + polynomial(subs(P.R2x{2}-P.R2x{3},ttt,sss)), P.R2x{1}];
                        dP.R2x{2} = [diff(P.R2x{2},sss) zeros(size(P.R2x{2}))];
                        dP.R2x{3} = [diff(P.R2x{3},sss) zeros(size(P.R2x{3}))];
                    end
                    if dim_new(3,2)~=0
                        dP.R2y{1} = diff(P.R2y{1},sss);
                        dP.R2y{2} = diff(P.R2y{2},sss);
                        dP.R2y{3} = diff(P.R2y{3},sss);
                    end
                    if dim_new(4,2)~=0
                        dP.R22{1,1} = [diff(P.R22{1,1},sss) + polynomial(subs(P.R22{2,1}-P.R22{3,1},ttt,sss)), P.R22{1,1}];
                        dP.R22{1,2} = [diff(P.R22{1,2},sss) + polynomial(subs(P.R22{2,2}-P.R22{3,2},ttt,sss)), P.R22{1,2}];
                        dP.R22{1,3} = [diff(P.R22{1,3},sss) + polynomial(subs(P.R22{2,3}-P.R22{3,3},ttt,sss)), P.R22{1,3}];
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
            if (strcmpi(opts,'exclude') || strcmpi(opts,'pure'))
                % In this case, our state [u0; ux; uy; u2] doesn't change
                P = clean_opvar(P,1e-12);
                if any(any(~isequal(P.Ryy{1},zeros(size(P.Ryy{1}))))) ||...
                     any(any(~isequal(P.Ry2{1},zeros(size(P.Ry2{1}))))) ||...
                       any(any(~isequal(P.R2y{1},zeros(size(P.R2y{1}))))) || ...
                         any(any(~isequal(P.R22{1,1},zeros(size(P.R22{1,1}))))) ||...
                           any(any(~isequal(P.R22{2,1},zeros(size(P.R22{2,1}))))) ||...
                             any(any(~isequal(P.R22{3,1},zeros(size(P.R22{3,1})))))
                    error('"Pure differentiation" of PI operator is only possible when there are no multiplier terms along the dimension of differentiation')
                else
                    dP.dim = dim_new;
                    
                    % Take the derivatives for the map to L2[y]
                    dP.Ry0 = diff(P.Ry0,sss);
                    
                    dP.Ryy{1} = polynomial(subs(P.Ryy{2}-P.Ryy{3},ttt,sss));
                    dP.Ryy{2} = diff(P.Ryy{2},sss);
                    dP.Ryy{3} = diff(P.Ryy{3},sss);
                    
                    dP.Ryx = diff(P.Ryx,sss);
                    
                    dP.Ry2{1} = polynomial(subs(P.Ry2{2}-P.Ry2{3},ttt,sss));
                    dP.Ry2{2} = diff(P.Ry2{2},sss);
                    dP.Ry2{3} = diff(P.Ry2{3},sss);
                    
                    % Take the derivatives for the map to L2[x,y]
                    dP.R20 = diff(P.R20,sss);
                    
                    dP.R2x{1} = diff(P.R2x{1},sss);
                    dP.R2x{2} = diff(P.R2x{2},sss);
                    dP.R2x{3} = diff(P.R2x{3},sss);
                    
                    dP.R2y{1} = polynomial(subs(P.R2y{2}-P.R2y{3},ttt,sss));
                    dP.R2y{2} = diff(P.R2y{2},sss);
                    dP.R2y{3} = diff(P.R2y{3},sss);
                    
                    dP.R22{1,1} = polynomial(subs(P.R22{1,2}-P.R22{1,3},ttt,sss));
                    dP.R22{2,1} = polynomial(subs(P.R22{2,2}-P.R22{2,3},ttt,sss));
                    dP.R22{3,1} = polynomial(subs(P.R22{3,2}-P.R22{3,3},ttt,sss));
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
                
                dP.Ryy{1} = [diff(P.Ryy{1},sss) + polynomial(subs(P.Ryy{2}-P.Ryy{3},ttt,sss)), P.Ryy{1}];
                dP.Ryy{2} = [diff(P.Ryy{2},sss) zeros(size(P.Ryy{2}))];
                dP.Ryy{3} = [diff(P.Ryy{3},sss) zeros(size(P.Ryy{3}))];
                
                dP.Ryx = diff(P.Ryx,sss);
                
                dP.Ry2{1} = [diff(P.Ry2{1},sss) + polynomial(subs(P.Ry2{2}-P.Ry2{3},ttt,sss)), P.Ry2{1}];
                dP.Ry2{2} = [diff(P.Ry2{2},sss) zeros(size(P.Ry2{2}))];
                dP.Ry2{3} = [diff(P.Ry2{3},sss) zeros(size(P.Ry2{3}))];
                
                % Take the derivatives for the map to L2[x,y]
                dP.R20 = diff(P.R20,sss);
                
                dP.R2x{1} = diff(P.R2x{1},sss);
                dP.R2x{2} = diff(P.R2x{2},sss);
                dP.R2x{3} = diff(P.R2x{3},sss);
                
                dP.R2y{1} = [diff(P.R2y{1},sss) + polynomial(subs(P.R2y{2}-P.R2y{3},ttt,sss)), P.R2y{1}];
                dP.R2y{2} = [diff(P.R2y{2},sss) zeros(size(P.R2y{2}))];
                dP.R2y{3} = [diff(P.R2y{3},sss) zeros(size(P.R2y{3}))];
                
                dP.R22{1,1} = [diff(P.R22{1,1},sss) + polynomial(subs(P.R22{1,2}-P.R22{1,3},ttt,sss)), P.R22{1,1}];
                dP.R22{2,1} = [diff(P.R22{2,1},sss) + polynomial(subs(P.R22{2,2}-P.R22{2,3},ttt,sss)), P.R22{2,1}];
                dP.R22{3,1} = [diff(P.R22{3,1},sss) + polynomial(subs(P.R22{3,2}-P.R22{3,3},ttt,sss)), P.R22{3,1}];
                dP.R22{1,2} = [diff(P.R22{1,2},sss) zeros(size(P.R22{1,2}))];
                dP.R22{2,2} = [diff(P.R22{2,2},sss) zeros(size(P.R22{2,2}))];
                dP.R22{3,2} = [diff(P.R22{3,2},sss) zeros(size(P.R22{3,2}))];
                dP.R22{1,3} = [diff(P.R22{1,3},sss) zeros(size(P.R22{1,3}))];
                dP.R22{2,3} = [diff(P.R22{2,3},sss) zeros(size(P.R22{2,3}))];
                dP.R22{3,3} = [diff(P.R22{3,3},sss) zeros(size(P.R22{3,3}))];
                
            end
        end
        P = dP; % Differentiate dP wrt remaining variables addvar
    end 
end

% Get rid of terms which are practically zero
tol = 1e-14;
dP = clean_opvar(P,tol);

% % Continue differentiation with respect to remaining variables
% if ~isempty(addvar) && exist('opts','var')
%     dP = diff_opvar2d(dP,addvar,opts);
% elseif ~isempty(addvar)
%     dP = diff_opvar2d(dP,addvar);
% end

end