function [dP] = diff(P,invar,deg,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dP = diff(P) calculates spatial derivative of (Pu)(s) with respect to
% variable var so that d(P*u) = dP*[u;du]
% 
% INPUT 
%   P: opvar2d object
%   invar: nx1 cell or polynomial
%           - Each element invar{i} or invar(i) should correspond to a
%             a variable appearing in the opvar object
%           - If invar(i) = x^p for some power p, differentiation wrt x
%             will be performed p times
%   deg:    nx1 array of integers, defining to what order P should be
%           differentiated with respect to each variable
%   opts:   Optional inputs (see OUTPUT for details)
% 
% OUTPUT 
%   dP: composition of P with a spatial derivative wrt variable invar:
%           - dP will be such that
%             (d/d invar)*P*[u0;u1] = dP*[u0;u1;u1_x]
%           - if opts=='pure', we assume (d/d invar)*P*u does not depend on
%             any derivatives of u, and dP will be such that
%             (d/d invar)*P*[u0;u1] = dP*[u0;u1]
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - diff
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
% DJ, 01/24/2026: Correction in check if opts='pure' is supported.

% % Error checking
if nargin>4
    error('At most 4 inputs are accepted')
end
if isempty(invar)
    dP = P;
    return
end
if nargin==2 
    if isnumeric(invar)
        deg = invar;
        invar = P.var1;
        opts = '';
    elseif ischar(invar)
        deg = 1;
        opts = invar;
        invar = P.var1;
    else
        deg = 1;
        opts = '';
    end
elseif nargin==3
    if isnumeric(invar)
        opts = deg;
        deg = invar;
        invar = P.var1;
    elseif isa(deg,'char')
        opts = deg;
    else
        opts = '';
    end
end
if ~all(size(deg)==size(invar)) && numel(deg)~=1
    error('Number of orders of derivatives should match the number of variables.')
else
    invar = invar.^deg;
end

% if ndims(invar)>=3 || ~any(size(invar)==1)
%     error('Differentiation variables must be input as an nx1 cell of strings, or nx1 polynomial')
% end
   

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
            error('The variable of differentiation must match the primary variable of the opvar object')
        end
        
        % % % Perform the differentiation % % %
        opvar dP;
        dP.I = P.I;  
        dP.var1 = P.var1;   dP.var2 = P.var2;
        dim_new = P.dim;
        if strcmp(P.var1(1).varname,difvar) || strcmp(P.var2(1).varname,difvar)
            % The differentiation variable is the first variable
            sss = P.var1(1);    ttt = P.var2(1);
            if (strcmpi(opts,'exclude') || strcmpi(opts,'pure'))
                % In this case, our state [u0; ux; uy; u2] doesn't change
                P = clean_opvar(P,1e-12);
                if any(any(~isequal(P.R.R0,zeros(size(P.R.R0)))))
                    error('"Pure differentiation" of PI operator is only possible when there are no multiplier terms along the dimension of differentiation')
                else
                    dP.dim = dim_new;
                    
                    % Take the derivatives for the map to L2[x]
                    if dim_new(2,1)~=0
                        if dim_new(1,2)~=0
                            dP.Q2 = diff(polynomial(P.Q2),sss);
                        end
                        if dim_new(2,2)~=0
                            dP.R.R0 = subs(polynomial(P.R.R1-P.R.R2),ttt,sss);
                            dP.R.R1 = diff(polynomial(P.R.R1),sss);
                            dP.R.R2 = diff(polynomial(P.R.R2),sss);
                        end
                    end
                end
            else
                % In this case, our state extends to [u0; u1; u1_x];
                dim_new(2,2) = 2*dim_new(2,2);  % state [ux] becomes [ux; ux_x];
                dP.dim = dim_new;
                
                % Take the derivatives for the map to L2[x]
                if dim_new(2,1)~=0
                    if dim_new(1,2)~=0
                        dP.Q2 = diff(polynomial(P.Q2),sss);
                    end
                    if dim_new(2,2)~=0
                        dP.R.R0 = [diff(polynomial(P.R.R0),sss) + subs(polynomial(P.R.R1-P.R.R2),ttt,sss), P.R.R0];
                        dP.R.R1 = [diff(polynomial(P.R.R1),sss) zeros(size(P.R.R1))];
                        dP.R.R2 = [diff(polynomial(P.R.R2),sss) zeros(size(P.R.R2))];
                    end
                end
            end            
        else
            % The differentiation variable is some other variable which the
            % operator may or may not depend on
            sss = polynomial(difvar);
            if (strcmpi(opts,'exclude') || strcmpi(opts,'pure'))
                % In this case, we assume the state [u0; u1] doesn't depend
                % on the variable wrt which we differentiate
                % d/dx *Pop*u = (d/dx*Pop)*u
                dP.dim = dim_new;

                dP.P = diff(polynomial(P.P,sss));
                dP.Q1 = diff(polynomial(P.Q1,sss));
                dP.Q2 = diff(polynomial(P.Q2,sss));
                dP.R.R0 = diff(polynomial(P.R.R0,sss));
                dP.R.R1 = diff(polynomial(P.R.R1,sss));
                dP.R.R2 = diff(polynomial(P.R.R2,sss));
                
            else
                % In this case, we assume the state does depend on the
                % variable, and we extend it to to [u0; u1; u1_x];
                dim_new(2,2) = 2*dim_new(2,2);  % state [u1] becomes [u1; u1_x];
                dP.dim = dim_new;

                dP.P = diff(polynomial(P.P,sss));
                dP.Q1 = [diff(polynomial(P.Q1,sss)), Q1];
                
                dP.Q2 = diff(polynomial(P.Q2,sss));
                dP.R.R0 = [diff(polynomial(P.R.R0,sss)),P.R.R0];
                dP.R.R1 = [diff(polynomial(P.R.R1,sss)),P.R.R1];
                dP.R.R2 = [diff(polynomial(P.R.R2,sss)),P.R.R2];                
            end
        end
        P = dP; % Differentiate dP wrt remaining variables addvar
    end 
end

% Get rid of terms which are practically zero
tol = 1e-14;
dP = clean_opvar(P,tol);

end