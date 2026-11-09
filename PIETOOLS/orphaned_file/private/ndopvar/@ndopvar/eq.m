function logval = eq(P1,P2,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logval = eq(P1,P2,tol) tests equality of two dopvars P1 and P2 with tolerance tol
% Date: 6/13/19
% Version: 1.0
% 
% INPUT
% P1, P2: ndopvar class objects
% tol: acceptable tolerance value. If max(P1-P2)<tol, then P1=P2
% 
% OUTPUT
% logval: returns 1 if the objects are equal, 0 if not equal
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - eq
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
% AT, 01/21/2026: Initial coding

if nargin<3
    tol=1e-14;
end    

if ~isa(P1,'ndopvar')&&(P1==0)
    P1 = 0*P2;
    % ndopvar P1
    % P1.I = P2.I;
    % P1.dim = P2.dim;
    % P1.vars = P2.vars;
elseif  ~isa(P2,'ndopvar')&&(P2==0) 
    P2 = 0*P1;
    % ndopvar P2
    % P2.I = P1.I;
    % P2.dim = P1.dim;
    % P2.vars = P1.vars;
elseif (isa(P1,'nopvar')&&isa(P2,'ndopvar'))||(isa(P1,'ndopvar')&&isa(P2,'nopvar'))||(isa(P1,'ndopvar')&&isa(P2,'ndopvar'))
    if any(any(P1.dom~=P2.dom))||any(P1.dim(:)~=P2.dim(:)) 
        error('Operators begin compared do not have same interval or have a mismatch in dimensions');
    end
elseif ~isa(P1,'ndopvar')|| ~isa(P2,'ndopvar')
    error('To check equality either both values must be ndopvar objects, or one of them have to be zero');
end


logval = true;
if any(P1.dim~=P2.dim)
    disp('Dopvars have different dimensions and hence cannot be equal');
    logval = false;
    return
end

% if the degrees are different convert to the same
if any(P1.deg(:)~=P2.deg(:))
    max_degree = max(P1.deg, P2.deg);
    P1p = change_degree(P1, max_degree);
    P2p = change_degree(P2, max_degree);
    P1 = P1p;
    P2 = P2p;
end

if isa(P1, 'ndopvar')
    Aop_dvarname = string(P1.dvarname);
    Aop_dvarname = strtrim(Aop_dvarname);
else
    Aop_dvarname = {};
end
if isa(P2, 'ndopvar')
    Bop_dvarname = string(P2.dvarname);
    Bop_dvarname = strtrim(Bop_dvarname);
else
    Bop_dvarname = {};
end

if numel(Aop_dvarname) ~= numel(Bop_dvarname) || ~isequal(Aop_dvarname,Bop_dvarname)
    dvars1 = Aop_dvarname; % convert array to char array
    dvars2 = Bop_dvarname;
    % numberOfCharacters = max(size(dvars1, 2), size(dvars2, 2));
    % dvars1 = pad(dvars1, numberOfCharacters); % pad with ' ' if needed
    % dvars2 = pad(dvars2, numberOfCharacters); % pad with ' ' if needed 
    if isempty(dvars1) || isempty(dvars2)
        common_dvar = [];
        full_dvars = [dvars1; dvars2];
    else
        common_dvar = intersect(dvars1, dvars2, 'rows');
        new_dvars   = setdiff(dvars2, common_dvar, 'rows');
        full_dvars = [dvars1; new_dvars];
    end
    P1 = change_dec_var(P1, full_dvars);
    P2 = change_dec_var(P2, full_dvars);
end


if logval
    for ii=1:numel(P1.C)

        diff = P1.C{ii}-P2.C{ii}; 
        % check if coefficients are small 
        if max(max(abs(diff))) < tol
            logval = true;
        else
            logval = false;
        end
        % stop checking if we 
        if logval == false
            break
        end
    end
    

end
end