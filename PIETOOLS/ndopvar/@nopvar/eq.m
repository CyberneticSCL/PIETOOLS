function logval = eq(P1,P2,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logval = eq(P1,P2,tol) tests equality of two dopvars P1 and P2 with tolerance tol
% Date: 6/13/19
% Version: 1.0
% 
% INPUT
% P1, P2: nopvar class objects
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

if ~isa(P1,'nopvar')&&(P1==0)
    P1 = 0*P2;
    % ndopvar P1
    % P1.dom = P2.dom;
    % P1.dim = P2.dim;
    % P1.vars = P2.vars;
elseif  ~isa(P2,'nopvar')&&(P2==0) 
    P2 = 0*P1;
    % ndopvar P2
    % P2.dom = P1.dom;
    % P2.dim = P1.dim;
    % P2.vars = P1.vars;
elseif ~isa(P1,'nopvar')|| ~isa(P2,'nopvar')
    error('To check equality either both values must be nopvar objects, or one of them have to be zero');
end



logval = true;
if any(P1.dim~=P2.dim)|| any(any(P1.dom~= P2.dom))
    disp('nopvars have different dimensions or domain and hence cannot be equal');
    logval = false;
    return
end

if ~all(strcmp(P1.vars(:, 1).varname, cell2mat(P2.vars(:, 1).varname)))
    disp('nopvars have different output variables and hence cannot be equal');
    logval = false;
    return
end



if any(P1.deg(:)~=P2.deg(:))
    disp('nopvar have different degree. Convert to the same degree and compare ==')
    logval = false;
    return
    % nopvar P1p P2p
    % P1p.dom = P1.dom;
    % P2p.dom = P2.dom;
    % 
    % P1p.dim = P1.dim;
    % P2p.dim = P1.dim;
    % 
    % P1p.deg = max(P1.deg, P2.deg); % common degree base
    % P2p.deg = max(P1.deg, P2.deg); 
    % 
    % P1p.vars = P1.vars;
    % P2p.vars = P2.vars;
    % 
    % % construct degmat for common basis of monomials
    % for dim = 1:DIM
    %     Zd1_out{dim} = monomials(P1.vars(dim, 1), P1.deg);
    %     Zd1_in{dim} = monomials(P1.vars(dim, 2), P1.deg);
    %     Zd2_out{dim} = monomials(P2.vars(dim, 1), P2.deg);
    %     Zd2_in{dim} = monomials(P2.vars(dim, 2), P2.deg);
    %     Zdfull_out{dim} = monomials(P2.vars(dim, 1), max(P1.deg, P2.deg));
    %     Zdfull_in{dim} = monomials(P2.vars(dim, 2), max(P1.deg, P2.deg));
    % end
    % 
    % % construct degmat for P1, P2 for rows and columns
    % % left_degmat1= kron(eye(P1.dim(1)), Zd11.degmat);
    % % right_degmat1= kron(eye(P1.dim(2)), Zd21.degmat);
    % % left_degmat2= kron(eye(P2.dim(1)), Zd12.degmat);
    % % right_degmat2= kron(eye(P2.dim(2)), Zd22.degmat);
    % % left_degmat_full = kron(eye(P2.dim(1)), Zd1_full.degmat);
    % % right_degmat_full= kron(eye(P2.dim(2)), Zd2_full.degmat);
    % % 
    % % % compare degmat for common basis and P1, P2 operators
    % % [~, Loc_l1] =ismember(left_degmat1,  left_degmat_full, 'rows');
    % % [~, Loc_r1] =ismember(right_degmat1, right_degmat_full, 'rows');
    % % [~, Loc_l2] =ismember(left_degmat2,  left_degmat_full, 'rows');
    % % [~, Loc_r2] =ismember(right_degmat2, right_degmat_full, 'rows');
    % 
    % % define indeces for nopvar
    % DIM = size(P1.vars, 1); % domain 
    % binStr = dec2base(0:(3^DIM-1), 3); % indeces for P.C
    % 
    % for iter_idx = 1:size(binStr, 1)
    %     substr = binStr(iter_idx, :);
    %     Qindx = str2num(substr')'+1; % array of .C indeces
    %     Qindx_cell = num2cell(Qindx); % cell array of .C indeces
    %     index_for_monom_in_theta = Qindx ~= 1; % include monomials for theta
    % 
    % end
    % for ii=1:numel(P1.C)
    %     % check if coefficients are small  
    %     C1p = P1.C{ii};
    %     C2p = P2.C{ii};
    %     n = max(size(C1p, 1), size(C2p, 1));
    %     m = max(size(C1p, 2), size(C2p, 2));
    %     P1p.C{ii} = sparse(n, m);
    %     P1p.C{ii}(Loc_l1, Loc_r1) = C1p;
    % 
    %     P2p.C{ii} = sparse(n, m);
    %     P2p.C{ii}(Loc_l2, Loc_r2) = C2p;
    % end
else
    P1p = P1;
    P2p = P2;
end


if logval % check only if Dopvar have the same dimensions
    for ii=1:numel(P1p.C)
        diff = P1p.C{ii} - P2p.C{ii};
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