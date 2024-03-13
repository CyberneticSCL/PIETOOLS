function [isgood_D,isgood_Dpar,deg] = checkdeg_lpi_eq_2d(Pop,Dop,deg,p_indx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [isgood_D,isgood_Dpar,deg] = checkdeg_lpi_eq_2d(Pop,Dop,deg,p_indx)
% checks whether all monomials appearing in the parameters of a PI operator
% (or variable) Pop also appear in the PI variable Dop, so that the
% constraint Pop == Dop can be imposed.
% 
% INPUT
%   Pop:    opvar2d or dopvar2d object, on which to impose the (LPI)
%            equality constraint Pop==Dop.
%   Dop:    dopvar2d object, used to enforce the (LPI) equality constraint.
%   deg:    Degrees for "poslpivar_2d" used to define Dop.
%   p_indx: (Optional) 1D array of indices ranging between 1 and 16,
%           indicating which of the defining fields of Pop and Dop to
%           compare:
%           1 - R00;   5 - R0x;   9 - R0y;   13 - R02;
%           2 - Rx0;   6 - Rxx;  10 - Rxy;   14 - Rx2;
%           3 - Ry0;   7 - Ryx;  11 - Ryy;   15 - Ry2;
%           4 - R20;   8 - R2x;  12 - R2y;   16 - R22;
%           Defaults to all elements except R00, (2:16);
%
% OUTPUT 
%   isgood_D:       Binary value indicating if Dop has all necessary monomials
%   isgood_Dpar:    4x4 Binary array indicating if each parameter in Dop
%                   has all necessary monomials
%   deg:            Possibly adjusted degree specifications for 
%                   constructing a new operator Dop that (hopefully)
%                   contains all monomials in Pop.
%                   If isgood_D==1, the output deg will match the input
%                   deg.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - checkdeg_lpi_ineq_2d
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
% Initial coding DJ - 03_19_2022 

if nargin<3
    error('Insufficient inputs provided')
elseif nargin==3
    p_indx = (2:16)';       % Just check all the polynomial parameters
end

if ~isa(Pop,'dopvar2d') && ~isa(Pop,'opvar2d')
    error('Input operator on which to impose lpi constraint must be of type ''opvar2d'' or ''dopvar2d''');
elseif ~isa(Dop,'dopvar2d')
    error('Dummy operator used to enforce lpi constraint must be of type ''dopvar2d''');
end

params = Pop.params;        % Extract the names of the parameters defining Pop (and Dop)
dim = Pop.dim;              % Extract the dimensions of the operator
isgood_Dpar = ones(4,4);    % Logical array indicating which params in Pop can be matched by Dop             

% Loop over each of the requested parameters, and check if the parameter in
% Dop has all necessary monomials to match the same parameter in Pop
for j=1:length(p_indx)
    p = p_indx(j);                  % Index of desired parameter in params
    PR = Pop.(params{p});           % Desired parameter in Pop
    DR = Dop.(params{p});           % Associated parameter in Dop
    [row,col] = ind2sub([4,4],p);   % Position of parameter in the opvar object
    
    if ~isempty(PR) && ~(dim(row,1)==0 || dim(col,1)==0) && isa(PR,'cell')
        % Certain parameters contain integral and multiplier contributions
        for k=1:numel(PR)
            PR{k} = dpvar(PR{k});
            DR{k} = dpvar(DR{k});
            PRdeg = PR{k}.degmat;                   % Degrees of monomials in parameter PRk
            DRdeg = DR{k}.degmat;                   % Degrees of monomials in parameter DRk
            deg_dif = setdiff(PRdeg,DRdeg,'rows');  % Monomials in PRk that are not in DRk
            isgood_DRk = isempty(deg_dif);          % 1 if the DRk has all monomial also in PRk
            isgood_Dpar(p) = isgood_Dpar(p) && isgood_DRk; % Is parameter p okay?
            
            % If DR is missing monomials, increase the maximal monomial
            % degree. NOTE: Only Rxx, Ryy and R22 are considered here.
            % Off-diagonal degrees are adjusted later
            if ~isgood_DRk && row==col
                if row==2
                    deg.dx{k} = deg.dx{k}+1;
                elseif row==3
                    deg.dy{k} = deg.dy{k}+1;
                elseif row==4
                    deg.d2{k} = deg.d2{k}+1;
                end
            end
            
        end
    elseif ~isempty(PR) && ~(dim(row,1)==0 || dim(col,1)==0)
        PR = dpvar(PR);
        DR = dpvar(DR);
        PRdeg = PR.degmat;                          % Degrees of monomials in parameter PR
        DRdeg = DR.degmat;                          % Degrees of monomials in parameter DR
        deg_dif = setdiff(PRdeg,DRdeg,'rows');      % Monomials in PR that are not in DR
        isgood_DR = isempty(deg_dif);               % 1 if the DR has all monomial also in PR
        isgood_Dpar(p) = isgood_Dpar(p) && isgood_DR;   % Is parameter p okay?
    end
end

% Check if all parameters in Pop can be matched by those in Dop
isgood_D = all(isgood_Dpar(:));

% If all diagoanl parameters (Rxx, Ryy, R22) are fine, but any off-diagonal 
% parameter (e.g. R2x) cannot be matched, we have to puzzle a bit with
% which degrees to increase...
if ~isgood_D && all(diag(isgood_Dpar))
    isgood_Dpar_tmp = isgood_Dpar;
    while any(~isgood_Dpar_tmp(:))
        isbad_Dpar_in = sum(~isgood_Dpar_tmp,1)';   % How many parameters mapping from (R,L2[x],L2[y],L2[x,y])' are insufficient?
        isbad_Dpar_out = sum(~isgood_Dpar_tmp,2);   % How many parameters mapping to (R,L2[x],L2[y],L2[x,y])' are insufficient?
        
        % Check which input/output dimension is in greatest need of
        % additional monomials
        [~,deg_add_indx] = max(isbad_Dpar_in(2:end) + isbad_Dpar_out(2:end));
        deg_add_indx = deg_add_indx + 1;
        
        % Increase the maximal degree of the determined monomials
        if deg_add_indx==2
            for k=1:numel(deg.dx)
                deg.dx{k} = deg.dx{k}+1;
            end
        elseif deg_add_indx==3
            for k=1:numel(deg.dy)
                deg.dy{k} = deg.dy{k}+1;
            end
        elseif deg_add_indx==4
            for k=1:numel(deg.d2)
                deg.d2{k} = deg.d2{k}+1;
            end
        end
        
        % We assume the new degree offers sufficient freedom for the
        % considered input/output dimension
        isgood_Dpar_tmp(deg_add_indx,:) = 1;
        isgood_Dpar_tmp(:,deg_add_indx) = 1;
    end
end

end