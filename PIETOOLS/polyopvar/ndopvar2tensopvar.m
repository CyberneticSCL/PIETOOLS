function Top = ndopvar2tensopvar(Pop,type)
% TOP = NDOPVAR2TENSOPVAR(POP) returns a 'tensopvar_new' object TOP
% representing the same operator as the input 'nopvar' or 'ndopvar' object
% POP
%
% INPUTS
% - Pop:    m x n 'nopvar' or 'ndopvar' object representing a PI operator;
%
% OUTPUTS
% - Top:    m x n 'tensopvar' object representing the same operator as the
%           input;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ndopvar2tensopvar
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
% DJ, 04/14/2026: Initial coding

% Check the input
if ~isa(Pop,'nopvar') && ~isa(Pop,'ndopvar') && ~isa(Pop,'cell')
    error("Operator must be specified as object of type 'nopvar' or 'ndopvar'.")
end
% By default assume pointwise product in row dimension, and Kronecker in
% the column dimension
if nargin<=1
    type = [false,true];
elseif isscalar(type)
    type = [type,true];
end

if ~isa(Pop,'cell')
    % Declare a tensor-PI operator involving just a single factor
    Top = tensopvar();
    Top.ops = {Pop};
    Top.type = type;
    % Set spatial variables and domain
    Top.vars = Pop.vars;
    Top.dom = Pop.dom;
    Top.depmat1 = ones(1,size(Top.vars,1));
    Top.depmat2 = ones(1,size(Top.vars,1));
    % Set the order of the operators
    Top.order = 1;
else
    % Assume a cell of terms and factors is specified
    Top = tensopvar();
    Top.ops = Pop;
    Top.type = type;
    % Set spatial variables and domain
    varname = cell(0,2);
    dom = zeros(0,2);
    for j=1:size(Pop,2)
        if ~isa(Pop{1,j},'nopvar') && ~isa(Pop{1,j},'ndopvar')
            error("Elements of the cell must all be of type 'nopvar' or 'ndopvar'")
        end
        % Establish in which variables the operator is defined
        vars_j = pvar2varname(Pop{1,j}.vars);
        dom_j = Pop{1,j}.dom;
        % Add to the full list
        varname = [varname; vars_j];
        dom = [dom; dom_j];
        % Remove redundant variables
        [~,idcs] = unique(varname(:,1),'stable');
        varname = varname(idcs,:);
        dom = dom(idcs,:);
    end
    Top.vars = polynomial(varname);
    Top.dom = dom;
    % Set the dependency arrays
    depmat1 = zeros(size(Pop,2),size(Top.vars,1));
    depmat2 = zeros(size(Pop,2),size(Top.vars,1));
    for j=1:size(Pop,2)
        % Establish in which variables the operator is defined
        var1_j = pvar2varname(Pop{1,j}.vars(:,1));
        var_idx = ismember(varname(:,1),var1_j);
        depmat1(j,var_idx) = 1;
        depmat2(j,var_idx) = depmat2(j,var_idx)+1;
    end
    Top.depmat1 = depmat1;
    Top.depmat2 = depmat2;
    % Set the order of the operators
    Top.order = 1:size(Pop,2);
end

end