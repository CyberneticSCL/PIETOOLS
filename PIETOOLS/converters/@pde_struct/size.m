function [n_eqs,n_terms] = size(PDE,objs,opts)
% [N_EQS,N_TERMS] = SIZE(PDE,OBJS,OPTS) returns the total number of
% equations and maximal number of terms in these equations for a given PDE
% structure.
%
% INPUT
% - PDE:        'pde_struct' object of which to determine the size.
% - objs:       'char', 'cell' or numeric object. 
%               If numeric, must be one of 1 and 2, specifying whether to 
%               return the number of equations "n_eqs" (1), or number of 
%               terms "n_terms" (2) in single output case.
%               If 'char', must be one of 'all', 'free', 'x', 'y', 'z',
%               'w', 'u', or 'BC', specifying for which object to return
%               the number of equations and terms. 
%               If 'cell', each element must be one of 'free', 'x', 'y',
%               'z', 'w', 'u' or 'BC', specifying for which objects to
%               return the number of equations and terms.
% - opts:       'char' object, set to either 'vec_size' or 'vec_size_tot'.
%
% OUTPUT
% - n_eqs:      Number of equations that appear in the PDE structure. If
%               "obj" is of type 'char', the output corresponds only to the
%               number of equations of type "obj". If opts='vec_size_tot',
%               the number of equations will account for the row-size of
%               each of these equations. If opts='vec_size', and a single
%               object is specified, the output with be an array with
%               number of elements matching the number of equations of the
%               type specified by "obj", and each element providing the
%               size of the vector defining this equation. If "obj" is a
%               cell specifying multiple objects, the output will be a cell
%               as well, with element i providing the output associated
%               to equation type obj{i}.
% 
% [n_eqs,n_terms] = size(PDE)    returns the total number of equations
% of any type ('free', 'x', 'y', 'z', and 'BC'), as well as the maximal
% number of terms that appears in any of these equations. The sizes do not
% account for the fact that terms may be vector-valued
%
% [n_eqs,n_terms] = size(PDE,obj) returns the total number equations of
% type "obj", if "obj" is one of 'free', 'x', 'y', 'z', or 'BC', as well as
% the maximal number of terms that appears in the object. If "obj" is one
% of 'u' or 'w', n_eqs will just correspond to the number of components of
% this type, and n_terms will be 0. If "obj" is a cell array specifying
% multiply object types, the outputs will be standard arrays of the same
% size as "obj", with each element specifying the total number of equations
% and maximum number of terms for the associated object.
%
% [n_eqs,n_terms] = size(PDE,obj,'vec_size_tot') returns the same output as
% [n_eqs,n_terms] = size(PDE,obj), but now with the output "n_eqs"
% providing the number of equations accounting for the vector-valued size
% of these equations. The output n_terms will be the same.
%
% [n_eqs,n_terms] = size(PDE,obj,'vec_size') returns an output "n_eqs" that
% is an array, with the number of elements corresponding to the number of
% equations of type "obj" in the PDE, and each element corresponding to the
% vector size of the terms that appear in each of these equations. If "obj"
% is a cell specifying multiple object types, "n_eqs" and "n_terms" will be
% cells as well, with each element yielding the specified output for the
% given object.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 06/23/2024

% % % Process the inputs
% % Check that the PDE is properly specified
if ~isa(PDE,'pde_struct')
    error("This size function should only be called for 'pde_struct' class objects.")
end

% % Check that the PDE object for which to determine the size is properly
% % specified.
use_all = false;    % Set to true to count all types of equations, without distinguishing.
use_dim = 0;        % Set to 1 or 2 to return only the number of equations, or number of terms, respectively.
use_cell = false;   % Set to true to return the output as a cell, which is necessary if 'objs' is a cell and 'vec_size' is used.
if nargin<2
    % Determine the size for every type of equation.
    if ~isempty(PDE.free)
        % If free terms are specified, the other objects should not contain
        % any equations, so we don't count them.
        objs = {'free'};
    else
        objs = {'x';'y';'z';'BC'};
    end
    use_all = true;
else
    if isnumeric(objs)
        % objs=1 or 2, specifying which of the outputs is desired
        if ~numel(objs)==1 || ~ismember(objs,[1;2])
            error("Second argument of 'size' can be set to either 1 or 2 to return the number of equations or maximal number of terms, respectively.")
        end
        use_dim = objs;
        if ~isempty(PDE.free)
            objs = {'free'};
        else
            objs = {'x';'y';'z';'BC'};
        end
        use_all = true;
    else
        if ~ischar(objs) && ~iscellstr(objs)
            error("Components of the PDE of which to determine the size should be specified as a cell of 'char' objects.")
        elseif (ischar(objs) && strcmpi(objs,'all')) || (iscellstr(objs) && numel(objs)==1 && strcmpi(objs{1},'all'))
            if ~isempty(PDE.free)
                % If free terms are specified, the other objects should not
                % contain any equations, so we don't count them.
                objs = {'free'};
            else
                objs = {'x';'y';'z';'BC'};
            end
            use_all = true;
        elseif any(~ismember(objs,{'free';'x';'y';'z';'u';'w';'BC'}))
            error("The PDE component of which to return the size should be one of 'free', 'x', 'y', 'z', 'u', 'w', or 'BC'.")
        elseif iscell(objs)
            use_cell = true;
        elseif ischar(objs)
            objs = {objs};
        end
    end
end

% % Check that the third input makes sense.
use_vec_size = false;
use_vec_size_tot = false;
if nargin>=3
    if ischar(opts) && (strcmpi(opts,'vec_size') || strcmpi(opts,'vec_size_tot'))
        use_vec_size = true;
    else
        error("The third argument should be a 'char' specifying either 'vec_size' or 'vec_size_tot'.")
    end
    if strcmpi(opts,'vec_size_tot')
        use_vec_size_tot = true;
    end
end



% % % Determine the size 
n_eqs = cell(length(objs),1);
n_terms = cell(length(objs),1);
for kk=1:numel(objs)
    obj = objs{kk};
    % % For now, keep track of the number of rows and number of terms in
    % each individual equation.
    n_eqs{kk} = ones(numel(PDE.(obj)),1);
    n_terms{kk} = zeros(numel(PDE.(obj)),1);
    for ii=1:numel(PDE.(obj))
        % Check the row-size of the equation, if desired.
        if use_vec_size && isfield(PDE.(obj){ii},'size')
            n_eqs{kk}(ii) = PDE.(obj){ii}.size;
        end
        % Check the number of terms that appear in the equation.
        if isfield(PDE.(obj){ii},'term')
            n_terms{kk}(ii) = max(numel(PDE.(obj){ii}.term));
        end        
    end
end



% % % Post-process the outputs.

% % If 'all' option is used, do not distinguish between objects of
% % different types, but combine information.
if use_all
    n_eqs = cell2mat(n_eqs);        n_terms = cell2mat(n_terms);
    if any(n_terms)
        n_eqs(n_terms==0) = 0;      % if no equation is specified, don't count it as equation.
    end
    n_eqs = {n_eqs};            n_terms = {n_terms};
end

% % If the individual vector sizes of the equations are not desired, sum
% % over these size to get the total size.
if ~use_vec_size || use_vec_size_tot
    % Sum over the sizes of each equation to get the total number of
    % equations of each type.
    n_eqs = cellfun(@(a) sum(a),n_eqs);
    % Determine the maximal number of terms that appear in the equations of
    % each type.
    n_terms = cellfun(@(a) max(a),n_terms);
end

% Convert cell to array.
if isa(n_eqs,'cell') && ~use_cell
    n_eqs = n_eqs{1};   n_terms = n_terms{1};
end

if nargin==1 && (nargout==1 || nargout==0)
    % For single input single output, return size as 1x2 array.
    n_eqs = [n_eqs,n_terms];
elseif nargout==1 && use_dim==2
    % If desired direction of size is set to 2, return only maximal number
    % of terms.
    n_eqs = n_terms;
end

end