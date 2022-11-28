function PIE = initialize(PIE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE = initialize(PIE) takes a "pie_struct" object and
% checks that all the necessary fields are appropriately specified, and
% assigns a default value to all the optional fields that have not been
% specified.
% 
% INPUT
%   PIE: "pie_struct" class object.
%
% OUTPUT
%   PIE: "pie_struct" class object describing the same system as the input.

% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Initial coding SS - 08/01/2022
% Adjust and add as class file function, DJ - 11/01/2022

% Initialize a spatial domain for the PIE.
if isempty(PIE.dom)
    if isa(PIE.T,'opvar') || isa(PIE.T,'opvar2d')
        % Use the domain specified for operator T
        PIE.dom = PIE.T.I;
    elseif ~isnan(PIE.dim)
        % Assume default domain [0,1]^D
        PIE.dom = [zeros(PIE.dim,1),ones(PIE.dim,1)];
    end
end
% Initialize spatial variables for the PIE.
if isempty(PIE.vars)
    if isa(PIE.T,'opvar') || isa(PIE.T,'opvar2d')
        % Use the domain specified for operator T
        PIE.vars = [PIE.T.var1,PIE.T.var2];
    elseif ~isnan(PIE.dim)
        % Assume default variables s theta
        if PIE.dim==1
            PIE.vars = polynomial({'s','theta'});
        elseif PIE.dim==2
            PIE.vars = polynomial({'ss1','tt1';'ss2','tt2'});
        elseif PIE.dim>=3
            error('PIEs in more than 2 dimensions are currently not supported.')
        end
    end
end
% Check the spatial dimension of the PIE.
if ~isnan(PIE.dim) && PIE.dim>2
    error('PIEs in more than 2 dimensions are currently not supported.')
end


% % The operators T and A are mandatory!
% Check that operator T has been appropriately specified.
if ~isa(PIE.T,'opvar') && ~isa(PIE.T,'opvar2d')
    error("PIE.T is a mandatory element of PIE structure and cannot be defaulted to zero");
else
    if PIE.dim==1 && isa(PIE.T,'opvar2d')
        warning("PI operators for a 1D PIE should be specified as opvar objects; attempting to convert T to opvar...")
        PIE.T = opvar2d2opvar(PIE.T);
    elseif PIE.dim==2 && isa(PIE.T,'opvar')
        warning("PI operators for a 2D PIE should be specified as opvar2d objects; attempting to convert T to opvar2d...")
        PIE.T = opvar2d2opvar(PIE.T);
    end
    if any(size(PIE.T.I)~=size(PIE.dom))
        error("The spatial dimension of the operator T does not match the spatial dimension of the PIE.");
    elseif any(any(PIE.T.I~=PIE.dom))
        error("The spatial domain of the operator T does not match the spatial domain of the PIE.");
    elseif ~all(all(isequal([PIE.T.var1,PIE.T.var2],PIE.vars)))
        error("The spatial variables of the operator T do not match those of the PIE.")
    end
    if any(PIE.T.dim(:,1)~=PIE.T.dim(:,2))
        error("PIE.T operator must be a map between symmetric spaces");
    end
    nx = PIE.T.dim(:,1);
end
% Check that operator A has been appropriately specified.
if ~isa(PIE.A,'opvar') && ~isa(PIE.A,'opvar2d')
    error("PIE.A is a mandatory element of PIE structure and cannot be defaulted to zero");
else
    if PIE.dim==1 && isa(PIE.A,'opvar2d')
        warning("PI operators for a 1D PIE should be specified as opvar objects; attempting to convert A to opvar...")
        PIE.A = opvar2d2opvar(PIE.A);
    elseif PIE.dim==2 && isa(PIE.A,'opvar')
        warning("PI operators for a 2D PIE should be specified as opvar2d objects; attempting to convert A to opvar2d...")
        PIE.A = opvar2d2opvar(PIE.A);
    end
    if any(size(PIE.A.I)~=size(PIE.dom))
        error("The spatial dimension of the operator A does not match the spatial dimension of the PIE.");
    elseif any(any(PIE.A.I~=PIE.dom))
        error("The spatial domain of the operator A does not match the spatial domain of the PIE.");
    elseif ~all(all(isequal([PIE.A.var1,PIE.A.var2],PIE.vars)))
        error("The spatial variables of the operator A do not match those of the PIE.")
    end
    if any(any(PIE.A.dim(:,1)~=PIE.T.dim))
        error("The dimensions of the operator A should match those of the operator T.");
    end
end

% Find remaining dimensions
nz = find_dim(PIE,1,{'C1';'D11';'D12'});
ny = find_dim(PIE,1,{'C2';'D21';'D22'});
nw = find_dim(PIE,2,{'Tw';'B1';'D11';'D21'});
nu = find_dim(PIE,2,{'Tu';'B2';'D12';'D22'});


% % Finally, loop over all operators, checking that their domains,
% % variables, and dimensions are appropriate.
opnames = {'T', 'Tw', 'Tu';
           'A', 'B1', 'B2';
           'C1','D11','D12';
           'C2','D21','D22'};
row_dims = {nx; nx; nz; ny};
row_obj = {'state','x'; 'state','x'; 'output','z'; 'output','z'};
col_dims = {nx; nw; nu};
col_obj = {'state','x'; 'input','w'; 'input','u'};

n_ops = numel(opnames);
for kk=1:n_ops
    [ii,jj] = ind2sub(size(opnames),kk);
    rdim = row_dims{ii};
    robj = [row_obj{ii,1}," ",row_obj{ii,2}];
    cdim = col_dims{jj};
    cobj = [col_obj{jj,1}," ",col_obj{jj,2}];
    opname = opnames{kk};
    opval = PIE.(opname);
    if (isa(opval,'opvar') || isa(opval,'opvar2d')) && ~(opval==0)
        % A non-zero operator has been specified
        % --> check that domain, vars, and dimensions are appropriate.
        if PIE.dim==1 && isa(opval,'opvar2d')
            warning(["PI operators for a 1D PIE should be specified as opvar objects; attempting to convert ",opname," to opvar..."])
            opval = opvar2d2opvar(opval);
        elseif PIE.dim==2 && isa(opval,'opvar')
            warning(["PI operators for a 2D PIE should be specified as opvar2d objects; attempting to convert ",opname," to opvar2d..."])
            opval = opvar2d2opvar(opval);
        end
        if any(size(opval.I)~=size(PIE.dom))
            error(["The spatial dimension of the operator ",opname," does not match the spatial dimension of the PIE."]);
        elseif any(any(opval.I~=PIE.dom))
            error(["The spatial domain of the operator ",opname," does not match the spatial domain of the PIE."]);
        elseif ~all(all(isequal([opval.var1,opval.var2],PIE.vars)))
            error(["The spatial variables of the operator ",opname," do not match those of the PIE."])
        end
        if any(opval.dim(:,1)~=rdim)
            error(["The output (row) dimensions of the operator ",opname," do not match the expected dimensions of the ",robj,"."]);
        elseif any(opval.dim(:,2)~=cdim)
            error(["The input (column) dimensions of the operator ",opname," do not match the expected dimensions of the ",cobj,"."]);
        end
    else
        % No (non-zero) operator has been specified
        % --> initialize a zero operator of appropriate dimensions.
        if PIE.dim<=1
            opval = opvar();
        else
            opval = opvar2d();
        end
        opval.var1 = PIE.vars(:,1);     opval.var2 = PIE.vars(:,2);
        opval.I = PIE.dom;              opval.dim = [rdim,cdim];
    end
    PIE.(opname) = opval;
end
end



%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
function dim = find_dim(PIE,idx,op_names)
% Establish row or column dimensions of the operators "op_names" in the PIE
% structure. Returns row dimensions if idx==1, returns column dimensions if
% idx==2. Only checks dimensions of nonzero operators. 

% Assume all zero dimensions to start.
dim = zeros(2^max(PIE.dim,1),1);

% Loop over the operators, updating dimension if the operator is nonzero.
n_ops = numel(op_names);
for kk=1:n_ops
    opname = op_names{kk};
    opval = PIE.(opname);
    if (isa(opval,'opvar') || isa(opval,'opvar2d')) && ~(opval==0)
        % A non-zero operator has been specified
        % --> extract the row dimensions
        if PIE.dim==1 && isa(opval,'opvar2d')
            warning(["PI operators for a 1D PIE should be specified as opvar objects; attempting to convert ",opname," to opvar..."])
            opval = opvar2d2opvar(opval);
        elseif PIE.dim==2 && isa(opval,'opvar')
            warning(["PI operators for a 2D PIE should be specified as opvar2d objects; attempting to convert ",opname," to opvar2d..."])
            opval = opvar2d2opvar(opval);
        end
        if all(dim==0)
            % If no dimensions were previously known, set the dimensions.
            dim = opval.dim(:,idx);
            opname_1 = opname;
        elseif ~all(opval.dim(:,idx)==dim)
            % If dimensions were already known, the value should match.
            if idx==1
                error(["The output (row) dimensions of the operator ",opname," should match those of the operator ",opname_1,"."]);
            else
                error(["The input (column) dimensions of the operator ",opname," should match those of the operator ",opname_1,"."]);
            end
        end
    end
end
end