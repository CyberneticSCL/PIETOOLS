function out = set(Pop,prop,val,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = set(obj,prop,val) takes in one opvar2d object obj, and returns a
% opvar2d object out with out.prop = val
%
% Version 2.0
% Date: 06/13/22
% 
% INPUT
% Pop:  opvar2d class object
% prop: char class object specifying a field name for the opvar2d
% val:  a value to be assigned to Pop.(prop)
% opts: optional char field that can be 'nosubs' to perform no substitution
%       of the variables in the parameters when adjusting the variables
%       Pop.var1 or Pop.var2, or can be 'nocheck' when setting a parameter
%       Pop.R.. to avoid performing a full check that the dimensions and
%       variables in the parameter match those of the operator.
%
% OUTPUT
% out:  opvar2d object with Pop.(prop)=val
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or D. Jagt at djagt@asu.edu
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 06/13/2022
% DJ, 08/09/2022 - Keep parameters polynomial after setting vars.


if nargin<4
    opts = [];
end
     
% % Call the appropriate subroutine to set the desired property
if strcmp(prop,'dom') || strcmp(prop,'I')
    out = Pop;
    out.I = val;    % I field has its own function in the class function file
    return
elseif strcmp(prop,'var1')
    out = set_var1(Pop,val,opts);
    return
elseif strcmp(prop,'var2')
    out = set_var2(Pop,val,opts);
    return
elseif strcmp(prop,'dim')
    out = Pop;
    out.dim = val;  % dim field has its own function in the class function file
    return
else
    Rparams = {'R00', 'R0x', 'R0y', 'R02';
               'Rx0', 'Rxx', 'Rxy', 'Rx2';
               'Ry0', 'Ryx', 'Ryy', 'Ry2';
               'R20', 'R2x', 'R2y', 'R22'};
    if ismember(prop,Rparams)
        out = set_param(Pop,prop,val,opts);
    else
        error(['The proposed "',prop,'" is not a settable field in ''opvar2d'' class objects'])        
    end
end

end



%%
function out = set_var1(Pop,val,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = set_var1(obj,val) takes in one opvar2d object Pop, and returns a
% opvar2d object out with out.var1 = val.
%
% Version 1.0
% Date: 06/13/22
% 
% INPUT
% Pop:  opvar2d class object;
% val:  a 2x1 pvar (polynomial class) object or cellstr, specifying the
%       (names of the) primary spatial variables of the operator;
% opts: (optional argument) char option 'nosubs'. If this option is
%       specified, the new variables will be assigned without further
%       adjustments (i.e. substitutions) of the operator parameters.
%
% OUTPUT
% out:  opvar2d object with Pop.var1 = val. Note that the parameters of Pop
%       will be adjusted to match the new variable names. That is, e.g.
%       out.Rxy = subs(Pop.Rxy,Pop.var1,val);
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or D. Jagt at djagt@asu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 06/13/2022
% DJ, 08/09/2022 - Keep parameters polynomial after substitution.

% If 'nosubs' option is passed, we assume the variables to be properly
% specified, and simply set them without performing any substitution
if nargin==3 && strcmp(opts,'nosubs')
    out = Pop;
    out.var1 = val;
    return
elseif nargin==3 && ~isempty(opts)
    error('Unrecognized option; third argument must be ''nosubs'' or empty.')
end

% Check that the variables are appropriately specified
if ~all(size(val)==[2,1])
    error('Variables should be specified as a 2x1 ''polynomial'' class object')
end
if iscellstr(val)
    val = polynomial(val);
end
if ~ispvar(val)
    error('Variables should be specified as a vector of 2 polynomial variables')
end
% Check that we get a unique set of variables
out = Pop;
var1 = Pop.var1;
if isequal(val(1),val(2))
    error('The variables along the x and y direction cannot match');
elseif any(isequal(val(1),out.var2)) || any(isequal(val(2),out.var2))
    error('The primary variables cannot match the dummy variables in the opvar2d object');
end
% If the desired variables do not match those currently present,
% replace the current variables with the desired ones
if any(~isequal(val,var1))
    Rparams = {'R00', 'R0x', 'R0y', 'R02';
               'Rx0', 'Rxx', 'Rxy', 'Rx2';
               'Ry0', 'Ryx', 'Ryy', 'Ry2';
               'R20', 'R2x', 'R2y', 'R22'};
    for k=1:numel(Rparams)
        PR = Pop.(Rparams{k});
        if isa(PR,'polynomial')
            out.(Rparams{k}) = polynomial(subs(PR,var1(~isequal(val,var1)),val(~isequal(val,var1))));
        elseif isa(PR,'cell')
            for l=1:numel(PR)
                PRR = PR{l};
                if isa(PRR,'polynomial')
                    PR{l} = polynomial(subs(PRR,var1(~isequal(val,var1)),val(~isequal(val,var1))));
                end
            end
            out.(Rparams{k}) = PR;
        end
    end
end
% Set the desired variables
out.var1 = val;

end



%%
function out = set_var2(Pop,val,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = set_var2(obj,val) takes in one opvar2d object Pop, and returns a
% opvar2d object out with out.var2 = val.
%
% Version 1.0
% Date: 06/13/22
% 
% INPUT
% Pop:  opvar2d class object;
% val:  a 2x1 pvar (polynomial class) object or cellstr, specifying the
%       (names of the) primary spatial variables of the operator;
% opts: (optional argument) char option 'nosubs'. If this option is
%       specified, the new variables will be assigned without further
%       adjustments (i.e. substitutions) of the operator parameters.
%
% OUTPUT
% out:  opvar2d object with Pop.var2 = val. Note that the parameters of Pop
%       will be adjusted to match the new variable names. That is, e.g.
%       out.R22{2,2} = subs(Pop.R22{2,2},Pop.var2,val);
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or D. Jagt at djagt@asu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 06/13/2022
% DJ, 08/09/2022 - Keep parameters polynomial after substitution.

% If 'nosubs' option is passed, we assume the variables to be properly
% specified, and simply set them without performing any substitution
if nargin==3 && strcmp(opts,'nosubs')
    out = Pop;
    out.var2 = val;
    return
elseif nargin==3 && ~isempty(opts)
    error('Unrecognized option; third argument must be ''nosubs'' or empty.')
end

% Check that the variables are appropriately specified
if ~all(size(val)==[2,1])
    error('Variables should be specified as a 2x1 polynomial')
end
if iscellstr(val)
    val = polynomial(val);
end
if ~ispvar(val)
    error('Variables should be specified as a vector of 2 polynomial variables')
end
% Check that we get a unique set of variables
out = Pop;
var2 = Pop.var2;
if isequal(val(1),val(2))
    error('The dummy variables along the x and y direction cannot match');
elseif any(isequal(val(1),out.var1)) || any(isequal(val(2),out.var1))
    error('The primary variables cannot match the dummy variables in the opvar2d object');
end
% If the desired dummy variables do not match those currently present,
% replace the current dummy variables with the desired ones
if any(~isequal(val,var2))
    RRparams = {'Rxx','R2x','Ryy','Ry2','Rx2','Ry2','R22'};
    for k=1:numel(RRparams)
        PR = Pop.(RRparams{k});
        for l=2:numel(PR)
            PRR = PR{l};
            if isa(PRR,'polynomial')
                PR{l} = polynomial(subs(PRR,var2(~isequal(val,var2)),val(~isequal(val,var2))));
            end
        end
        out.(RRparams{k}) = PR;
    end
end
% Set the desired dummy variables
out.var2 = val;

end



%%
function out = set_param(Pop,PR_name,PR_val,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = set_param(obj,PR_name,PR_val) takes in one opvar2d object Pop, and
% a opvar2d object out with out.(PR_name) = PR_val, if possible.
%
% Version 1.0
% Date: 06/13/22
% 
% INPUT
% Pop:      opvar2d class object;
% PR_name:  char object, specifying one of the following parameter names:
%                     {'R00', 'R0x', 'R0y', 'R02';
%                      'Rx0', 'Rxx', 'Rxy', 'Rx2';
%                      'Ry0', 'Ryx', 'Ryy', 'Ry2';
%                      'R20', 'R2x', 'R2y', 'R22'};
% PR_val:   double or polynomial class object if the parameter is one of
%           {'R00','R0x','R0y','R02','Rx0','Rxy','Ry0','Ryx','R20'}.
%           cell of double or polynomial class objects otherwise. PR_val
%           should specify the desired value of the parameter
%           Pop.(PR_name);
% opts:     (optional argument) char object 'nocheck', used to assign the
%           desired parameter without checking for the dimensions and
%           variables to match.
%
% OUTPUT
% out:  opvar2d object with Pop.(PR_name) = PR_val. 
%       The dimensions and variables of PR_val should match the expect
%       dimensions and variables as established from the operator.
%       If the dimensions do not match, but all other parameters in the
%       same row and/or column as PR_name in Pop are zeroes, these
%       0-parameters will be expanded/compressed to match the new
%       dimensions of the operator, as determined by the dimensions of
%       PR_val. Any other parameters will not be adjusted.
%       If the parameter introduces new or unexpected variables, checks are
%       performed to make sure these new variables do not conflict with
%       current variables, and the fields out.var1 and out.var2 are updated
%       accordingly (if sufficient variable information is available).
%       In any other case of conflicting dimensions or variables, an error
%       will be produced.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or D. Jagt at djagt@asu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 13/06/2022


% List of parameters that may be adjusted
Rparams = {'R00', 'R0x', 'R0y', 'R02';
           'Rx0', 'Rxx', 'Rxy', 'Rx2';
           'Ry0', 'Ryx', 'Ryy', 'Ry2';
           'R20', 'R2x', 'R2y', 'R22'};

% % If no check is requested, simply assign the parameter without checking
% % the variables and dimensions to match.
if strcmp(opts,'nocheck')
    out = Pop;
    [ii_loc,jj_loc] = find(ismember(Rparams,PR_name));
    out.(Rparams{ii_loc,jj_loc}) = PR_val;
    % If the desired parameter is a dpvar, convert the operator to a
    % dopvar2d object
    if isa(PR_val,'dpvar')
        out = opvar2dopvar2d(out);
    end
    return
end

% Otherwise, perform a relatively simple check to see if the dimensions and
% variables of the new parameter make sense with those of the operator.
[ii_loc,jj_loc] = find(ismember(Rparams,PR_name));
if isempty(PR_val)
    % If the parameter to replace is defined by a cell, assume the user
    % wants all of the cell elements to be empty.
    if isa(Pop.(Rparams{ii_loc,jj_loc}),'cell')
        PR_val = cell(size(Pop.(Rparams{ii_loc,jj_loc})));
    end
    % If the proposed parameter and current parameter are both empty,
    % adjust the dimensions of the proposed parameter.
    if any([Pop.dim(ii_loc,1),Pop.dim(jj_loc,2)]==0)
        if isa(PR_val,'cell')
            for k=1:numel(PR_val)
                PR_val{k} = polynomial(zeros([Pop.dim(ii_loc,1),Pop.dim(jj_loc,2)]));
            end
        else
            PR_val = polynomial(zeros([Pop.dim(ii_loc,1),Pop.dim(jj_loc,2)]));
        end
        ismatch_dim = true;
        ismatch_vars = true;
    else
        % If the current parameter is not empty, we'll have to see if we
        % can reconcile the dimensions of the proposed parameter with the
        % current shape of the opvar2d object.
        ismatch_dim = false;
        ismatch_vars = false;
    end
elseif isa(PR_val,'double')
    % Check that the dimensions of the parameter match those of the
    % operator.
    ismatch_dim = all(size(PR_val)==[Pop.dim(ii_loc,1),Pop.dim(jj_loc,2)]);
    ismatch_vars = true;    % No variables is always fine
elseif isa(PR_val,'polynomial') || isa(PR_val,'dpvar')
    if isa(Pop.(Rparams{ii_loc,jj_loc}),'cell')
        error(['The proposed parameter "',Rparams{ii_loc,jj_loc},'" should be specified as a cell.'])
    end
    % Check that the dimensions of the parameter match those of the
    % operator.
    ismatch_dim = all(size(PR_val)==[Pop.dim(ii_loc,1),Pop.dim(jj_loc,2)]);
    % Check which variables may appear in the proposed parameter.
    Param_vars = {polynomial([]),          Pop.var1(1), Pop.var1(2),    Pop.var1;
                  Pop.var1(1),    polynomial([]),       Pop.var1,       polynomial([]);
                  Pop.var1(2),    Pop.var1,             polynomial([]), polynomial([]);
                  Pop.var1,       polynomial([]),       polynomial([]), polynomial([])};
    Param_vars = Param_vars{ii_loc,jj_loc};    
    % Check if the variables in the parameter match the expected variables.
    if isempty(PR_val.varname)
        ismatch_vars = true;
    else
        ismatch_vars = all(ismember(PR_val.varname,Param_vars.varname));
    end
elseif isa(PR_val,'cell')
    % Check that the parameters are appropriately specified.
    if ~isa(Pop.(Rparams{ii_loc,jj_loc}),'cell')
        error(['The proposed parameter "',Rparams{ii_loc,jj_loc},'" should be specified as an object of type "double" or "polynomial".'])
    elseif numel(PR_val)~=numel(Pop.(Rparams{ii_loc,jj_loc}))
        error(['The number of elements of the specified parameter do not match the number of elements of "Pop.',Rparams{ii_loc,jj_loc},'".'])
    else
        reshape(PR_val,size(Pop.(Rparams{ii_loc,jj_loc})));
    end
    % Establish which variables the parameter may depend on.
    if ii_loc==2 && jj_loc==2
        Param_vars = {Pop.var1(1);[Pop.var1(1);Pop.var2(1)];[Pop.var1(1);Pop.var2(1)]};
    elseif ii_loc==3 && jj_loc==3
        Param_vars = {Pop.var1(2),[Pop.var1(2);Pop.var2(2)],[Pop.var1(2);Pop.var2(2)]};
    elseif ii_loc==2 || jj_loc==2
        Param_vars = {Pop.var1; [Pop.var1;Pop.var2(1)]; [Pop.var1;Pop.var2(1)]};
    elseif ii_loc==3 || jj_loc==3
        Param_vars = {Pop.var1, [Pop.var1;Pop.var2(2)], [Pop.var1;Pop.var2(2)]};
    else
        Param_vars = {Pop.var1,               [Pop.var1;Pop.var2(2)], [Pop.var1;Pop.var2(2)];
                      [Pop.var1;Pop.var2(1)], [Pop.var1;Pop.var2],    [Pop.var1;Pop.var2];
                      [Pop.var1;Pop.var2(1)], [Pop.var1;Pop.var2],    [Pop.var1;Pop.var2]};
    end
    % For each of the elements of the parameter, make sure the dimensions
    % match those of the operator, and the variables match the expected
    % variables.
    ismatch_dim = true;
    ismatch_vars = true;
    for k=1:numel(PR_val)
        ismatch_dim = all(size(PR_val{k})==[Pop.dim(ii_loc,1),Pop.dim(jj_loc,2)]);
        if isa(PR_val{k},'polynomial') || isa(PR_val{k},'dpvar')
            ismatch_vars = ismatch_vars && all(ismember(PR_val{k}.varname,Param_vars{k}.varname));
        elseif ~isa(PR_val{k},'double')
            error(['Each element of the proposed parameter "',Rparams{ii_loc,jj_loc},'" should be specified as an object of type "double" or "polynomial".'])
        end
        if ~ismatch_dim || ~ismatch_vars
            break
        end        
    end
else
    error(['Parameters should be specified as objects of type "double" or "polynomial".'])
end

% If the dimensions and variables of the parameter make sense, we can set
% the new parameter, and return.
if ismatch_dim && ismatch_vars
    out = Pop;
    out.(Rparams{ii_loc,jj_loc}) = PR_val;
    % If the desired parameter is a dpvar, convert the operator to a
    % dopvar2d object
    if isa(PR_val,'dpvar')
        out = opvar2dopvar2d(out);
    end
    return
end

% % Otherwise, perform some more extensive checks:
% % To set the new parameter, we first determine the dimensions and spatial
% % variables of the current parameters. If a parmeter is all zeroes, we
% % will allow its dimension to be adjusted (if necessary) to match the
% dimensions of the new variable.

dim_out = zeros(size(Rparams));         % Output dimensions of the parameters    
dim_in = zeros(size(Rparams));          % Input dimensions of the parameters
iszero_param = false(size(Rparams));    % Logical array indicating if parameters are all zero
varname_list = cell(0,1);               % Full list of variable names that appear
varname_indx = zeros(0,4);              % logical array indicating to which variable each varname might correspond

for k=1:numel(Rparams)
    [ii,jj] = ind2sub(size(Rparams),k);
    PR = Pop.(Rparams{k});
    if isempty(PR)
        nr = size(PR,1);    nc = size(PR,2);
        dim_out(k) = nr;    iszero_param(k) = true;
        dim_in(k) = nc;
    elseif isa(PR,'double')
        nr = size(PR,1);    nc = size(PR,2);
        dim_out(k) = nr;    dim_in(k) = nc;
        if all(all(PR==0))
            iszero_param(k) = true;
        end
    elseif isa(PR,'polynomial')
        nr = size(PR,1);    nc = size(PR,2);
        dim_out(k) = nr;    dim_in(k) = nc;
        %varname_cell{k} = PR.varname;
        if all(all(isequal(PR,0)))
            iszero_param(k) = true;
        end
        varname_list = [varname_list; PR.varname];
        if ii==4 || jj==4 || (ii==2 && jj==3) || (ii==3 && jj==2)   % The parameter may depend on x and y
            varname_indx = [varname_indx; repmat([1,1,0,0],length(PR.varname),1)];
        elseif ii==2 || jj==2                                       % The parameter only depends on x
            varname_indx = [varname_indx; repmat([1,0,0,0],length(PR.varname),1)];
        elseif ii==3 || jj==3                                       % The parameter only depends on y
            varname_indx = [varname_indx; repmat([0,1,0,0],length(PR.varname),1)];
        end
        
    elseif isa(PR,'cell')
        % If PR corresponds to multiple parameters, check that their
        % dimensions match, and allow these dimensions to be modified only
        % if all of the parameters are zero
        dim_out_k = zeros(size(PR));            dim_in_k = zeros(size(PR));
        iszero_param_k = false(size(PR));
        for l=1:numel(PR)
            PRR = PR{l};
            [ii_l,jj_l] = ind2sub(size(PR),l);
            if isempty(PRR)
                nr = size(PRR,1);   nc = size(PRR,2);
                dim_out_k(l) = nr;
                dim_in_k(l) = nc;   iszero_param_k(l) = true;
            elseif isa(PRR,'double')
                nr = size(PRR,1);   nc = size(PRR,2);
                dim_out_k(l) = nr;  dim_in_k(l) = nc;
                if all(all(PRR==0))
                    iszero_param_k(l) = true;
                end
            elseif isa(PRR,'polynomial')
                nr = size(PRR,1);   nc = size(PRR,2);
                dim_out_k(l) = nr;  dim_in_k(l) = nc;
                if all(all(isequal(PRR,0)))
                    iszero_param_k(l) = true;
                end
                varname_list = [varname_list; PRR.varname];
                if ii==4 || jj==4
                    if ii_l==1 && jj_l==1
                        varname_indx = [varname_indx; repmat([1,1,0,0],length(PRR.varname),1)];
                    elseif jj_l==1
                        varname_indx = [varname_indx; repmat([1,1,1,0],length(PRR.varname),1)];
                    elseif ii_l==1
                        varname_indx = [varname_indx; repmat([1,1,0,1],length(PRR.varname),1)];
                    else
                        varname_indx = [varname_indx; repmat([1,1,1,1],length(PRR.varname),1)];
                    end
                elseif ii==2 || jj==2
                    if l==1
                        varname_indx = [varname_indx; repmat([1,0,0,0],length(PRR.varname),1)];
                    else
                        varname_indx = [varname_indx; repmat([1,0,1,0],length(PRR.varname),1)];
                    end
                elseif ii==3 || jj==3
                    if l==1
                        varname_indx = [varname_indx; repmat([0,1,0,0],length(PRR.varname),1)];
                    else
                        varname_indx = [varname_indx; repmat([0,1,0,1],length(PRR.varname),1)];
                    end
                else
                    warning(['Parameter "',Rparams{k},'" should not be a cell, how did this happen?']);
                end
            else
                warning(['Element ',num2str(l)',' of field "',Rparams{k},'" is not appropriately specified'])
            end
        end
        % Parameters that are not zero have a fixed dimension, which may
        % not be modified
        fixed_dim_out_k = dim_out_k(~iszero_param_k);
        fixed_dim_in_k = dim_in_k(~iszero_param_k);
        if isempty(fixed_dim_out_k)
            if all((dim_out_k(2:end)-dim_out_k(1:end-1))==0)
                dim_out(k) = dim_out_k(1);
            else
                dim_out(k) = nan;
            end
        elseif ~all((fixed_dim_out_k(2:end)-fixed_dim_out_k(1:end-1))==0)
            warning(['The row dimensions of parameters "',Rparams{k},'" do not match!'])
        else
            dim_out(k) = fixed_dim_out_k(1);
        end
        if isempty(fixed_dim_in_k)
            if all((dim_in_k(2:end)-dim_in_k(1:end-1))==0)
                dim_in(k) = dim_in_k(1);
            else
                dim_in(k) = nan;
            end
        elseif ~all((fixed_dim_in_k(2:end)-fixed_dim_in_k(1:end-1))==0)
            warning(['The column dimensions of parameters "',Rparams{k},'" do not match!'])
        else
            dim_in(k) = fixed_dim_in_k(1);
        end
        % Only allow the dimensions to be modified if all parameters PR{l}
        % are zero.
        iszero_param(k) = all(iszero_param_k(:));
    else
        error(['The parameter "Pop.',Rparams{k},'" is not properly specified.'])
    end
end
        
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % %
% % Given the input and output dimensions, establish dimensions of the
% % operator, allowing dimensions associated to only zero-parameters to be
% % modified.
Pdim = zeros(4,2);
for l=1:size(Pdim,1)
    % Row dimensions
    fixed_dim_out_j = dim_out(l,~iszero_param(l,:));
    if isempty(fixed_dim_out_j)
        if all((dim_out(l,2:end)-dim_out(l,1:end-1))==0)
            Pdim(l,1) = dim_out(l,1);
        else
            Pdim(l,1) = nan;
        end
    elseif ~all((fixed_dim_out_j(2:end)-fixed_dim_out_j(1:end-1))==0)
        warning('The row dimensions of parameters mapping to R^n do not match')
    else
        Pdim(l,1) = fixed_dim_out_j(1);
    end
    % Column dimensions
    fixed_dim_in_j = dim_in(~iszero_param(:,l),l);
    if isempty(fixed_dim_in_j)
        if all((dim_in(2:end,l)-dim_in(1:end-1,l))==0)
            Pdim(l,2) = dim_in(1,l);
        else
            Pdim(l,2) = nan;
        end
    elseif ~all((fixed_dim_in_j(2:end)-fixed_dim_in_j(1:end-1))==0)
        warning('The column dimensions of parameters mapping from R^n do not match')
    else
        Pdim(l,2) = fixed_dim_in_j(1);
    end
end
% Which dimensions correspond only to zero-parameters, i.e. are flexible to
% change?
isflex_Pdim = [all(iszero_param,2), all(iszero_param,1)'];


% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % %
% % Next, check the variables:
% Assign to each varname an associated integer varnum
[varname_unique_Pop,vindx,varnum_Pop] = unique(varname_list);
% Here "varname_unique = varname_list(vindx);" and 
% "varname_list = varname_unique(varnum);"
varnum_unique_Pop = varnum_Pop(vindx);     
if ~isempty(varnum_unique_Pop)
    max_varnum_Pop = max(varnum_unique_Pop);
else
    max_varnum_Pop = 0;
end
% Build a vartable:
% The first column of vartab indicates the considered variable number.
% The remaining columns are logical indices indicating which spatial
% variables this varnum may be associated to.
vartab_Pop = [varnum_Pop,varname_indx];
[vartab_Pop,~] = unique(vartab_Pop,'rows');
nvars_Pop = size(vartab_Pop,2)-1;

% If the new parameter introduces a new row to the vartab, we set
% adjust_vartab = true;
adjust_vartab = false;


% % % --------------------------------------------------------------- % % %
% % Having determined the dimensions and variables of the current operator,
% % establish the dimensions and variables of the new parameter.

% First, determine which parameter to adjust, and extract the associated
% dimensions in Pop.
[ii_loc,jj_loc] = find(ismember(Rparams,PR_name));
dim_ij = [Pdim(ii_loc,1),Pdim(jj_loc,2)];

% % We distinguish two cases: PR is a single parameter, or PR is a cell of
% % parameters (e.g. Pop.Rxx{i}).
if isa(PR_val,'double') || isa(PR_val,'polynomial') || isa(PR_val,'dpvar')
    % First case: PR is a single parameter
    if isa(Pop.(Rparams{ii_loc,jj_loc}),'cell')
        error(['Parameters "',Rparams{ii_loc,jj_loc},'" should be specified as a cell.'])
    end
    % If the desired parameter is a dpvar, convert the operator to a
    % dopvar2d object
    if isa(PR_val,'dpvar')
        isdopvar = true;
    else
        isdopvar = false;
    end
    
    % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % %
    % % Add the variables information of the parameter to the vartable
    if (isa(PR_val,'polynomial') || isa(PR_val,'dpvar')) && ~isempty(PR_val.varname)
        var_indx_val = [];      % binary values indicating which combination of variables a varname corresponds to
        if ((ii_loc==1 && jj_loc==2) || (ii_loc==2 && jj_loc==1))   % The parameter depends just on x
            if length(PR_val.varname)>1
                error(['The parameter "',Rparams{ii_loc,jj_loc},'" can vary in at most one variable.'])
            elseif ~isempty(PR_val.varname) && ~all(ismember(PR_val.varname,Pop.var1(1).varname))
                var_indx_val = [1,0,0,0];       % The variable name can only correspond to the x-variable
            end
        elseif ((ii_loc==1 && jj_loc==3) || (ii_loc==3 && jj_loc==1))   % The parameter depends just on y
            if length(PR_val.varname)>1
                error(['The parameter "',Rparams{ii_loc,jj_loc},'" can vary in at most one variable.'])
            elseif ~isempty(PR_val.varname) && ~all(ismember(PR_val.varname,Pop.var1(2).varname))
                var_indx_val = [0,1,0,0];       % The variable name can only correspond to the y-variable
            end
        else    % The parameter may depend on both x and y
            if length(PR_val.varname)>2
                error(['The parameter "',Rparams{ii_loc,jj_loc},'" can vary in at most two variables.'])
            elseif ~isempty(PR_val.varname) && ~all(ismember(PR_val.varname,Pop.var1.varname))
                var_indx_val = repmat([1,1,0,0],[length(PR_val.varname),1]);   % The variable name can correspond to both the x- and y-variable
            end
        end
        
        if ~isempty(var_indx_val)
            % % Update the variable information
            [varname_unique_PR,vindx_PR,varnum_PR] = unique(PR_val.varname);
            varnum_PR_new = varnum_PR;
            for k=1:size(varnum_PR(vindx_PR),1)
                vnum = varnum_PR(vindx_PR(k));
                % Check which (if any) variable name in varname_unique is
                % associated to the variable PR_varname_unique(vnum).
                PR_varmatch = ismember(varname_unique_Pop,varname_unique_PR(vnum));
                if any(PR_varmatch)
                    % If the variable in PR_varname_unique(vnum) already appears in
                    % in the operator, adjust the variable number of the variable
                    % in the parameter to match that in the operator.
                    varnum_PR_new(varnum_PR==vnum) = varnum_unique_Pop(PR_varmatch);
                else
                    % Otherwise, adjust the parameter varnum to assure it is
                    % distinct from the varnums of the variables in the
                    % operator.
                    varnum_PR_new(varnum_PR==vnum) = vnum + max_varnum_Pop;
                end
            end
            varnum_PR = varnum_PR_new;
            varnum_unique_PR = varnum_PR(vindx_PR);
            
            % Update the list of unique variable names and associated varnums
            varname_unique_Pop = [varname_unique_Pop; varname_unique_PR(~ismember(varnum_unique_PR,varnum_unique_Pop))];
            varnum_unique_Pop = [varnum_unique_Pop; varnum_unique_PR(~ismember(varnum_unique_PR,varnum_unique_Pop))];
            
            % Update the vartab
            vartab_PR = [varnum_PR, var_indx_val];
            ismember_valvar = ismember(vartab_PR,vartab_Pop,'rows');
            vartab_Pop = [vartab_Pop; vartab_PR(~ismember_valvar,:)];
            adjust_vartab = any(~ismember_valvar);
        end
    end
    % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % %
    % % Check that the dimensions make sense
    if all(size(PR_val)==dim_ij)
        % If all dimensions agree, just assign the parameter
        out = Pop;
        out.(Rparams{ii_loc,jj_loc}) = PR_val;
    else
        % If some dimension does not agree, check if this dimension
        % is flexible, i.e. corresponds to only 0-parameters in Pop
        nr = size(PR_val,1);           nc = size(PR_val,2);
        Pdim_new = Pdim;
        Pdim_new(isnan(Pdim_new)) = 0;
        Pdim_new(ii_loc,1) = nr;    Pdim_new(jj_loc,2) = nc;
        dim_match = Pdim_new==Pdim;
        
        if all(dim_match(~isflex_Pdim))
            % If the incorrect dimensions are flexible, expand/compress
            % 0-parameters to the new dimensions, and assign the new
            % parameter
            out = opvar2d([],Pdim_new,Pop.I,Pop.var1,Pop.var2);
            for k=1:numel(Rparams)
                [ii,jj] = ind2sub(size(Rparams),k);
                if (ii~=ii_loc || jj~=jj_loc) && ~iszero_param(ii,jj)
                    out.(Rparams{k}) = Pop.(Rparams{k});
                elseif ii==ii_loc && jj==jj_loc
                    out.(Rparams{k}) = PR_val;
                end
            end
        else
            error(['Unable to perform assignment because the size of the proposed parameter is ',num2str(nr),'-by-',num2str(nc),' but should be ',num2str(dim_ij(1)),'-by-',num2str(dim_ij(2)),'.'])
        end
    end
    % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % %
    
elseif isa(PR_val,'cell')
    % Second case: PR is a cell of parameters
    if ~isa(Pop.(Rparams{ii_loc,jj_loc}),'cell') || ~numel(Pop.(Rparams{ii_loc,jj_loc}))==numel(PR_val)
        error(['Parameters "',Rparams{ii_loc,jj_loc},'" should be specified as a cell of size ',num2str(size(Pop.Rparams{ii_loc,jj_loc},1))','-by-',num2str(size(Pop.Rparams{ii_loc,jj_loc},2)),'.'])
    end
    PR_val = reshape(PR_val,size(Pop.(Rparams{ii_loc,jj_loc})));
    nr = dim_ij(1);     nc = dim_ij(2);
    
    % % Check if all parameters have the same dimensions, and
    % % expand/compress zero arrays if necessary.
    % % Also check if all parameters have the same variables, and the
    % % number of variables makes sense.
    nonzero_idx = numel(PR_val)+1;
    var_list_PR = cell(0,1);           % List of varnames that appear in the different parameters
    var_indx_PR = zeros(0,4);   % Table indicating whether each varname may correspond to variable [x, y, tt, nu]
    for l=1:numel(PR_val)
        PRR = PR_val{l};
        if isa(PRR,'dpvar')
            isdopvar = true;
        else
            isdopvar = false;
        end
        
        % Collect the names of the variables that appear in the
        % parameter, keeping track of which variables these names
        % might belong too
        if (isa(PRR,'polynomial') || isa(PRR,'dpvar')) && ~isempty(PRR.varname)
            nPRR_vars = length(PRR.varname);
            if ii_loc==2 && jj_loc==2   % The parameters Pop.Rxx{l} can vary only in x and tt
                if l==1
                    if nPRR_vars>1
                        error(['The first parameter of P.',Rparams{ii_loc,jj_loc},' cannot vary in more than 1 variable.']);
                    elseif any(ismember(PRR.varname,[Pop.var1(2).varname;Pop.var2.varname]))
                        error(['The first parameter of P.',Rparams{ii_loc,jj_loc},' cannot vary in a dummy variable along the x-direction. Please adjust the variable names of the operator or of the proposed parameter.']);
                    elseif ~all(ismember(PRR.varname,Pop.var1(1).varname))
                        var_indx_PR = [var_indx_PR; repmat([1,0],[nPRR_vars,1])];
                    end
                else
                    if nPRR_vars>2
                        error(['The second and third parameters of P.',Rparams{ii_loc,jj_loc},' cannot vary in more than 2 variables.']);
                    elseif any(ismember(PRR.varname,[Pop.var1(2).varname;Pop.var2(2).varname]))
                        error(['The parameters P.',Rparams{ii_loc,jj_loc},' cannot vary along the y-direction. Please adjust the variable names of the operator or of the proposed parameter.']);
                    elseif ~all(ismember(PRR.varname,[Pop.var1(1).varname;Pop.var2(1).varname]))
                        var_indx_PR = [var_indx_PR; repmat([1,1],[nPRR_vars,1])];
                    end
                end
                PR_var_col_include = logical([1,0,1,0]);    % Which columns in the full var_indx are PR_var_indx associated to?
            elseif ii_loc==3 && jj_loc==3   % The parameters Pop.Ryy{l} can vary only in y and nu
                if l==1
                    if nPRR_vars>1
                        error(['The first parameter of P.',Rparams{ii_loc,jj_loc},' cannot vary in more than 1 variable.']);
                    elseif any(ismember(PRR.varname,[Pop.var1(1).varname;Pop.var2.varname]))
                        error(['The first parameter of P.',Rparams{ii_loc,jj_loc},' cannot vary in a dummy variable along the y-direction. Please adjust the variable names of the operator or of the proposed parameter.']);
                    elseif ~all(ismember(PRR.varname,Pop.var1(2).varname))
                        var_indx_PR = [var_indx_PR; repmat([1,0],[nPRR_vars,1])];
                    end
                else
                    if nPRR_vars>2
                        error(['The second and third parameters of P.',Rparams{ii_loc,jj_loc},' cannot vary in more than 2 variables.']);
                    elseif any(ismember(PRR.varname,[Pop.var1(1).varname;Pop.var2(1).varname]))
                        error(['The parameters P.',Rparams{ii_loc,jj_loc},' cannot vary along the x-direction. Please adjust the variable names of the operator or of the proposed parameter.']);
                    elseif ~all(ismember(PRR.varname,[Pop.var1(2).varname;Pop.var2(2).varname]))
                        var_indx_PR = [var_indx_PR; repmat([1,1],[nPRR_vars,1])];
                    end
                end
                PR_var_col_include = logical([0,1,0,1]);    % Which columns in the full var_indx are PR_var_indx associated to?
            elseif ii_loc==2 || jj_loc==2       % The parameters Pop.Rx2{l} and Pop.R2x{l} can vary only in x, y and tt
                if l==1
                    if nPRR_vars>2
                        error(['The first parameter of P.',Rparams{ii_loc,jj_loc},' cannot vary in more than 2 variable.']);
                    elseif any(ismember(PRR.varname,Pop.var2.varname))
                        error(['The first parameter of P.',Rparams{ii_loc,jj_loc},' cannot vary in a dummy variable. Please adjust the variable names of the operator or of the proposed parameter.']);
                    elseif ~all(ismember(PRR.varname,Pop.var1.varname))
                        var_indx_PR = [var_indx_PR; repmat([1,1,0],[nPRR_vars,1])];
                    end
                else
                    if nPRR_vars>3
                        error(['The second and third parameters of P.',Rparams{ii_loc,jj_loc},' cannot vary in more than 3 variables.']);
                    elseif any(ismember(PRR.varname,Pop.var2(2).varname))
                        error(['The parameters P.',Rparams{ii_loc,jj_loc},'{2} and P.',Rparams{ii_loc,jj_loc},'{3} cannot vary in a dummy variable along the y-direction. Please adjust the variable names of the operator or of the proposed parameter.']);                    
                    elseif ~all(ismember(PRR.varname,[Pop.var1.varname;Pop.var2(1).varname]))
                        var_indx_PR = [var_indx_PR; repmat([1,1,1],[nPRR_vars,1])];
                    end
                end
                PR_var_col_include = logical([1,1,1,0]);    % Which columns in the full var_indx are PR_var_indx associated to?
            elseif ii_loc==3 || jj_loc==3       % The parameters Pop.Ry2{l} and Pop.R2y{l} can vary only in x, y and nu
                if l==1
                    if nPRR_vars>2
                        error(['The first parameter of P.',Rparams{ii_loc,jj_loc},' cannot vary in more than 2 variable.']);
                    elseif any(ismember(PRR.varname,Pop.var2.varname))
                        error(['The first parameter of P.',Rparams{ii_loc,jj_loc},' cannot vary in a dummy variable. Please adjust the variable names of the operator or of the proposed parameter.']);
                    elseif ~all(ismember(PRR.varname,Pop.var1.varname))
                        var_indx_PR = [var_indx_PR; repmat([1,1,0],[nPRR_vars,1])];
                    end
                else
                    if nPRR_vars>3
                        error(['The second and third parameters of P.',Rparams{ii_loc,jj_loc},' cannot vary in more than 3 variables.']);
                    elseif any(ismember(PRR.varname,Pop.var2(1).varname))
                        error(['The parameters P.',Rparams{ii_loc,jj_loc},'{2} and P.',Rparams{ii_loc,jj_loc},'{3} cannot vary in a dummy variable along the x-direction. Please adjust the variable names of the operator or of the proposed parameter.']);
                    elseif ~all(ismember(PRR.varname,[Pop.var1.varname;Pop.var2(2).varname]))
                        var_indx_PR = [var_indx_PR; repmat([1,1,1],[nPRR_vars,1])];
                    end
                end
                PR_var_col_include = logical([1,1,0,1]);    % Which columns in the full var_indx are PR_var_indx associated to?
            else                                % The parameters Pop.R22{l} can vary in x, y, tt and nu
                [ii_l,jj_l] = ind2sub(size(PR_val),l);
                if ii_l==1 && jj_l==1
                    if nPRR_vars>2
                        error(['The parameter P.',Rparams{ii_loc,jj_loc},'{',num2str(ii_l),',',num2str(jj_l),'} cannot vary in more than 2 variables.']);
                    elseif any(ismember(PRR.varname,Pop.var2.varname))
                        error(['The first parameter of P.',Rparams{ii_loc,jj_loc},' cannot vary in a dummy variable. Please adjust the variable names of the operator or of the proposed parameter.']);
                    elseif ~all(ismember(PRR.varname,Pop.var1.varname))
                        var_indx_PR = [var_indx_PR; repmat([1,1,0,0],[nPRR_vars,1])];
                    end
                elseif ii_l==1 || jj_l==1
                    if nPRR_vars>3
                        error(['The parameter P.',Rparams{ii_loc,jj_loc},'{',num2str(ii_l),',',num2str(jj_l),'} cannot vary in more than 3 variables.']);
                    elseif jj_l==1 && any(ismember(PRR.varname,Pop.var2(2).varname))
                        error(['The parameters P.',Rparams{ii_loc,jj_loc},'{2,1} and P.',Rparams{ii_loc,jj_loc},'{3,1} cannot vary in a dummy variable along the y-direction. Please adjust the variable names of the operator or of the proposed parameter.']);
                    elseif jj_l==1 && ~all(ismember(PRR.varname,[Pop.var1.varname;Pop.var2(1).varname]))
                        var_indx_PR = [var_indx_PR; repmat([1,1,1,0],[nPRR_vars,1])];
                    elseif ii_l==1 && any(ismember(PRR.varname,Pop.var2(1).varname))
                        error(['The parameters P.',Rparams{ii_loc,jj_loc},'{1,2} and P.',Rparams{ii_loc,jj_loc},'{1,3} cannot vary in a dummy variable along the x-direction. Please adjust the variable names of the operator or of the proposed parameter.']);
                    elseif ii_l==1 && ~all(ismember(PRR.varname,[Pop.var1.varname;Pop.var2(2).varname]))
                        var_indx_PR = [var_indx_PR; repmat([1,1,0,1],[nPRR_vars,1])];
                    end
                else
                    if nPRR_vars>4
                        error(['The parameter P.',Rparams{ii_loc,jj_loc},'{',num2str(ii_l),',',num2str(jj_l),'} cannot vary in more than 4 variables.']);
                    elseif ~all(ismember(PRR.varname,[Pop.var1.varname;Pop.var2.varname]))
                        var_indx_PR = [var_indx_PR; repmat([1,1,1,1],[nPRR_vars,1])];
                    end
                end
                PR_var_col_include = logical([1,1,1,1]);    % Which columns in the full var_indx_Pop are var_indx_PR associated to?
            end
            if ~isempty(var_indx_PR)
                var_list_PR = [var_list_PR; PRR.varname];   % A full list of variable names that appear in the new parameter
            end
        end
        
        % Check the dimensions of each parameter PR{l}. If the dimensions
        % do not match, but the parameter is all zeroes, expand/compress
        % the parameter to a 0-parameter of appropriate dimensions
        if isempty(PRR) || ...
                (isa(PRR,'double') && all(all(PRR==0))) || ...
                (isa(PRR,'polynomial') && all(all(isequal(PRR,0)))) || ...
                (isa(PRR,'dpvar') && all(all(PRR.C==0)))
            if nonzero_idx<=numel(PR_val)
                PR_val{l} = polynomial(zeros(nr,nc));
            end
        elseif (isa(PRR,'double') || isa(PRR,'polynomial') || isa(PRR,'dpvar'))
            if nonzero_idx==numel(PR_val)+1
                nonzero_idx = min(nonzero_idx,l);
                nr = size(PRR,1);   nc = size(PRR,2);
            elseif any(size(PRR)~=[nr,nc])
                error('Dimensions of the proposed parameters do not match');
            end
        else
            error('Parameters of opvar2d object must be specified as ''double'' or ''polynomial'' class objects');
        end
    end
    % Expand/compress zero arrays that we missed
    for l=1:nonzero_idx-1
        PR_val{l} = polynomial(zeros(nr,nc));
    end
    
    % % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - % %
    % % Having checked the dimensions and variables of each parameter
    % % PR{l}, add the new variable information to vartab.
    if ~isempty(var_list_PR)
        % Establish a unique set of variable names in PR.
        [varname_unique_PR,vindx_PR,varnum_PR] = unique(var_list_PR);
        varnum_PR_new = varnum_PR;
        for k=1:size(varnum_PR(vindx_PR),1)
            vnum = varnum_PR(vindx_PR(k));
            % Check which (if any) variable name in varname_unique is
            % associated to the variable PR_varname_unique(vnum).
            PR_varmatch = ismember(varname_unique_Pop,varname_unique_PR(vnum));
            if any(PR_varmatch)
                % If the variable in PR_varname_unique(vnum) already appears in
                % in the operator, adjust the variable number of the variable
                % in the parameter to match that in the operator.
                varnum_PR_new(varnum_PR==vnum) = varnum_unique_Pop(PR_varmatch);
            else
                % Otherwise, adjust the parameter varnum to assure it is
                % distinct from the varnums of the variables in the
                % operator.
                varnum_PR_new(varnum_PR==vnum) = vnum + max_varnum_Pop;
            end
        end
        varnum_PR = varnum_PR_new;
        varnum_unique_PR = varnum_PR(vindx_PR);
        
        % Update the list of unique variable names and associated varnums
        varname_unique_Pop = [varname_unique_Pop; varname_unique_PR(~ismember(varnum_unique_PR,varnum_unique_Pop))];
        varnum_unique_Pop = [varnum_unique_Pop; varnum_unique_PR(~ismember(varnum_unique_PR,varnum_unique_Pop))];
        
        % Build a vartable for PR:
        % The first column of vartab indicates the considered variable name in
        % varname_unique. The remaining columns are logical indices indicating
        % which spatial variables this varname may be associated to
        vartab_PR = [varnum_PR, var_indx_PR];
        [vartab_PR,~] = unique(vartab_PR,'rows');
        
        % Include contributions from variables that do not appear in
        % PR_vartab (but may appear in Pop) as zero.
        vartab_PR_tmp = zeros(size(vartab_PR,1),nvars_Pop+1);
        vartab_PR_tmp(:,[true,PR_var_col_include]) = vartab_PR;
        
        % Add the new parameter variable information to the current variable
        % information on the operator.
        ismember_valvar = ismember(vartab_PR_tmp,vartab_Pop,'rows');
        vartab_Pop = [vartab_Pop; vartab_PR_tmp(~ismember_valvar,:)];
        adjust_vartab = any(~ismember_valvar);
    end
    
    
    % % Finally, check that the dimensions of the parameter match those of
    % % the operator.
    if all([nr,nc]==dim_ij)
        out = Pop;
        out.(Rparams{ii_loc,jj_loc}) = PR_val;
    else
        % If some dimension does not agree, check if this dimension
        % is flexible, i.e. corresponds to only 0-parameters in Pop
        Pdim_new = Pdim;
        Pdim_new(isnan(Pdim_new)) = 0;
        Pdim_new(ii_loc,1) = nr;    Pdim_new(jj_loc,2) = nc;
        dim_match = Pdim_new==Pdim;
        
        if all(dim_match(~isflex_Pdim))
            % If the incorrect dimensions are flexible, expand/compress
            % 0-parameters to the new dimensions, and assign the new
            % parameter
            out = opvar2d([],Pdim_new,Pop.I,Pop.var1,Pop.var2);
            for k=1:numel(Rparams)
                [ii,jj] = ind2sub(size(Rparams),k);
                if (ii~=ii_loc || jj~=jj_loc) && ~iszero_param(ii,jj)
                    out.(Rparams{k}) = Pop.(Rparams{k});
                elseif ii==ii_loc && jj==jj_loc
                    out.(Rparams{k}) = PR_val;
                end
            end
        else
            error(['Unable to perform assignment because the size of the proposed parameter is ',num2str(nr),'-by-',num2str(nc),' but should be ',num2str(dim_ij(1)),'-by-',num2str(dim_ij(2)),'.'])
        end
    end    
else
    error('Parameters of opvar2d object must be specified as ''double'' or ''polynomial'' class objects, or cells (where applicable).');
end

% % If new variable information is available, run through the vartable
% % again, aiming to establish varnums (and thus variable names) associated
% % to individual variables (rather than combinations).
if adjust_vartab
    nrows_Pop_vartab_old = 0;
    nrows_Pop_vartab_new = size(vartab_Pop,1);
    for rep = 1:2^nvars_Pop-1
        % If the vartab appears not to change anymore, stop repeating
        if nrows_Pop_vartab_new==nrows_Pop_vartab_old
            break
        else
            nrows_Pop_vartab_old = size(vartab_Pop,1);
        end
        % Loop over each possible combination of variables
        for k=1:2^nvars_Pop-1
            % Establish a combination k of variables, and determine which
            % rows in vartab are associated to this combination
            var_cols_k = logical(str2num(dec2bin(k,nvars_Pop)')');
            contains_all_var_k = all(vartab_Pop(:,[false,var_cols_k]),2);
            excludes_any_var_notk = ~any(vartab_Pop(:,[false,~var_cols_k]),2);
            contains_only_var_k = contains_all_var_k & excludes_any_var_notk;
            
            % Extract variable indices associated to the variable combination k
            varnum_k = vartab_Pop(contains_only_var_k,1);
            varnum_k = unique(varnum_k);
            
            % Check that the number of variable names does not exceed the maximum
            % number of variables that may appear in the combination k
            if length(varnum_k)>sum(var_cols_k)
                error(['The number of spatial variables in the proposed operator exceeds the allowed amount of variables. '...
                       'Please make sure the variables in the proposed parameter match those already present in the operator, and that the total number of distinct variables does not exceed 4.']);
            elseif isempty(varnum_k)    % If we have no variable index, just move on to the next combination
                continue
            end
            
            % Make sure the current variable names also appear in any other
            % combination kl of variables that includes the combination k
            for l=1:2^(sum(~var_cols_k))-1        % Loop over variable combinations that exclude combination k
                var_cols_l = logical(str2num(dec2bin(l,sum(~var_cols_k))')');
                var_cols_kl = var_cols_k;
                var_cols_kl(~var_cols_k) = var_cols_l;                
                
                % Only add rows to the table that do not currently appear
                vartab_l = [varnum_k, repmat(var_cols_kl,[length(varnum_k),1])];
                ismember_var_l = ismember(vartab_l,vartab_Pop,'rows');
                vartab_Pop = [vartab_Pop; vartab_l(~ismember_var_l,:)];
            end
            vartab_Pop = unique(vartab_Pop,'rows');
            
            % Try to establish other combinations of variables l by excluding the
            % current combination k from combinations kl
            for l=2^(sum(~var_cols_k))-1 :-1:1    % Loop over variable combinations that exclude combination k
                var_cols_l = logical(str2num(dec2bin(l,sum(~var_cols_k))')');
                var_cols_kl = var_cols_k;
                var_cols_kl(~var_cols_k) = var_cols_l;
                
                % If we have n variable names for a combination kl of n variables,
                % we can get rid of the known variables k in the combination to
                % establish variable names for the unknown variables l
                contains_all_var_kl = all(vartab_Pop(:,[false,var_cols_kl]),2);
                excludes_any_var_notkl = ~any(vartab_Pop(:,[false,~var_cols_kl]),2);
                contains_only_var_kl = contains_all_var_kl & excludes_any_var_notkl;
                vartab_kl = vartab_Pop(contains_only_var_kl,:);
                if size(vartab_kl,1)==sum(var_cols_kl) && length(varnum_k)==sum(var_cols_k)
                    isnot_kvar = ~ismember(vartab_kl(:,1),varnum_k);
                    if sum(isnot_kvar)~=(size(vartab_kl,1)-length(varnum_k))
                        error('Variables in the proposed operator cannot act as both variables along the x-direction and the y-direction, nor as both primary and dummy variables.')
                    end
                    vartab_l = vartab_kl(isnot_kvar,:);
                    vartab_l(:,[false,var_cols_k]) = false;
                    % Only add rows to the table that do not currently appear
                    ismember_var_l = ismember(vartab_l,vartab_Pop,'rows');
                    vartab_Pop = [vartab_Pop; vartab_l(~ismember_var_l,:)];
                end
            end
        end
        nrows_Pop_vartab_new = size(vartab_Pop,1);
    end
    
    % % Using the new variable information, (re)set the variables of the
    % % output
    % Set the x-variable
    varnum_var1_1 = unique(vartab_Pop(vartab_Pop(:,2)&~any(vartab_Pop(:,[3,4,5]),2),1));
    if isempty(varnum_var1_1)
        var1_1 = Pop.var1(1);   % If we have insufficient information, use the old variable
    elseif length(varnum_var1_1)>1
        error('The x-variable in the proposed parameter does not match the current x-variable of the operator.')
    else
        varname_var1_1 = varname_unique_Pop(varnum_unique_Pop==varnum_var1_1);
        var1_1 = polynomial(varname_var1_1);
    end
    
    % Set the y-variable
    varnum_var1_2 = unique(vartab_Pop(vartab_Pop(:,3)&~any(vartab_Pop(:,[2,4,5]),2),1));
    if isempty(varnum_var1_2)
        var1_2 = Pop.var1(2);
    elseif length(varnum_var1_2)>1
        error('The y-variable in the proposed parameter does not match the current y-variable of the operator.')
    else
        varname_var1_2 = varname_unique_Pop(varnum_unique_Pop==varnum_var1_2);
        var1_2 = polynomial(varname_var1_2);
    end
    
    % Set the dummy x-variable
    varnum_var2_1 = unique(vartab_Pop(vartab_Pop(:,4)&~any(vartab_Pop(:,[2,3,5]),2),1));
    if isempty(varnum_var2_1)
        var2_1 = Pop.var2(1);
    elseif length(varnum_var2_1)>1
        error('The dummy-x-variable in the proposed parameter does not match the current dummy-x-variable of the operator.')
    else
        varname_var2_1 = varname_unique_Pop(varnum_unique_Pop==varnum_var2_1);
        var2_1 = polynomial(varname_var2_1);
    end
    
    % Set the dummy y-variable
    varnum_var2_2 = unique(vartab_Pop(vartab_Pop(:,5)&~any(vartab_Pop(:,[2,3,4]),2),1));
    if isempty(varnum_var2_2)
        var2_2 = Pop.var2(2);
    elseif length(varnum_var2_2)>1
        error('The dummy-y-variable in the proposed parameter does not match the current dummy-y-variable of the operator.')
    else
        varname_var2_2 = varname_unique_Pop(varnum_unique_Pop==varnum_var2_2);
        var2_2 = polynomial(varname_var2_2);
    end
    
    var1 = [var1_1; var1_2];
    var2 = [var2_1; var2_2];
    % Check that the new set of variables makes sense.
    if length(unique([var1.varname;var2.varname]))~=(length(var1)+length(var2))
        error(['The variables in the proposed parameters conflict with the currently set variables of the operator.',...
                ' Please adjust the variables in the parameter, or the variables "Pop.var1" and "Pop.var2" of the operator.'])
    end
    
    % Finally, some checks:
    if isempty(varnum_var1_1) && isempty(varnum_var1_2)
        varnum_var1 = unique(vartab_Pop(all(vartab_Pop(:,[2,3]),2)&~any(vartab_Pop(:,[4,5]),2),1));
        if ~isempty(varnum_var1)
            varname_var1 = varname_unique_Pop(ismember(varnum_unique_Pop,varnum_var1));
            if length(varnum_var1)>2
                error('The proposed operator contains more (primary) spatial variables than allowed.')
            elseif length(varnum_var1)==2 && ~all(ismember(varname_var1,var1.varname))
                % If one of the variable names is already associated to a
                % variable, link the other variable name with the other
                % variable.
                var_match = ismember(varname_var1,var1.varname);
                if any(ismember(varname_var1(var_match),var1(1).varname))
                    var1(2) = polynomial(varname_var1(~var_match));
                elseif any(ismember(varname_var1(var_match),var1(2).varname))
                    var1(1) = polynomial(varname_var1(~var_match));
                else
                    warning(['It is unclear which variable in your operator correspond to which spatial dimension. Please explicitly specify the primary variables in your operator as e.g. "P.var1 = [x;y]" and the dummy variables as e.g. "P.var2 = [theta;nu]".'])
                end
            end
        end
    end
    if isempty(varnum_var1_1) && isempty(varnum_var2_1)
        varnum_var_1 = unique(vartab_Pop(all(vartab_Pop(:,[2,4]),2)&~any(vartab_Pop(:,[3,5]),2),1));
        if ~isempty(varnum_var_1)
            varname_var_1 = varname_unique_Pop(ismember(varnum_unique_Pop,varnum_var_1));
            if length(varnum_var_1)>2
                error('The proposed operator contains more spatial variables along the x-direction than allowed.')
            elseif length(varnum_var_1)==2 && ~all(ismember(varname_var_1,[var1(1).varname;var2(1).varname]))
                % If one of the variable names is already associated to a
                % variable, link the other variable name with the other
                % variable.
                var_match = ismember(varname_var_1,[var1(1).varname;var2(1).varname]);
                if ismember(varname_var_1(var_match),var1(1).varname)
                    var2(1) = polynomial(varname_var_1(~var_match));
                elseif ismember(varname_var_1(var_match),var2(1).varname)
                    var1(1) = polynomial(varname_var_1(~var_match));
                else
                    warning(['It is unclear which variable along the x-direction is the primary and dummy variable. Please explicitly specify the primary variables in your operator as e.g. "P.var1 = [x;y]" and the dummy variables as e.g. "P.var2 = [theta;nu]".'])
                end
            end
        end
    end
    if isempty(varnum_var1_2) && isempty(varnum_var2_2)
        varnum_var_2 = unique(vartab_Pop(all(vartab_Pop(:,[3,5]),2)&~any(vartab_Pop(:,[2,4]),2),1));
        if ~isempty(varnum_var_2)
            varname_var_2 = varname_unique_Pop(ismember(varnum_unique_Pop,varnum_var_2));
            if length(varnum_var_2)>2
                error('The proposed operator contains more spatial variables along the y-direction than allowed.')
            elseif length(varnum_var_2)==2 && ~all(ismember(varname_var_2,[var1(2).varname;var2(2).varname]))
                % If one of the variable names is already associated to a
                % variable, link the other variable name with the other
                % variable.
                var_match = ismember(varname_var_2,[var1(2).varname;var2(2).varname]);
                if any(ismember(varname_var_2(var_match),var1(2).varname))
                    var2(2) = polynomial(varname_var_2(~var_match));
                elseif any(ismember(varname_var_2(var_match),var2(2).varname))
                    var1(2) = polynomial(varname_var_2(~var_match));
                else
                    warning(['It is unclear which variable along the y-direction is the primary and dummy variable. Please explicitly specify the primary variables in your operator as e.g. "P.var1 = [x;y]" and the dummy variables as e.g. "P.var2 = [theta;nu]".'])
                end
            end
        end
    end
else
    % If we have no new variable information, just use the current
    % variables.
    var1 = Pop.var1;
    var2 = Pop.var2;
end

% Set the variables of the output
out = set_var1(out,var1,'nosubs');
out = set_var2(out,var2,'nosubs');

% If any parameter is a dpvar class object, convert the operator to a
% dopvar2d class object
if isdopvar
    out = opvar2dopvar2d(out);
end

end