function Pop_out = subs_vars_op(Pop_in,vars_old,vars_new)
% POP_OUT = SUBS_VARS_OP(POP_IN,VARS_OLD,VARS_NEW) takes an operator
% "Pop_in" specified as 'opvar' or 'opvar2d' object, and replaces variables
% "vars_old" appearing in the operator with variables "vars_new".
%
% INPUT
% - Pop_in:     'opvar', 'dopvar', 'opvar2d', or 'dopvar2d' object
%               specifying an operator (variable).
% - vars_old:   nx1 array of type 'polynomial', with each element
%               specifying a single polynomial variable that appears in the
%               parameters defining "Pop_in". Each variable can
%               correspond to a spatial variable (Pop_in.var1) or dummy
%               variable (Pop_in.var2), but this is not necessary.
% - vars_new:   nx1 array of type 'polynomial', with each element
%               specifying a single polynomial variable with which to
%               replace the associated variables from "vars_old" in the
%               parameters defining "Pop_in".
%
% OUTPUT
% - Pop_out:    Object of same type as "Pop_in", specifying the same
%               operator (variable), but now with all instance of the
%               variables in "vars_old" replaced with the variables in
%               "vars_new", in all parameters. If any of the variables in
%               "vars_old" appear in 'Pop_in.var1' or 'Pop_in.var2', then
%               those instances are replaced as well.
%
% NOTES
% The input "vars_new" must correspond to individual polynomial variables;
% substitution with polynomial expressions is not supported.
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subs_vars_op
%
% Copyright (C)2024  PIETOOLS Team
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
% 11/29/2024, DJ: Initial coding;

% % % Check the inputs
% Check that the operator is properly specified.
use_2D = false;
if isa(Pop_in,'opvar2d') || isa(Pop_in,'dopvar2d')
    use_2D = true;
elseif ~isa(Pop_in,'opvar') && ~isa(Pop_in,'dopvar')
    error('Operator must be specified as object of type ''opvar'', ''dopvar'', ''opvar2d'', or ''dopvar2d''.')
end
% Check that the old variables are properly specified.
if ischar(vars_old)
    vars_old = polynomial({vars_old});
elseif iscellstr(vars_old)
    vars_old = polynomial(vars_old);
elseif ~isa(vars_old,'polynomial') || ~ispvar(vars_old)
    error('Variables to substitute should be specified as array of type ''polynomial'', with each element corresponding to a single variable.')
end
% Check that the new variables are properly specified.
if ischar(vars_new)
    vars_new = polynomial({vars_new});
elseif iscellstr(vars_new)
    vars_new = polynomial(vars_new);
elseif ~isa(vars_new,'polynomial') || ~ispvar(vars_new)
    error('Variables to substitute should be specified as array of type ''polynomial'', with each element corresponding to a single variable.')
end
if any(size(vars_old)~=size(vars_new))
    error('Number of original variables should match number of new variables.')
end
vars_old = vars_old(:);
vars_new = vars_new(:);


% % % Perform the actual substitution.
Pop_out = Pop_in;
if ~use_2D
    % 1D case --> parameters are stored in the following fields:
    fnames1 = {'P','Q1','Q2'};
    fnames2 = {'R0','R1','R2'};
    for jj=1:numel(vars_old)
        % Extract the old and new variables.
        varname_old = vars_old(jj).varname;
        varname_new = vars_new(jj).varname;
        if strcmp(varname_old{1},varname_new{1})
            continue
        end
        % Replace primary spatial variable if applicable.
        if strcmp(Pop_out.var1.varname{1},varname_old{1})
            Pop_out.var1 = vars_new(jj);
        end
        % Replace dummy variable if applicable.
        if strcmp(Pop_out.var2.varname{1},varname_old{1})
            Pop_out.var2 = vars_new(jj);
        end
        % Substitute variable in parameters
        for fname = fnames1
            Rjj = Pop_out.(fname{1});
            if isa(Rjj,'polynomial') || isa(Rjj,'dpvar')
                var_idx = ismember(Rjj.varname,varname_old);
                Rjj.varname(var_idx) = varname_new;
            end
            Pop_out.(fname{1}) = Rjj;
        end
        for fname = fnames2
            % Paramerers from 'fnames2' are stored under Pop_out.R
            Rjj = Pop_out.R.(fname{1});
            if isa(Rjj,'polynomial') || isa(Rjj,'dpvar')
                var_idx = ismember(Rjj.varname,varname_old);
                Rjj.varname(var_idx) = varname_new;
            end
            Pop_out.R.(fname{1}) = Rjj;
        end
    end
else
    % 2D case --> parameters are stored in the following fields:
    fnames1 = {'R00','R0x','R0y','R02','Rx0','Rxy','Ry0','Ryx','R20'};
    fnames2 = {'Rxx','Rx2','Ryy','Ry2','R2x','R2y','R22'};
    for jj=1:numel(vars_old)
        % Extract the old and new variables.
        varname_old = vars_old(jj).varname;
        varname_new = vars_new(jj).varname;
        if strcmp(varname_old{1},varname_new{1})
            continue
        end
        % Replace primary spatial variable if applicable.
        var1_idx = ismember(Pop_out.var1.varname,varname_old);
        Pop_out.var1.varname(var1_idx) = varname_new;
        % Replace dummy variable if applicable.
        var2_idx = ismember(Pop_out.var2.varname,varname_old);
        Pop_out.var2.varname(var2_idx) = varname_new;
        % Substitute variable in parameters
        for fname = fnames1
            Rjj = Pop_out.(fname{1});
            if isa(Rjj,'polynomial') || isa(Rjj,'dpvar')
                var_idx = ismember(Rjj.varname,varname_old);
                Rjj.varname(var_idx) = varname_new;
            end
            Pop_out.(fname{1}) = Rjj;
        end
        for fname = fnames2
            % Parameters from 'fnames2' are stored in cells
            Rjj = Pop_out.(fname{1});
            for kk=1:numel(Rjj)
                Rjj_kk = Rjj{kk};
                if isa(Rjj_kk,'polynomial') || isa(Rjj_kk,'dpvar')
                    var_idx = ismember(Rjj_kk.varname,varname_old);
                    Rjj_kk.varname(var_idx) = varname_new;
                end
                Rjj{kk} = Rjj_kk;
            end
            Pop_out.(fname{1}) = Rjj;
        end
    end
end

end