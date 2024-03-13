function Pop_2d = opvar2opvar2d(Pop_1d,dom_new,var1_new,var2_new,use_space_idx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pop_2d = opvar2opvar2d(Pop_1d,dom_new,var1_new,var2_new,use_space_idx) 
% converts an opvar object to a corresponding dopvar2d object, with 2d
% domain and variables as determined by the additional inputs.
% 
% INPUT
% - Pop_1d:     opvar object, representing an operator mapping
%               R^n0 x L2^n1[s1] --> R^m0 x L2^m1[s1].
% - dom_new:    1x2 array specifying the domain of the new variable
%               s2 in the opvar2d object.
% - var1_new:   1x1 pvar (polynomial class) object specifying the new 
%               primary variable s2 in the opvar2d object.
% - var2_new:   1x1 pvar (polynomial class) object specifying the new 
%               dummy variable theta2 in the opvar2d object.
% - use_space_idx:  Scalar 1 or 2, indicating whether the new spatial
%                   variable s2 will be the first variable, or the second
%                   variable, i.e. if the opvar2d object will map
%                       R^n0 x L2^n1[s1] x L2^0[s2] x L2^0[s1,s2]
%                           --> R^m0 x L2^m1[s1] x L2^0[s2] x L2^0[s1,s2];
%                   or
%                       R^n0 x L2^0[s2] x L2^n1[s1] x L2^0[s2,s1]
%                           --> R^0m0 x L2^0[s2] x L2^m1[s1] x L2^0[s2,s1];
% 
% OUTPUT
% - Pop_2d:     opvar_2d object, representing an operator mapping
%                   R^n0 x L2^n1[s1] x L2^0[s2] x L2^0[s1,s2]
%                       --> R^m0 x L2^m1[s1] x L2^0[s2] x L2^0[s1,s2];
%               if use_space_idx==1, or
%                   R^n0 x L2^0[s2] x L2^n1[s1] x L2^0[s2,s1]
%                       --> R^0m0 x L2^0[s2] x L2^m1[s1] x L2^0[s2,s1];
%               if use_space_idx==2.
%               The new variable s2 will have a variable name as specified
%               by var1_new, and exist on the interbal defined by dom_new.
%
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu
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
% Initial coding DJ, - 08/12/2022.

% % % Process the inputs.

% Check if we're converting an opvar or dopvar object.
is_dopvar = false;
if isa(Pop_1d,'dopvar')
    is_dopvar = true;
elseif ~isa(Pop_1d,'opvar')
    error(['Input must be of class ''opvar'' or ''dopvar''.'])
end

% Check that the new domain is properly specified.
if nargin>=2
    if ~all(size(dom_new)==[1,2])
        error(['The domain of the new spatial variable should be specified as a 1x2 array']);
    elseif dom_new(1)>=dom_new(2)
        error(['The upper boundary of the domain of the new spatial variable should be strictly greater than the lower boundary.'])
    end
else
    % If no domain is provided for the new variable, assume the same domain
    % as for the current variable.
    dom_new = Pop_1d.I;
end

use_space = false(2,1);
% Check that the variables and "use_space_idx" have been properly
% specified.
if nargin==3
    % Convert variable names to polynomial variables.
    if iscellstr(var1_new)
        var1_new = polynomial(var1_new);
    end
    % Split primary and dummy variable if specified in a single array.
    if ~ispvar(var1_new) || length(var1_new)>2
        error(['The new spatial variable should be specified as a 1x2 pvar (polynomial class) array.'])
    elseif length(var1_new)==2
        var2_new = var1_new(2);
        var1_new = var1_new(1);
    end
elseif nargin==4
    % Convert variable names to polynomial variables.
    if iscellstr(var1_new)
        var1_new = polynomial(var1_new);
    end
    if iscellstr(var2_new)
        var2_new = polynomial(var2_new);
    end
    if ~ispvar(var2_new) && (~ispvar(var1_new) || length(var1_new)~=2)
        error(['The new primary and dummy variables should be specified as a 1x1 pvar (polynomial class) objects.'])
    elseif ~ispvar(var2_new) 
        % If var1_new includes both primary and dummy variable, we assume
        % var2_new is the "use_space_idx".
        use_space_idx = var2_new;
        var2_new = var1_new(2);
        var1_new = var1_new(1);
        if ~ismember(use_space_idx,[1;2;'x';'y'])
            error(['The desired function space index should be either 1 or 2'])
        elseif ismember(use_space_idx,[1;'x'])
            use_space(1) = true;
        else
            use_space(2) = true;
        end
    else
        if ~ispvar(var1_new) || length(var1_new)~=1 || length(var2_new)~=1
            error(['The new primary and dummy variables should be specified as a 1x1 pvar (polynomial class) objects.'])
        end
    end
    % Check that the new variables have been properly specified.
    if ~ispvar(var1_new) || length(var1_new)>2
        error(['The new spatial variable should be specified as a 1x2 pvar (polynomial class) array.'])
    elseif length(var1_new)==2
        % Split primary and dummy variable if specified in a single array.
        var2_new = var1_new(2);
        var1_new = var1_new(1);
    end
elseif nargin==5
    % Convert variable names to polynomial variables.
    if iscellstr(var1_new)
        var1_new = polynomial(var1_new);
    end
    if iscellstr(var2_new)
        var2_new = polynomial(var2_new);
    end
    % Check that the new variables have been properly specified.
    if ~ispvar(var1_new) || length(var1_new)~=1 || ~ispvar(var2_new) || length(var2_new)~=1
        error(['The new primary and dummy spatial variables should be specified as 1x1 pvar (polynomial class) objects.'])
    end
    % Check that the "use_space_idx" has been properly specified.
    if ~ismember(use_space_idx,[1;2;'x';'y'])
        error(['The desired function space index should be either 1 or 2'])
    elseif ismember(use_space_idx,[1;'x'])
        use_space(1) = true;
    else
        use_space(2) = true;
    end
end
% If the user has not specified whether the opvar object maps functions in
% the first or second variable of the new opvar2d object, assume the first
% variable.
if ~any(use_space)
    use_space(1) = true;
end
% If no primary variable has been specified, create a new one.
if ~exist('var1_new','var')
    if use_space(1)
        var1_new = polynomial({[Pop_1d.var1.varname{1},'_2']});
    else
        var1_new = polynomial({[Pop_1d.var1.varname{1},'_1']});
    end
end
% If no dummy variable has been specified, create a new one.
if ~exist('var2_new','var')
    if use_space(1)
        var2_new = polynomial({[Pop_1d.var2.varname{1},'_2']});
    else
        var2_new = polynomial({[Pop_1d.var2.varname{1},'_1']});
    end
end



% % % Construct the new opvar2d object.
if use_space(1)
    % The new variable is the second spatial variable.
    var1 = [Pop_1d.var1; var1_new];
    var2 = [Pop_1d.var2; var2_new];
    dom = [Pop_1d.I; dom_new];
    dim = [Pop_1d.dim; [0,0;0,0]];
    
    if is_dopvar
        Pop_2d = dopvar2d([],dim,dom,var1,var2);
    else
        Pop_2d = opvar2d([],dim,dom,var1,var2);
    end
    Pop_2d.R00 = Pop_1d.P;
    Pop_2d.R0x = Pop_1d.Q1;
    Pop_2d.Rx0 = Pop_1d.Q2;
    Pop_2d.Rxx{1} = Pop_1d.R.R0;
    Pop_2d.Rxx{2} = Pop_1d.R.R1;
    Pop_2d.Rxx{3} = Pop_1d.R.R2;
    
else
    % The new variable is the first spatial variable.
    var1 = [var1_new; Pop_1d.var1];
    var2 = [var2_new; Pop_1d.var2];
    dom = [dom_new; Pop_1d.I];
    dim = [Pop_1d.dim(1,:); [0,0]; Pop_1d.dim(2,:); [0,0]];

    if is_dopvar
        Pop_2d = dopvar2d([],dim,dom,var1,var2);
    else
        Pop_2d = opvar2d([],dim,dom,var1,var2);
    end
    Pop_2d.R00 = Pop_1d.P;
    Pop_2d.R0y = Pop_1d.Q1;
    Pop_2d.Ry0 = Pop_1d.Q2;
    Pop_2d.Ryy{1} = Pop_1d.R.R0;
    Pop_2d.Ryy{2} = Pop_1d.R.R1;
    Pop_2d.Ryy{3} = Pop_1d.R.R2;
end

end