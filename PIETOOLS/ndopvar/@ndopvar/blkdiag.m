function Fdiag = blkdiag(varargin)
% F = blkdiag(varargin) generates a block diagonal ndopvar object using
% varargs as the blocks on the diagonal
% 
% INPUTS
% - varargin:   'ndopvar' class objects.
% 
% OUTPUTS
% - Fdiag:      'ndopvar' object with each of the parameters definined as the
%               block diagonal concatenation of the input operators,  
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%W
% AT, 02/18/2026: Initial coding;

if nargin==1
    Fdiag=varargin{1};
else % sequentially creating a block diagonal matrix
% Extract the operators
    a = varargin{1};    b = varargin{2};
    
    % Currently support some matrix-opvar concatenation
    if ~isa(b,'ndopvar') && ~isa(a, 'ndopvar')
        error("Currently supported only for ndopvar or nopvar");
        % Only supported if a in fact corresponds to a matrix as well
    end
    
    % Check that domain and variables match
    if any(any(a.dom~=b.dom))|| any(~strcmp(a.vars(:, 1).varname,b.vars(:, 1).varname)) || any(~strcmp(a.vars(:, 2).varname,b.vars(:, 2).varname))
        error('Operators being concatenated have different intervals or different independent variables');
    end
    % Check that the output dimensions match
    a.dim = a.dim;  b.dim = b.dim;  % correction to make components have consistent dimensions 8/27-ss
    if any(a.dim(2)~=b.dim(2))
        error('Cannot concatenate vertically: Input dimensions of opvar objects do not match')
    end  
    
    % convert to the same degree
    if any(a.deg~=b.deg)
        max_degree = max(a.deg, b.deg);
        Aop_new = change_degree(a, max_degree);
        Bop_new = change_degree(b, max_degree);
        a = Aop_new;
        b = Bop_new;
    end
    
    % if operators have different decvar, convert to the same decvar
    if isa(a, 'ndopvar')
        Aop_dvarname = a.dvarname;
    else
        Aop_dvarname = {};
    end
    if isa(b, 'ndopvar')
        Bop_dvarname = b.dvarname;
    else
        Bop_dvarname = {};
    end
    if numel(Aop_dvarname) ~= numel(Bop_dvarname) || ~isequal(Aop_dvarname,Bop_dvarname)
        % choose common dvar
        dvars1 = string(Aop_dvarname); % convert array to char array
        dvars2 = string(Bop_dvarname);
        % numberOfCharacters = max(size(dvars1{1}, 2), size(dvars2{1}, 2));
        % dvars1 = pad(dvars1, numberOfCharacters); % pad with ' ' if needed
        % dvars2 = pad(dvars2, numberOfCharacters); % pad with ' ' if needed 
    
        if isempty(dvars1) || isempty(dvars2)
            common_dvar = [];
            full_dvars = [dvars1; dvars2];
        else
            common_dvar = intersect(dvars1, dvars2, 'rows');
            new_dvars   = setdiff(dvars2, common_dvar, 'rows');
            full_dvars = char([dvars1; new_dvars]);
        end
        a = change_dec_var(a, full_dvars); % convert to the same dec var
        b = change_dec_var(b, full_dvars); % convert to the same dec var
    end
    
    % Finally, let's actually concatenate
    Fdiag = ndopvar(); % empty nopvar
    Fdiag.dom = a.dom;
    Fdiag.deg = a.deg;
    Fdiag.vars = a.vars;
    Fdiag.dvarname = a.dvarname;
    N = size(a.dom,1);
    Fdiag.C = cell([3*ones(1,N),1]);
    for ii=1:numel(a.C)
        Fdiag.C{ii} = blkdiag(a.C{ii}, b.C{ii});
    end
    
    if nargin>2 % repeat when there are more than two block elements
        Fdiag = blkdiag(Fdiag,varargin{3:end});
    end
end
end