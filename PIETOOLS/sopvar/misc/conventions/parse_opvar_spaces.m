function [meta,sp_out,sp_in] = parse_opvar_spaces(dims,spaces,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [META,SP_OUT,SP_IN] = PARSE_OPVAR_SPACES(DIMS,SPACES,DOM) reads the
% description of the output and input spaces of a container into the
% container metadata.
%
% INPUT
% - dims:   struct with fields 'out' (M x 1) and 'in' (N x 1), the component
%           counts; or an M x 1 array for input spaces equal to the output
%           spaces. A scalar is expanded to every space; an empty dims.out
%           (dims.in) gives no output (input) spaces, an operator with no
%           rows (columns);
% - spaces: scalar struct with fields 'out' (1 x M cell) and 'in' (1 x N
%           cell) of 'cellstr' variable names, one per space, an empty entry
%           being the finite-dimensional space R^q; 'in' may be omitted. Or
%           a 1 x M cell for input spaces equal to the output spaces. A plain
%           'cellstr', a char or a 'polynomial' is ONE space; an empty
%           input is one space R^q. A struct array is refused;
% - dom:    nv x 2 domains in the order of the SORTED registry, one 1 x 2
%           row for every variable, or a scalar struct with fields 'vars'
%           (cellstr, char or 'polynomial') and 'dom' pairing names with
%           rows;
%
% OUTPUT
% - meta:   struct with fields, in this order,
%             vars       1 x nv cellstr, the sorted union of all spaces;
%             dom        nv x 2, dom(d,:) the interval of vars{d};
%             space_out  M x nv logical, row i marking the variables of
%                        output space i;
%             space_in   N x nv logical, likewise for the input spaces;
%             dim_out    M x 1 component counts;
%             dim_in     N x 1 component counts;
%           the metadata of an 'opvar' or 'dopvar' container, less the
%           decision variable list;
% - sp_out: 1 x M cell of 1 x n_i cellstr, output space i by name, in the
%           order given;
% - sp_in:  1 x N cell, likewise for the input spaces.
%
% NOTES
% The variables of a space are those of L_2[s]: a block from input space j
% to output space i shares S3 = s^i n s^j and has S2, S1 as created and
% lost variables. The registry is sorted because a block indexes its
% parameter cell over S3 in sorted order; with a sorted registry every
% space, and the S3 of every pair, is already in that order.
%
% Errors are raised in the order spaces, dimensions, reserved names,
% domain. Error messages carry no identifier.
%
% Cost: O(M + N) spaces and set operations on the registry, one name per
% spatial direction; nothing depends on decision variables.
%
% See also CHECK_RESERVED_NAMES, LPIVAR, OPQUADVAR, DOPVAR,
% DERIVE_OPVAR_META.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - parse_opvar_spaces
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
% Initial coding MMP, 09/30/2026

% Space lists and component counts of both sides
[sp_out,sp_in,d_out,d_in] = parse_spaces(dims,spaces);
M = numel(sp_out);      N = numel(sp_in);
% The registry: the sorted union of the variables of all spaces, with no
% reserved suffix, and its domains
vars = [sp_out{:}, sp_in{:}];       if isempty(vars),   vars = cell(1,0);   end
vars = reshape(unique(vars),1,[]);     nv = numel(vars);
check_reserved_names(vars);
dom = parse_dom(dom,vars);
% Masks of the variables of each space over the registry
so = false(M,nv);       si = false(N,nv);
for i = 1:M,    so(i,:) = ismember(vars,sp_out{i});    end
for j = 1:N,    si(j,:) = ismember(vars,sp_in{j});     end
meta = struct('vars',{vars},'dom',dom,'space_out',so,'space_in',si,...
    'dim_out',d_out(:),'dim_in',d_in(:));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [so,si,dout,din] = parse_spaces(dims,spaces)
% [SO,SI,DOUT,DIN] = PARSE_SPACES(DIMS,SPACES) reads the output and input
% space lists and component counts from the accepted forms of DIMS and
% SPACES.
% INPUTS  dims, spaces: as for the primary function.
% OUTPUT  so, si: 1 x M and 1 x N cells of cellstr; dout, din: M x 1 and
%         N x 1 counts.
% Space lists: a scalar struct with 'out' (and optionally 'in'), or one
% list for both sides
if isstruct(spaces)
    if ~isscalar(spaces)
        error("'spaces' should be a single struct, not an array; struct() with "...
              +"a cell value returns an array, so assign the fields instead.")
    end
    if ~isfield(spaces,'out')
        error("A 'struct' spaces argument needs a field 'out'.")
    end
    % Input spaces default to the output spaces
    so = norm_spaces(spaces.out);
    if isfield(spaces,'in'),    si = norm_spaces(spaces.in);
    else,                       si = so;
    end
else
    so = norm_spaces(spaces);   si = so;
end
% Component counts: a struct with 'out' (and optionally 'in'), or one
% array for both sides
if isstruct(dims)
    dout = reshape(dims.out,[],1);
    if isfield(dims,'in'),  din = reshape(dims.in,[],1);
    else,                   din = dout;
    end
else
    dout = reshape(dims,[],1);  din = dout;
end
% Empty dimensions specify no output (input) spaces at all, i.e. no rows
% (columns) of the container, whatever the spaces specified
if isempty(dout),   so = cell(1,0);     end
if isempty(din),    si = cell(1,0);     end
% A scalar count applies to every space; one positive integer per space
M = numel(so);      N = numel(si);
if isscalar(dout),  dout = repmat(dout,M,1);    end
if isscalar(din),   din = repmat(din,N,1);      end
if numel(dout)~=M || numel(din)~=N
    error("Dimensions should have one entry per space: %d output and %d input.",M,N)
end
if any([dout;din]<1) || any([dout;din]~=round([dout;din]))
    error("Component counts should be positive integers.")
end
end


function sp = norm_spaces(sp)
% SP = NORM_SPACES(SP) normalizes a space list: a plain cellstr, char or
% pvar is ONE space; otherwise a cell with one cellstr per space.
% INPUT   sp: the list in any accepted form.
% OUTPUT  sp: 1 x M cell of 1 x n_i cellstr.
if isempty(sp),     sp = {cell(1,0)};   return,     end
if iscellstr(sp) || ischar(sp) || isa(sp,'polynomial'),  sp = {sp};  end
% Anything else must be a cell of spaces. Without this, brace indexing a
% non-cell below fails with the brace-indexing error of MATLAB.
if ~iscell(sp)
    error("Spaces should be specified as a cell of 'cellstr' objects.")
end
% Each space as a 1 x n cellstr: a 'polynomial' gives its variable names,
% a char one name, an empty numeric R^q; no repeated names
sp = reshape(sp,1,[]);
for k = 1:numel(sp)
    sk = sp{k};
    if isa(sk,'polynomial'),            sk = sk.varname(:)';    end
    if ischar(sk),                      sk = {sk};              end
    if isnumeric(sk) && isempty(sk),    sk = cell(1,0);         end
    if ~iscellstr(sk)
        error("Space %d should be a 'cellstr' of variable names; an empty "...
              +"one is the finite-dimensional space R^q.",k)
    end
    % A row of distinct names
    sk = reshape(sk,1,[]);
    if numel(unique(sk))~=numel(sk)
        error("Space %d repeats a variable name.",k)
    end
    sp{k} = sk;
end
end


function dom = parse_dom(dom,vars)
% DOM = PARSE_DOM(DOM,VARS) reads the domains: nv x 2 in registry order,
% one 1 x 2 row for all, or a struct pairing names with intervals.
% INPUTS  dom: in any accepted form; vars: the sorted registry.
% OUTPUT  dom: nv x 2 double.
nv = numel(vars);
% A scalar struct pairs names with rows: take the rows in registry order,
% one for every registry variable
if isa(dom,'struct')
    if ~isfield(dom,'vars') || ~isfield(dom,'dom')
        error("A 'struct' domain should have fields 'vars' and 'dom'.")
    end
    if ~isscalar(dom)
        error("A 'struct' domain should be a single struct, not an array; "...
              +"struct() with a cell value returns an array.")
    end
    % Names as a 1 x n cellstr, one row [a,b] of dom.dom per name
    dvars = dom.vars;
    if isa(dvars,'polynomial'),  dvars = dvars.varname(:)';   end
    if ischar(dvars),            dvars = {dvars};             end
    dvars = reshape(dvars,1,[]);
    if size(dom.dom,1)~=numel(dvars) || size(dom.dom,2)~=2
        error("A 'struct' domain should pair each name in 'vars' with a row of 'dom'.")
    end
    % The row of each registry variable, in registry order
    [tf,loc] = ismember(vars,dvars);
    if ~all(tf)
        error("No domain was given for the variable '"+string(vars{find(~tf,1)})+"'.")
    end
    dom = dom.dom(loc,:);
end
% No variables: an empty domain
if nv==0
    dom = zeros(0,2);
    return
end
% An nv x 2 array, or one row expanded to every variable, with
% dom(d,1) < dom(d,2)
if size(dom,2)~=2
    error("Domains should be specified as an nv x 2 array.")
end
if size(dom,1)==1 && nv~=1
    dom = repmat(dom,nv,1);
elseif size(dom,1)~=nv
    error("Domains should be an nv x 2 array for the %d registry variables.",nv)
end
if any(dom(:,2)<=dom(:,1))
    error("Each domain should satisfy dom(d,1) < dom(d,2).")
end
end
