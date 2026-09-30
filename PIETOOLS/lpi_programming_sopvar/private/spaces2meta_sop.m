function meta = spaces2meta_sop(dims,spaces,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% META = SPACES2META_SOP(DIMS,SPACES,DOM) turns the space description that
% 'lpivar_cdopvar' and 'poscopvar' take into container metadata (vars, dom,
% space_out, space_in, dim_out, dim_in), for the constructors of this folder.
%
% INPUTS (the conventions of 'lpivar_cdopvar')
% - dims:   struct with fields 'out' (M x 1) and 'in' (N x 1), or an M x 1
%           array for the same spaces in and out; a scalar is expanded;
% - spaces: struct with fields 'out' (1 x M cell) and 'in' (1 x N cell) of
%           cellstr variable names, one per space, {} being R^q; or a
%           1 x M cell for the same spaces in and out; a plain cellstr is
%           one space;
% - dom:    nv x 2 in the order of the SORTED registry, one 1 x 2 row for
%           every variable, or a struct with fields 'vars' and 'dom';
%
% OUTPUTS
% - meta:   the container metadata, registry sorted.
%
% NOTES
% The three parsers below are copies of the local functions of
% 'lpivar_cdopvar' (MMP, 09/25/2026), unchanged, so this folder accepts
% exactly what that routine accepts; they are local functions there and
% could not be called, and editing that routine was outside this change.
%
% Cost: O(M + N) spaces and set operations on the registry, one entry per
% spatial direction; nothing depends on decision variables.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - spaces2meta_sop
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
% Initial coding MMP, 09/29/2026

[sp_out,sp_in,d_out,d_in] = parse_spaces(dims,spaces);
M = numel(sp_out);      N = numel(sp_in);
vars = reshape(unique([sp_out{:}, sp_in{:}]),1,[]);     nv = numel(vars);
% '_int' and '_dum' suffixes are the classes' internal variable names.
is_reserved = ~cellfun(@isempty,regexp(vars,'_(int|dum)$','once'));
if any(is_reserved)
    error("Spatial variable names may not end in '_int' or '_dum'; "...
          +"'"+string(vars{find(is_reserved,1)})+"' does.")
end
dom = parse_dom(dom,vars);
so = false(M,nv);       si = false(N,nv);
for i = 1:M,    so(i,:) = ismember(vars,sp_out{i});    end
for j = 1:N,    si(j,:) = ismember(vars,sp_in{j});     end
meta = struct('vars',{vars},'dom',dom,'space_out',so,'space_in',si,...
    'dim_out',d_out(:),'dim_in',d_in(:));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copied unchanged from 'lpivar_cdopvar' (see NOTES above).
function [so,si,dout,din] = parse_spaces(dims,spaces)
% Output and input space lists, and component counts.
if isstruct(spaces)
    if ~isscalar(spaces)
        error("'spaces' should be a single struct, not an array; struct() with "...
              +"a cell value returns an array, so assign the fields instead.")
    end
    if ~isfield(spaces,'out')
        error("A 'struct' spaces argument needs a field 'out'.")
    end
    so = norm_spaces(spaces.out);
    if isfield(spaces,'in'),    si = norm_spaces(spaces.in);
    else,                       si = so;
    end
else
    so = norm_spaces(spaces);   si = so;
end
M = numel(so);      N = numel(si);
if isstruct(dims)
    dout = reshape(dims.out,[],1);
    if isfield(dims,'in'),  din = reshape(dims.in,[],1);
    else,                   din = dout;
    end
else
    dout = reshape(dims,[],1);  din = dout;
end
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
% A plain cellstr is ONE space; otherwise a cell of cellstr, one per space.
if isempty(sp),     sp = {cell(1,0)};   return,     end
if iscellstr(sp) || ischar(sp) || isa(sp,'polynomial'),  sp = {sp};  end
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
    sk = reshape(sk,1,[]);
    if numel(unique(sk))~=numel(sk)
        error("Space %d repeats a variable name.",k)
    end
    sp{k} = sk;
end
end


function dom = parse_dom(dom,vars)
% As 'copquadvar': nv x 2 in registry order, one 1 x 2 row for all, or a
% struct pairing names with intervals.
nv = numel(vars);
if isa(dom,'struct')
    if ~isfield(dom,'vars') || ~isfield(dom,'dom')
        error("A 'struct' domain should have fields 'vars' and 'dom'.")
    end
    if ~isscalar(dom)
        error("A 'struct' domain should be a single struct, not an array; "...
              +"struct() with a cell value returns an array.")
    end
    dvars = dom.vars;
    if isa(dvars,'polynomial'),  dvars = dvars.varname(:)';   end
    if ischar(dvars),            dvars = {dvars};             end
    dvars = reshape(dvars,1,[]);
    if size(dom.dom,1)~=numel(dvars) || size(dom.dom,2)~=2
        error("A 'struct' domain should pair each name in 'vars' with a row of 'dom'.")
    end
    [tf,loc] = ismember(vars,dvars);
    if ~all(tf)
        error("No domain was given for the variable '"+string(vars{find(~tf,1)})+"'.")
    end
    dom = dom.dom(loc,:);
end
if nv==0
    dom = zeros(0,2);
    return
end
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
