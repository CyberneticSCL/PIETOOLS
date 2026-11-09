function Pop = rand_cdopvar(spaces,dims,dom,deg,ndec,dnsty,occ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POP = RAND_CDOPVAR(SPACES,DIMS,DOM,DEG,NDEC,DNSTY,OCC) generates a random
% 'cdopvar' container over the given output and input spaces.
%
% INPUTS
% - spaces, dims, dom, deg, dnsty, occ: as for 'rand_copvar';
% - ndec:   number of decision variables, or an M x N array giving a
%           per-block count. A scalar gives every block the same count and
%           hence the same name list, which is the normal case;
%
% OUTPUTS
% - Pop:    'cdopvar' object whose blocks are 'sdopvar';
%
% NOTES
% Blocks come from 'rand_sdopvar', which takes a COUNT and generates its own
% names, so each block arrives with its own list object. The container
% constructor then collapses them: with a scalar ndec the lists agree in
% content, so it detects that and hands every block the one array, which is
% the cheap path and the state a real operator is in. Passing an M x N ndec
% forces the genuine merge instead, which is how to exercise that path
% deliberately - it is the expensive one, an O(q log q) setdiff over names.
%
% For a container with NO decision variables use 'rand_copvar'; a 'cdopvar'
% with ndec = 0 in every block is legal but is just a fixed operator in the
% decision container, which is only worth building to test that case.
%
% See also RAND_COPVAR, RAND_SDOPVAR, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - rand_cdopvar
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
% Initial coding MMP, 09/18/2026
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  rand_mdopvar -> rand_cdopvar, rand_mopvar -> rand_copvar.
%                  File was 'rand_mdopvar.m'.

if ~isstruct(spaces) || ~isfield(spaces,'out') || ~isfield(spaces,'in')
    error('rand_cdopvar:badSpaces',...
        "Spaces must be a struct with fields 'out' and 'in'.")
end
M = numel(spaces.out);      N = numel(spaces.in);
if nargin<6 || isempty(dnsty),  dnsty = 0.4;        end
if nargin<7 || isempty(occ),    occ = true(M,N);    end
if isscalar(deg),   deg = struct('out',deg,'in',deg);   end
if isscalar(ndec),  ndec = ndec*ones(M,N);          end
if ~isequal(size(ndec),[M,N])
    error('rand_cdopvar:badNdec','ndec must be scalar or %dx%d.',M,N)
end
dim_out = expand(dims.out,M,'dims.out');
dim_in  = expand(dims.in ,N,'dims.in');

if M>0 && N>0
    if any(~any(occ,2))
        error('rand_cdopvar:emptyRow',...
            'occ leaves row(s) %s empty; each row needs one populated block.',...
            mat2str(find(~any(occ,2))'))
    end
    if any(~any(occ,1))
        error('rand_cdopvar:emptyColumn',...
            'occ leaves column(s) %s empty; each column needs one populated block.',...
            mat2str(find(~any(occ,1))))
    end
end

C = cell(M,N);
for i = 1:M
    for j = 1:N
        if ~occ(i,j),   continue,   end
        vs = struct('out',{spaces.out{i}},'in',{spaces.in{j}});
        ds = struct('out',dom_rows(dom,vs.out,spaces),...
                    'in' ,dom_rows(dom,vs.in ,spaces));
        gs = struct('out',deg.out*ones(1,numel(vs.out)),...
                    'in' ,deg.in *ones(1,numel(vs.in )));
        C{i,j} = rand_sdopvar([dim_out(i),dim_in(j)],vs,ds,gs,ndec(i,j),dnsty);
    end
end
Pop = cdopvar(C);

end

% ------------------------------------------------------------------------
function v = expand(v,n,nm)
if isscalar(v),     v = v*ones(n,1);    end
if numel(v)~=n
    error('rand_cdopvar:badDims','%s must have %d entries.',nm,n)
end
v = v(:);
end

function d = dom_rows(dom,v,spaces)
% One domain row per variable of this space, in this block's variable order;
% see the note in 'rand_copvar' on why a multi-row dom is looked up by name.
if size(dom,1)==1
    d = repmat(dom,numel(v),1);
    return
end
if ~isfield(spaces,'vars')
    error('rand_cdopvar:noVarList',...
        ['A dom with more than one row needs spaces.vars, a cellstr giving '...
         'the variable each row belongs to.'])
end
d = zeros(numel(v),2);
for k = 1:numel(v)
    idx = find(strcmp(spaces.vars,v{k}),1);
    if isempty(idx)
        error('rand_cdopvar:unknownVar',...
            'Variable ''%s'' is not listed in spaces.vars.',v{k})
    end
    d(k,:) = dom(idx,:);
end
end
