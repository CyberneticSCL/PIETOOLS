function info = verify_copvar_meta(P,cls,okclasses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO = VERIFY_COPVAR_META(P,CLS,OKCLASSES) checks that the blocks of a
% container agree with its metadata and with each other, for the 'verify'
% methods of 'copvar' and 'cdopvar'.
%
% INPUTS
% - P:          'copvar' or 'cdopvar' object;
% - cls:        its class, 'copvar' or 'cdopvar', named in the messages;
% - okclasses:  cellstr of block classes the container admits, {'sopvar'}
%               or {'sopvar','sdopvar'}, as for 'derive_copvar_meta';
%
% OUTPUTS
% - info.true:  1 if P is consistent, 0 otherwise;
% - info.flags: cell array of messages, one per inconsistency, empty when
%               P is consistent;
%
% NOTES
% CHECKED: metadata shapes; no all-empty row or column; per block, that it
% is of an admitted class with 1x2 dims, that its output side matches its
% row and its input side matches its column, that it puts every variable on
% the registry's domain, and for a 'cdopvar' that an 'sdopvar' block carries
% the container's Zd.
%
% This is the report-mode counterpart of 'derive_copvar_meta', which throws
% on the same defects while deriving the metadata of a new container: here
% the metadata already exists, and each block is checked against it (derive
% checks each block against the first block of its row or column).
%
% The space checks are SET equality of the names, not 'isequal' on the
% lists: the block constructors enforce vars.out = [S2,S3] and vars.in =
% [S3,S1], so blocks in one row order the same output set differently
% whenever their input spaces differ, and blocks in one column likewise.
%
% Cost: O(M*N) block visits, spatial set operations on the registry, and
% O(q) string comparisons per 'sdopvar' block guarded by the size test in
% 'isequal'. P is a container read in a plain function, so each P.name read
% goes through the overloaded subsref (about 13 us, 'eq_copvar'): every
% property is read once, into a struct. No decision variable list is copied
% unless it is stored as a row.
%
% See also VERIFY, DERIVE_COPVAR_META, COPVAR, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - verify_copvar_meta
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
% Initial coding MMP, 09/30/2026. The body of @copvar/verify (MMP,
%                09/17/2026) and @cdopvar/verify (MMP, 09/17/2026), which
%                were one algorithm in two copies, differing only in the
%                admitted block classes and the Zd check; both methods now
%                call this. Changes from those copies, none in a verdict or
%                a message: the admitted classes are an argument, each
%                property of P is read once into the struct 'm', and the
%                Zd lists are compared as columns without copying a column.

if ~isa(P,cls)
    error('verify:badInput','Input must be a %s object.',cls)
end
hasZd = strcmp(cls,'cdopvar');

% One read per property: P.(name) in this plain function is the overloaded
% subsref. Not 'metadata(P)', whose struct() call would build a struct
% ARRAY from a malformed cell-valued property, which is what 'verify' is
% for; a direct read keeps any value as it is.
fn = {'vars','dom','space_out','space_in','dim_out','dim_in'};
for k = 1:numel(fn),    m.(fn{k}) = P.(fn{k});    end
C = P.C;
if hasZd
    Zd = P.Zd;
    if ~iscolumn(Zd),   Zd = Zd(:);     end
end

flags = cell(0,1);
M = size(C,1);      N = size(C,2);      % @copvar/size, @cdopvar/size
vars = m.vars;      nv = numel(vars);
occ = ~cellfun(@isempty,C);

% % % Metadata shapes. The per-block checks index into this metadata, so
% % % report shape problems and stop rather than cascade consequences.
want = {'dom',[nv,2]; 'space_out',[M,nv]; 'space_in',[N,nv];
        'dim_out',[M,1]; 'dim_in',[N,1]};
for k = 1:size(want,1)
    got = size(m.(want{k,1}));
    if ~isequal(got,want{k,2})
        flags{end+1,1} = sprintf('%s is %s; expected %s.',...
            want{k,1},mat2str(got),mat2str(want{k,2})); %#ok<AGROW>
    end
end
if ~isempty(flags)
    info.true = 0;      info.flags = flags;     return
end

% % % No structurally empty row or column: such a row has no block to say
% % % what its space and dimension are.
for i = find(~any(occ,2))'
    flags{end+1,1} = sprintf('Row %d contains no populated block.',i); %#ok<AGROW>
end
for j = find(~any(occ,1))
    flags{end+1,1} = sprintf('Column %d contains no populated block.',j); %#ok<AGROW>
end

% % % Per-block agreement with the container.
for i = 1:M
    s_i = sort(vars(m.space_out(i,:)));
    for j = 1:N
        if ~occ(i,j),   continue,   end
        b = C{i,j};
        s_j = sort(vars(m.space_in(j,:)));
        if ~any(cellfun(@(c) isa(b,c),okclasses)) || numel(b.dims)~=2
            flags{end+1,1} = sprintf(...
                'Block (%d,%d) is a ''%s'' with %d matrix dimensions; want %s with 2.',...
                i,j,class(b),numel(b.dims),strjoin(okclasses,' or ')); %#ok<AGROW>
            continue
        end
        if b.dims(1)~=m.dim_out(i)
            flags{end+1,1} = sprintf('Block (%d,%d) has output dimension %d; row %d is %d.',...
                i,j,b.dims(1),i,m.dim_out(i)); %#ok<AGROW>
        end
        if b.dims(2)~=m.dim_in(j)
            flags{end+1,1} = sprintf('Block (%d,%d) has input dimension %d; column %d is %d.',...
                i,j,b.dims(2),j,m.dim_in(j)); %#ok<AGROW>
        end
        if ~isequal(sort(b.vars.out(:))',s_i(:)')
            flags{end+1,1} = sprintf('Block (%d,%d) maps into {%s}; row %d is {%s}.',...
                i,j,strjoin(sort(b.vars.out(:))',','),i,strjoin(s_i(:)',',')); %#ok<AGROW>
        end
        if ~isequal(sort(b.vars.in(:))',s_j(:)')
            flags{end+1,1} = sprintf('Block (%d,%d) maps out of {%s}; column %d is {%s}.',...
                i,j,strjoin(sort(b.vars.in(:))',','),j,strjoin(s_j(:)',',')); %#ok<AGROW>
        end
        % Domains against the registry, not against another block: a
        % variable given two domains is invisible from inside one block.
        flags = check_dom(flags,i,j,b.vars.out,b.dom.out,vars,m.dom,'output');
        flags = check_dom(flags,i,j,b.vars.in ,b.dom.in ,vars,m.dom,'input');
        % Decision variables; 'isequal' compares sizes first, so a length
        % mismatch is settled without comparing q names. A 'sopvar' block
        % has none to check.
        if hasZd && isa(b,'sdopvar')
            zb = b.Zd;
            if ~iscolumn(zb),   zb = zb(:);     end
            if ~isequal(zb,Zd)
                flags{end+1,1} = sprintf(['Block (%d,%d) carries a decision variable '...
                    'list of length %d, differing from the container''s %d.'],...
                    i,j,numel(b.Zd),numel(Zd)); %#ok<AGROW>
            end
        end
    end
end

info.true = double(isempty(flags));
info.flags = flags;

end


% ========================================================================
function flags = check_dom(flags,i,j,v,d,vars,dom,side)
% Compare a block's domains, listed in its own variable order, against the
% container registry.
v = v(:)';
if size(d,1)~=numel(v)
    flags{end+1,1} = sprintf('Block (%d,%d) lists %d %s variables against %d domain rows.',...
        i,j,numel(v),side,size(d,1));
    return
end
for k = 1:numel(v)
    idx = find(strcmp(vars,v{k}),1);
    if isempty(idx)
        flags{end+1,1} = sprintf('Block (%d,%d) uses %s variable ''%s'', not in the registry.',...
            i,j,side,v{k}); %#ok<AGROW>
    elseif ~isequal(d(k,:),dom(idx,:))
        flags{end+1,1} = sprintf(['Block (%d,%d) puts %s variable ''%s'' on '...
            '[%g,%g]; the registry says [%g,%g].'],...
            i,j,side,v{k},d(k,1),d(k,2),dom(idx,1),dom(idx,2)); %#ok<AGROW>
    end
end
end
