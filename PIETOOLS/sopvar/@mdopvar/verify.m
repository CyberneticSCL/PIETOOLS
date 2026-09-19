function info = verify(Mop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO = VERIFY(MOP) checks that the blocks of a 'mdopvar' agree with its
% container metadata and with each other.
%
% OUTPUTS
% - info.true:   1 if MOP is consistent, 0 otherwise;
% - info.flags:  cell array of messages, one per inconsistency, empty when
%                MOP is consistent;
%
% CHECKED: metadata shapes; no all-empty row or column; per block, that its
% output side matches its row and its input side matches its column, that it
% puts every variable on the registry's domain, and that an 'sdopvar' block
% carries the container's Zd.
%
% The space checks are SET equality of the names, not 'isequal' on the
% lists: the block constructors enforce vars.out = [S2,S3] and vars.in =
% [S3,S1], so blocks in one row order the same output set differently
% whenever their input spaces differ, and blocks in one column likewise.
% Comparing with 'isequal' would reject legitimate containers.
%
% The constructor already rejects anything reported here, so a 'mdopvar'
% built through 'mdopvar(C)' always verifies. VERIFY is for objects reached
% another way: properties assigned directly, the two-argument trusted
% constructor, or a block replaced in place via P.C{i,j} = ....
%
% Cost: O(M*N) block visits, spatial set operations on the registry, and
% O(q) string comparisons per 'sdopvar' block guarded by an O(1) length test.
%
% See also MDOPVAR, MOPVAR, SIZE.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - verify(mdopvar)
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
% Initial coding MMP, 09/17/2026. Split out of @mopvar/verify, which until
%                this date checked both cases; that file keeps the history
%                of MP's 01/19/2026 draft. This copy differs from it only by
%                admitting 'sdopvar' blocks and by checking the decision
%                variable list.

if ~isa(Mop,'mdopvar')
    error('verify:badInput','Input must be a mdopvar object.')
end

flags = cell(0,1);
[M,N] = size(Mop);
C = Mop.C;      vars = Mop.vars;        nv = numel(vars);
occ = ~cellfun(@isempty,C);

% % % Metadata shapes. The per-block checks index into this metadata, so
% % % report shape problems and stop rather than cascade consequences.
want = {'dom',[nv,2]; 'space_out',[M,nv]; 'space_in',[N,nv];
        'dim_out',[M,1]; 'dim_in',[N,1]};
for k = 1:size(want,1)
    got = size(Mop.(want{k,1}));
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
    s_i = sort(vars(Mop.space_out(i,:)));
    for j = 1:N
        if ~occ(i,j),   continue,   end
        b = C{i,j};
        s_j = sort(vars(Mop.space_in(j,:)));
        if (~isa(b,'sopvar') && ~isa(b,'sdopvar')) || numel(b.dims)~=2
            flags{end+1,1} = sprintf(...
                'Block (%d,%d) is a ''%s'' with %d matrix dimensions; want sopvar or sdopvar with 2.',...
                i,j,class(b),numel(b.dims)); %#ok<AGROW>
            continue
        end
        if b.dims(1)~=Mop.dim_out(i)
            flags{end+1,1} = sprintf('Block (%d,%d) has output dimension %d; row %d is %d.',...
                i,j,b.dims(1),i,Mop.dim_out(i)); %#ok<AGROW>
        end
        if b.dims(2)~=Mop.dim_in(j)
            flags{end+1,1} = sprintf('Block (%d,%d) has input dimension %d; column %d is %d.',...
                i,j,b.dims(2),j,Mop.dim_in(j)); %#ok<AGROW>
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
        flags = check_dom(flags,i,j,b.vars.out,b.dom.out,vars,Mop.dom,'output');
        flags = check_dom(flags,i,j,b.vars.in ,b.dom.in ,vars,Mop.dom,'input');
        % Decision variables; length first, so a mismatch is usually settled
        % without comparing q names. A 'sopvar' block has none to check.
        if isa(b,'sdopvar') && (numel(b.Zd)~=numel(Mop.Zd) || ~isequal(b.Zd(:),Mop.Zd(:)))
            flags{end+1,1} = sprintf(['Block (%d,%d) carries a decision variable '...
                'list of length %d, differing from the container''s %d.'],...
                i,j,numel(b.Zd),numel(Mop.Zd)); %#ok<AGROW>
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

