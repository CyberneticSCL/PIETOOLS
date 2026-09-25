classdef (InferiorClasses={?polynomial,?sopvar,?sdopvar,?copvar}) cdopvar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CDOPVAR  Container for PI operators between concatenated mixed L2 spaces
% that may DEPEND ON DECISION VARIABLES. Sec. 8 of the sopvar document.
% An M x N grid of 'sdopvar' or 'sopvar' blocks,
%
%   P.C{i,j} : L_2^{p_j}[s^j] -> L_2^{q_i}[s^i],   (P*x)_i = sum_j P.C{i,j}*x_j
%
% For a purely fixed container use 'copvar'. The two stand to each other as
% 'sopvar' does to 'sdopvar': 'copvar' is closed under composition,
% 'cdopvar' is not, since the product of two operators affine in the
% decision variables is quadratic in them.
%
% BLOCKS NEED NOT ALL BE 'sdopvar'. A mixed container is not an oddity, it
% is reached by ordinary arithmetic: adding a fixed operator to a decision
% operator that has a zero block leaves a 'sopvar' block among 'sdopvar'
% ones, because a zero block passes through addition untouched. Such a
% container is legal, verifies, and composes. That is also why a container's
% decision-ness cannot be read off its block classes one at a time - it is a
% property of the container, which is what this class name states.
%
% ROW IS OUTPUT, COLUMN IS INPUT. Sec. 8 also writes the blocks inline as
% "P_{i,j}: L_2^{p_i}[s^i] -> L_2^{q_j}[s^j]", putting i on the input side;
% that contradicts its own display of P acting on [x_1;...;x_M] and its
% partition s_1^{ij} = s^j/s^i (input only), s_2^{ij} = s^i/s^j (output
% only). The partition is the correct reading and is implemented here.
%
% CLASS properties
% - P.C:         M x N cell of 'sdopvar' or 'sopvar'. A cell holding [] is a
%                STRUCTURALLY ZERO block: it stores nothing and its space
%                and dimension come from the metadata below. Sec. 8 has no
%                convention for an absent block, but a PIE is mostly zero;
% - P.vars:      1 x nv sorted cellstr registry of every spatial variable;
% - P.dom:       nv x 2, dom(k,:) = [a,b] the domain of P.vars{k};
% - P.space_out: M x nv logical, row i is the mask of s^i in P.vars;
% - P.space_in:  N x nv logical, row j is the mask of s^j;
% - P.dim_out:   M x 1 component counts q_i;
% - P.dim_in:    N x 1 component counts p_j;
% - P.Zd:        the SINGLE decision variable list shared by every
%                'sdopvar' block; cell(0,1) when there are none.
%
% Metadata placement, row and column consistency, and the Sec. 8.3.1
% decision not to store any row- or column-synchronized representation are
% all as documented in 'copvar'; this class differs from it only by
% admitting 'sdopvar' blocks and by carrying P.Zd. The validation itself is
% shared, in 'derive_copvar_meta'.
%
% DECISION VARIABLES ARE A CLASS INVARIANT. Every 'sdopvar' block carries
% the same P.Zd; the constructor merges them when they differ. Sec. 8.1
% calls this secondary; it is not. Zd is a cellstr of NAMES and measures
% 88-99% of a block's storage over q = 1e4..1e6. Copy-on-write shares one
% list across blocks for free (16 blocks at q = 2e5: 0 MB extra), but blocks
% built independently share nothing (322 MB for the same 16), and one
% differing name makes each pairwise operation pay an O(q log q) setdiff
% over names: 'plus' at q = 1e6 costs 0.063 s on a common list against
% 1.566 s on lists differing in one name. Holding the list here makes the
% sharing structural rather than accidental. It is also the truth for the
% operators that matter: a positive operator is Z*QZ for one Gram Q, so
% every block is affine in the same vec(Q).
%
% COST. Construction is O(M*N) block visits plus O(q) string comparisons per
% 'sdopvar' block for the Zd check. Nothing densifies or concatenates along
% the decision variable axis. Class methods pass validated metadata to the
% two-argument constructor, so the O(M*N*q) check is paid once at
% construction, not once per operation.
%
% See also COPVAR, VERIFY, SIZE, PLUS, CTRANSPOSE, MTIMES, SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - cdopvar
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
% Initial coding MMP, 09/17/2026. Split out of 'copvar', which until this
%                date held both cases in one class, so that the containers
%                mirror the sopvar/sdopvar pair and 'class(P)' states
%                whether decision variables are present. 'copvar' declared
%                ?sdopvar inferior from 09/11/2026, which was the basis for
%                the single-class reading; that declaration is now inert
%                there and the decision case lives here.
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  derive_mopvar_meta -> derive_copvar_meta. File was
%                  'mdopvar.m'.

    properties
        C = {};                     % M x N cell of blocks; [] = zero block
        vars = cell(1,0);           % 1 x nv global variable registry
        dom = zeros(0,2);           % nv x 2 domains, row k for vars{k}
        space_out = false(0,0);     % M x nv mask of s^i in vars
        space_in = false(0,0);      % N x nv mask of s^j in vars
        dim_out = zeros(0,1);       % M x 1 component counts q_i
        dim_in = zeros(0,1);        % N x 1 component counts p_j
        Zd = cell(0,1);             % the one decision variable list
    end

    methods
        function P = cdopvar(varargin)
            % P = CDOPVAR(C) builds a container from an M x N cell of
            % 'sdopvar'/'sopvar' blocks, deriving and validating all
            % metadata and putting every block on one decision variable
            % list. [] is a zero block.
            %
            % P = CDOPVAR(Q) promotes a 'copvar' Q, keeping its blocks and
            % metadata and leaving Zd empty.
            %
            % P = CDOPVAR(C,meta) trusts 'meta' and skips validation, used
            % by the class methods, which already know their result's
            % metadata and must not repeat the O(M*N*q) Zd check.
            %
            % CDOPVAR P1 P2 declares empty containers in the caller.
            if nargin==0
                return
            end
            if nargout==0 && (ischar(varargin{1}) || isstring(varargin{1}))
                for i = 1:nargin
                    assignin('caller',char(varargin{i}),cdopvar());
                end
                return
            end
            if isa(varargin{1},'copvar')
                Q = varargin{1};
                meta = metadata(Q);     meta.Zd = cell(0,1);
                P.C = Q.C;
                P.vars      = meta.vars;            P.dom       = meta.dom;
                P.space_out = meta.space_out;       P.space_in  = meta.space_in;
                P.dim_out   = meta.dim_out;         P.dim_in    = meta.dim_in;
                P.Zd        = meta.Zd;
                return
            end
            C = varargin{1};
            if ~iscell(C)
                error('cdopvar:badBlocks','Blocks must be an M x N cell array.')
            end
            if nargin>=2
                meta = varargin{2};
            else
                meta = derive_copvar_meta(C,'cdopvar',{'sopvar','sdopvar'});
                [C,meta.Zd] = unify_dvars(C);
            end
            P.C = C;
            P.vars      = meta.vars;            P.dom       = meta.dom;
            P.space_out = meta.space_out;       P.space_in  = meta.space_in;
            P.dim_out   = meta.dim_out;         P.dim_in    = meta.dim_in;
            P.Zd        = meta.Zd;
        end

        function meta = metadata(P)
            % META = METADATA(P) returns the validated metadata of P, for
            % handing to the two-argument constructor.
            meta = struct('vars',{P.vars},'dom',P.dom,...
                'space_out',P.space_out,'space_in',P.space_in,...
                'dim_out',P.dim_out,'dim_in',P.dim_in,'Zd',{P.Zd});
        end
    end
end


% ========================================================================
function [C,Zd] = unify_dvars(C)
% Put every 'sdopvar' block on one decision variable list and return it.
%
% The numel guard is O(1) and settles the usual mismatch without touching
% the names; 'isequal' is O(q) string comparisons and is only reached once
% the lengths already agree.

dec = find(cellfun(@(b) isa(b,'sdopvar'),C(:)))';
Zd = cell(0,1);
if isempty(dec)
    return
end
Zd = C{dec(1)}.Zd(:);
lists = cell(1,numel(dec));
same = true;
for k = 1:numel(dec)
    lists{k} = C{dec(k)}.Zd(:);
    if numel(lists{k})~=numel(Zd) || ~isequal(lists{k},Zd)
        same = false;
    end
end
if ~same
    % One merge for all blocks, not pairwise against a growing list:
    % pairwise re-synchronization profiled at 32% of 'possopvar' at three
    % spatial variables. 'unique(...,''stable'')' over the concatenation
    % gives the same list iterated setdiff would, both keeping first
    % occurrences in order.
    Zd = unique(vertcat(lists{:}),'stable');
end
for k = 1:numel(dec)
    if ~same
        C{dec(k)} = setdvars(C{dec(k)},Zd);     % remaps the rows of B
    end
    % Assign the one array itself, so all blocks share its storage.
    % 'setdvars' stores Zd(:).', a fresh q-length pointer array per block,
    % which costs 8*q bytes each (24.4 MB over 16 blocks at q = 2e5); this
    % makes the arrays identical instead, measured at 0 MB of extra memory.
    % It also normalizes the orientation to a column.
    C{dec(k)}.Zd = Zd;
end
end
