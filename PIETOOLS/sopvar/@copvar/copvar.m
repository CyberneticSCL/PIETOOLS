classdef (InferiorClasses={?polynomial,?sopvar,?sdopvar,?dpvar}) copvar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPVAR  Container for FIXED PI operators between concatenated mixed L2
% spaces. An M x N grid of 'sopvar' blocks,
%
%   P.C{i,j} : L_2^{p_j}[s^j] -> L_2^{q_i}[s^i],   (P*x)_i = sum_j P.C{i,j}*x_j
%
% For a container that carries decision variables, use 'cdopvar', whose
% blocks may be 'sdopvar'. The two stand to each other as 'sopvar' does to
% 'sdopvar': 'copvar' is closed under composition, 'cdopvar' is not. A
% 'dpvar' operand (gam*P, gam - T'*T, [-gam, D']) dispatches to the methods
% of this class, and the result is a 'cdopvar'; so is the concatenation of
% a 'copvar' with an 'sdopvar' block.
%
% ROW IS OUTPUT, COLUMN IS INPUT: block (i,j) maps L_2^{p_j}[s^j] into
% L_2^{q_i}[s^i]; its input-only variables are s_1^{ij} = s^j/s^i and its
% output-only variables s_2^{ij} = s^i/s^j.
%
% CLASS properties
% - P.C:         M x N cell of 'sopvar'. A cell holding [] is a
%                STRUCTURALLY ZERO block: it stores nothing and its space
%                and dimension come from the metadata below;
% - P.vars:      1 x nv sorted cellstr registry of every spatial variable;
% - P.dom:       nv x 2, dom(k,:) = [a,b] the domain of P.vars{k};
% - P.space_out: M x nv logical, row i is the mask of s^i in P.vars;
% - P.space_in:  N x nv logical, row j is the mask of s^j;
% - P.dim_out:   M x 1 component counts q_i;
% - P.dim_in:    N x 1 component counts p_j;
% - (protected) leftCommonBasis, leftCommonC, rightCommonBasis,
%                rightCommonC: row- and column-common forms of the blocks;
%                no method of the class fills them, so they stay empty.
%
% Blocks still carry their own vars/dom/dims, so they stay valid standalone
% operators and no block method changes. The container holds the
% authoritative copy for two reasons: it is the only place a cross-block
% conflict is visible, since a block check sees one block at a time and so
% catches only a variable given two domains by the two SIDES of one block,
% not by two different blocks; and it determines the zero blocks, which
% have no metadata of their own.
%
% The blocks are left intact, so the container-level checks compare masks
% over P.vars but the per-block spatial set operations remain: each block
% computes vars_S1/S2/S3 by setdiff/intersect, and a block product
% re-partitions by name per pair. The container-level domain check is
% likewise ADDITIONAL to the per-block one, not a replacement.
%
% ROW AND COLUMN CONSISTENCY. Every populated block in row i maps into
% L_2^{q_i}[s^i], every block in column j out of L_2^{p_j}[s^j]. The space
% check is SET equality, not 'isequal' on the name lists: the block
% constructors enforce vars.out = [S2,S3] with S2^{ij} = s^i/s^j, so two
% blocks in one row order the same variable set differently whenever their
% input spaces differ.
%
% NOT STORED: any row- or column-synchronized auxiliary representation.
% ZL and ZR are metadata - together under 0.1 MB of a 66.85 MB
% 'sdopvar' block at q = 1e6 - so sharing them across a row saves nothing,
% while forcing blocks onto the row-union basis multiplies their coefficient
% column count, the union being taken per variable and then tensored.
% Measured over the subset lattice at 2 components per space, total
% coefficient columns grow x2.3 (nv=2, deg 1), x5.9 (nv=3, deg 2) and x15.3
% (nv=4, deg 3) under synchronization.
%
% COST. Construction is O(M*N) block visits. The spatial set operations run
% on the registry, one entry per direction. Class methods pass validated
% metadata to the two-argument constructor rather than re-deriving it.
%
% INPUTS
% - C:      M x N cell of 'sopvar' blocks, [] for a zero block;
% - meta:   (optional) struct with fields vars, dom, space_out, space_in,
%           dim_out, dim_in, trusted as is, validation skipped;
%
% OUTPUTS
% - P:      'copvar' container. Called with no output, as 'copvar P1 P2',
%           it declares empty containers P1, P2 in the caller;
%
% See also CDOPVAR, VERIFY, SIZE, PLUS, CTRANSPOSE, MTIMES, SOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - copvar
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
% MP, 01/15/2026: Initial coding

    properties
        C = {};                     % M x N cell of blocks; [] = zero block
        vars = cell(1,0);           % 1 x nv global variable registry
        dom = zeros(0,2);           % nv x 2 domains, row k for vars{k}
        space_out = false(0,0);     % M x nv mask of s^i in vars
        space_in = false(0,0);      % N x nv mask of s^j in vars
        dim_out = zeros(0,1);       % M x 1 component counts q_i
        dim_in = zeros(0,1);        % N x 1 component counts p_j
    end
    properties(Access=protected)
        % Common row-basis representation.
        leftCommonBasis = cell(0,1);    % M x 1, one common ZL per row
        leftCommonC     = cell(0,0);    % M x N sopvar blocks in row-common form
    
        % Common column-basis representation.
        rightCommonBasis = cell(0,1);   % N x 1, one common ZR per column
        rightCommonC     = cell(0,0);   % M x N sopvar blocks in column-common form
    end

    methods
        function P = copvar(varargin)
            % P = COPVAR(C) builds a container from an M x N cell of
            % 'sopvar' blocks, deriving and validating all metadata. [] is a
            % zero block.
            %
            % P = COPVAR(C,meta) trusts 'meta' and skips validation, used by
            % the class methods, which already know their result's metadata.
            %
            % COPVAR P1 P2 declares empty containers in the caller.

            % No input: an empty container
            if nargin==0
                return
            end
            % Command form: declare each named container in the caller
            if nargout==0 && (ischar(varargin{1}) || isstring(varargin{1}))
                for i = 1:nargin
                    assignin('caller',char(varargin{i}),copvar());
                end
                return
            end
            C = varargin{1};
            if ~iscell(C)
                error('copvar:badBlocks','Blocks must be an M x N cell array.')
            end
            % Trust the metadata if given; otherwise derive and validate it
            if nargin>=2
                meta = varargin{2};
            else
                meta = derive_copvar_meta(C,'copvar',{'sopvar'});
            end
            % Store the blocks and the metadata
            P.C = C;
            P.vars      = meta.vars;            P.dom       = meta.dom;
            P.space_out = meta.space_out;       P.space_in  = meta.space_in;
            P.dim_out   = meta.dim_out;         P.dim_in    = meta.dim_in;
        end

        function meta = metadata(P)
            % META = METADATA(P) returns the validated metadata of P, for
            % handing to the two-argument constructor.
            meta = struct('vars',{P.vars},'dom',P.dom,...
                'space_out',P.space_out,'space_in',P.space_in,...
                'dim_out',P.dim_out,'dim_in',P.dim_in);
        end
    end
end
