function [objs,T,Zd,ZL,ZR] = sync_basis(objs,shared_Zd)                     % MMP, 09/26/2026
% function [objs,T,Zd,ZL,ZR] = sync_basis(objs)                             % MMP, 09/26/2026 (was)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [objs,T,Zd,ZL,ZR] = sync_basis(objs) puts N sdopvar objects on a common
% decision variable list and a common pair of monomial bases, in one pass
% over all N rather than N-1 pairwise merges.
%
% INPUTS
% objs: 1 x N cell of 'sdopvar' objects. They must already agree on vars and
%       dom; this routine reconciles only Zd, ZL and ZR;
% shared_Zd: (optional) true is the CALLER's guarantee that each operand's  % MMP, 09/26/2026
%       Zd holds the same entries as objs{1}.Zd, e.g. because all were      % MMP, 09/26/2026
%       built from one list. The O(q)-per-operand comparison is skipped.    % MMP, 09/26/2026
%       Default false;                                                      % MMP, 09/26/2026
%
% OUTPUTS
% objs: the same operators, with the rows of params.B reordered onto the
%       merged decision variable list Zd (zero rows for variables the
%       operand does not use);
% T:    1 x N cell. T{k} maps vec(C) of operand k onto the merged monomial
%       bases, so the merged coefficient of parameter i is
%           T{k}*objs{k}.params.A{i}    and    objs{k}.params.B{i}*T{k}.'
%       T{k} is EMPTY whenever that map is the identity, so a caller can
%       skip the multiply instead of paying for a no-op;
% Zd:   merged decision variable list;
% ZL,ZR: merged monomial bases, as column vectors per direction;
%
% NOTES
% Replaces the prologue that @sdopvar/plus, eq, horzcat and vertcat each
% carried inline, and makes it N-ary. Three things are gained over calling
% the pairwise version N-1 times:
%
%   1. The decision variable lists are merged with a single
%      unique(...,'stable') over the concatenation of all N: the same list
%      as iterating Zd = [Zd; setdiff(Zd_k,Zd,...)], since both keep first
%      occurrences in order, at one sort of the whole set instead of N-1
%      setdiffs against a growing accumulator (the dominant pairwise cost).
%   2. Each operand is embedded into the FINAL merged basis exactly once.
%      Chaining pairwise merges re-maps operand 1 through N-1 successive
%      transforms, and each of those acts on the growing accumulated
%      operator rather than the original operand.
%   3. An operand already on the merged bases gets T{k} = [], so the
%      identity multiply is skipped rather than performed.
%
% The bases are unioned per direction with 'union' and forced to columns.
% MATLAB's 'union' returns a ROW when both arguments are rows, and a monomial
% basis that comes back row-oriented breaks plus, minus and eq silently; see
% the CANONICAL MULTIPLIER FORM note in sdopvar.m.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
% Initial coding MMP, 09/07/2026
% MMP, 09/26/2026: (1) Optional input 'shared_Zd', the caller's guarantee
%                  that the decision lists are equal, skips the O(q) isequal
%                  per operand (1.24 of 1.36 s in the 2-D container Hinf
%                  build, 147 calls on 3.8e5 names); 'plus_batch' forwards
%                  it. (2) Zd is compared as stored when already a column:
%                  'P.Zd(:)' copies the q-entry cell. Same output.
% MMP, 09/29/2026: The positions of each operand's monomials in the merged
%                  bases are read directly by the local 'basis_positions'
%                  (ismember per direction, Kronecker strides) instead of
%                  through the selection matrices of 'UnionBasisMonomials'
%                  and find(); same positions, flat in q (3-D heat build:
%                  1.83 -> 0.62 s). A monomial missing from the merged
%                  basis now errors. The @sdopvar copy of
%                  'UnionBasisMonomials' (MMP, SS, AT, DJ, 06/08/2026) lost
%                  its only caller and is in the dead-code folder.
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): log entries
%                  condensed, a measurement dropped from the NOTES, the
%                  BEGIN note shortened; code unchanged.

N = numel(objs);
if N==0
    T = {};  Zd = cell(0,1);  ZL = {};  ZR = {};   return
end

% % % Decision variables: one merged list for all N.
% The common case is that every operand already carries the same list, and
% then no comparison beyond isequal is needed.
% Zd = objs{1}.Zd(:);                                                       % MMP, 09/26/2026 (was)
% '(:)' on a property copies the whole cell (O(q)); a column needs none.    % MMP, 09/26/2026
Zd = objs{1}.Zd;                                                            % MMP, 09/26/2026
if ~iscolumn(Zd),   Zd = Zd(:);     end                                     % MMP, 09/26/2026
allsame = true;
% Equal by the caller's guarantee: no O(q) comparison needed. A length      % MMP, 09/26/2026
% mismatch is a broken guarantee, and merging on it would pair the rows of  % MMP, 09/26/2026
% params.B with the wrong names, so it errors rather than proceeding.       % MMP, 09/26/2026
if nargin>=2 && shared_Zd                                                   % MMP, 09/26/2026
    for k = 2:N                                                             % MMP, 09/26/2026
        if numel(objs{k}.Zd)~=numel(Zd)                                     % MMP, 09/26/2026
            error("'shared_Zd' was set, but operand %d has %d decision "+...
                  "variables against %d.",k,numel(objs{k}.Zd),numel(Zd))    % MMP, 09/26/2026
        end                                                                 % MMP, 09/26/2026
    end                                                                     % MMP, 09/26/2026
else                                                                        % MMP, 09/26/2026
for k = 2:N
%   if ~isequal(objs{k}.Zd(:),Zd),    allsame = false;    break;    end     % MMP, 09/26/2026 (was)
    Zk = objs{k}.Zd;                                                        % MMP, 09/26/2026
    if ~iscolumn(Zk),   Zk = Zk(:);     end                                 % MMP, 09/26/2026
    if ~isequal(Zk,Zd),    allsame = false;    break;    end                % MMP, 09/26/2026
end
end                                                                         % MMP, 09/26/2026
if ~allsame
    lists = cell(1,N);
    for k = 1:N,    lists{k} = objs{k}.Zd(:);    end
    Zd = unique(vertcat(lists{:}),'stable');
end
if ~allsame
    % Only the rows of params.B move, so an operand whose list is already Zd
    % is left alone by ChangeDecVar's own isequal guard.
    for k = 1:N
        objs{k} = ChangeDecVar(objs{k},Zd);
    end
else
    for k = 1:N,    objs{k}.Zd = Zd;    end
end

% % % Monomial bases: union over all N, one direction at a time.
ZL = objs{1}.ZL;    ZR = objs{1}.ZR;
for k = 2:N
    for i = 1:numel(ZL)
        ZL{i} = reshape(union(ZL{i},objs{k}.ZL{i}),[],1);
    end
    for i = 1:numel(ZR)
        ZR{i} = reshape(union(ZR{i},objs{k}.ZR{i}),[],1);
    end
end

% % % One embedding per operand, straight onto the merged bases.
NLnew = prod([cellfun(@numel,ZL),1]);
NRnew = prod([cellfun(@numel,ZR),1]);
T = cell(1,N);
for k = 1:N
    if isequal(objs{k}.ZL,ZL) && isequal(objs{k}.ZR,ZR)
        % Already on the merged bases: leave T{k} empty so the caller skips
        % the remapping entirely rather than performing an identity one.
        continue
    end
    % BEGIN MMP, 09/29/2026: aL(a) is the position in ZL of monomial a of
    % objs{k}.ZL, read directly; likewise aR. The embedding into the merged
    % basis is a pure scatter, stored as the index vector T{k}.idx below
    % rather than as a matrix; see 'apply_basis_map'.
%   [~,CL] = UnionBasisMonomials(objs{k}.ZL,ZL);                            % MMP, 09/29/2026 (was)
%   [~,CR] = UnionBasisMonomials(objs{k}.ZR,ZR);                            % MMP, 09/29/2026 (was)
%   aL = zeros(size(CL,1),1);   [rL,cL] = find(CL);     aL(rL) = cL;        % MMP, 09/29/2026 (was)
%   aR = zeros(size(CR,1),1);   [rR,cR] = find(CR);     aR(rR) = cR;        % MMP, 09/29/2026 (was)
    aL = basis_positions(objs{k}.ZL,ZL);                                    % MMP, 09/29/2026
    aR = basis_positions(objs{k}.ZR,ZR);                                    % MMP, 09/29/2026
    % END MMP, 09/29/2026
    mk = objs{k}.dims(1);       nk = objs{k}.dims(2);

    % Entry (matrix row i, ZL index a, matrix column j, ZR index b) sits at
    % vec position (i-1)*NLk+a + ((j-1)*NRk+b-1)*mk*NLk, the monomial index
    % inner on both axes. Its destination is the same expression with NLk,
    % NRk replaced by the merged sizes and a,b by aL(a),aR(b). Building
    % rowmap and colmap first keeps this one reshape per axis.
    rowmap = reshape(aL + (0:mk-1)*NLnew,[],1);
    colmap = reshape(aR + (0:nk-1)*NRnew,[],1);
    T{k} = struct('idx',reshape(rowmap+(colmap.'-1)*(mk*NLnew),[],1), ...
                  'nC', mk*NLnew*nk*NRnew);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 09/29/2026
function a = basis_positions(Zs,Zf)                                         % MMP, 09/29/2026
% A = BASIS_POSITIONS(ZS,ZF): A(j) is the position in the monomial vector   % MMP, 09/29/2026
% kron(ZF{1},...,ZF{N}) of monomial j of kron(ZS{1},...,ZS{N}), first       % MMP, 09/29/2026
% variable outermost, as in 'UnionBasisMonomials' and the class bases.      % MMP, 09/29/2026
% Requires ZS{i} within ZF{i} for each direction i, which holds for the     % MMP, 09/29/2026
% unions formed in 'sync_basis'; otherwise it errors, since a missing       % MMP, 09/29/2026
% monomial would get position 0. Cost: the monomial axis only, not q.       % MMP, 09/29/2026
a = 1;                                                                      % MMP, 09/29/2026
for i = 1:numel(Zs)                                                         % MMP, 09/29/2026
    [tf,p] = ismember(Zs{i}(:),Zf{i}(:));                                   % MMP, 09/29/2026
    if ~all(tf)                                                             % MMP, 09/29/2026
        error('sync_basis:basisNotContained',['A monomial of an operand '...
              'basis is missing from the merged basis.'])                   % MMP, 09/29/2026
    end                                                                     % MMP, 09/29/2026
    % (a-1)*numel(Zf{i}) + p is the merged position of (earlier directions, % MMP, 09/29/2026
    % direction i); direction i is the inner index, hence the transpose.    % MMP, 09/29/2026
    a = reshape(((a(:)-1)*numel(Zf{i}) + p(:).').',[],1);                   % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026
