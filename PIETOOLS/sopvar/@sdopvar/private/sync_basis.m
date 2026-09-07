function [objs,T,Zd,ZL,ZR] = sync_basis(objs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [objs,T,Zd,ZL,ZR] = sync_basis(objs) puts N sdopvar objects on a common
% decision variable list and a common pair of monomial bases, in one pass
% over all N rather than N-1 pairwise merges.
%
% INPUT
% objs: 1 x N cell of 'sdopvar' objects. They must already agree on vars and
%       dom; this routine reconciles only Zd, ZL and ZR;
%
% OUTPUT
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
%      unique(...,'stable') over the concatenation of all N. That is
%      provably the same list as iterating Zd = [Zd; setdiff(Zd_k,Zd,...)],
%      which is the order 'CombineDecisionBasis' produces, because both keep
%      first occurrences in order; but it costs one sort of the whole set
%      instead of N-1 setdiffs against a growing accumulator. This is the
%      dominant cost: a setdiff over decision variable NAMES accounted for
%      roughly 90% of pairwise 'plus'/'horzcat' time at q = 21504.
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

N = numel(objs);
if N==0
    T = {};  Zd = cell(0,1);  ZL = {};  ZR = {};   return
end

% % % Decision variables: one merged list for all N.
% The common case is that every operand already carries the same list, and
% then no comparison beyond isequal is needed.
Zd = objs{1}.Zd(:);
allsame = true;
for k = 2:N
    if ~isequal(objs{k}.Zd(:),Zd),    allsame = false;    break;    end
end
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
    % ZL contains objs{k}.ZL, so the union below IS ZL and CL embeds
    % operand k into it; likewise CR. Reusing 'UnionBasisMonomials' keeps
    % the monomial ordering identical to the pairwise code it replaces.
    [~,CL] = UnionBasisMonomials(objs{k}.ZL,ZL);
    [~,CR] = UnionBasisMonomials(objs{k}.ZR,ZR);

    % CL and CR are selection matrices: exactly one 1 per row, in the
    % column this operand's monomial occupies in the merged basis. So
    % kron(CR',CL') would be an nC_new x nC_k sparse holding one 1 per
    % column and nothing else -- a pure scatter. Read the positions off CL
    % and CR and store the scatter as an index vector instead of building
    % that matrix; see 'apply_basis_map'.
    aL = zeros(size(CL,1),1);   [rL,cL] = find(CL);     aL(rL) = cL;
    aR = zeros(size(CR,1),1);   [rR,cR] = find(CR);     aR(rR) = cR;
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
