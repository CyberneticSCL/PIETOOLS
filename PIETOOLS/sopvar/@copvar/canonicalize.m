function outcopvar = canonicalize(incopvar)
% canonicalize
%
% Build and store the two auxiliary common-basis representations of a
% copvar:
%
%   1. row-common form:
%        all blocks in row i have the same left basis ZL;
%
%   2. column-common form:
%        all blocks in column j have the same right basis ZR.
%
% These representations are independent. No attempt is made to put the
% blocks simultaneously on common left and right bases.
%
% The native block representation incopvar.C is preserved unchanged.
%
% The auxiliary representations are stored in
%
%   outcopvar.leftCommonBasis
%   outcopvar.leftCommonC
%
%   outcopvar.rightCommonBasis
%   outcopvar.rightCommonC
%
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  inmopvar -> incopvar, outmopvar -> outcopvar. Moved from
%                  @mopvar/ with the class. Only identifiers changed; the
%                  routine is otherwise as committed in 0f4cb3d3.

outcopvar = incopvar;

% ---------------------------------------------------------------
% Build row-common representation.
%
% commonRowBasis returns a copvar whose C blocks have been rewritten
% onto common left bases, while leaving their right bases unchanged.
% ---------------------------------------------------------------
L = commonRowBasis(incopvar);

outcopvar.leftCommonBasis = L.leftCommonBasis;
outcopvar.leftCommonC     = L.C;

% ---------------------------------------------------------------
% Build column-common representation.
%
% IMPORTANT: start again from incopvar, not from L. This keeps the
% column-common representation independent of the row-common one.
% ---------------------------------------------------------------
R = commonColumnBasis(incopvar);

outcopvar.rightCommonBasis = R.rightCommonBasis;
outcopvar.rightCommonC     = R.C;

end
function Q = commonRowBasis(P)
% commonRowBasis
%
% Q = commonRowBasis(P) rewrites each block row of the copvar P so that
% all populated blocks in row i use the same left monomial basis ZL.
%
% For
%
%   P.C{i,j} = ZL_ij(s)' * C_ij * ZR_ij(theta),
%
% the result has
%
%   Q.C{i,j} = ZL_i(s)' * Cbar_ij * ZR_ij(theta),
%
% where ZL_i is a common basis containing the left basis of every
% populated block in row i.
%
% The represented operator is unchanged. Only zero rows are inserted
% into the coefficient matrices.
%
% The common basis for row i is stored in
%
%   Q.leftCommonBasis{i}.
%
% Row consistency is already guaranteed by the copvar constructor, so
% no checks on vars.out, dom.out, or dim_out are performed here.

[M,N] = size(P.C);

Cnew = P.C;
Zcommon = cell(M,1);

for i = 1:M

    % Find the populated blocks in this row.
    jj = find(~cellfun(@isempty,P.C(i,:)));

    if isempty(jj)
        % Structurally zero row. There is no block basis to synchronize.
        Zcommon{i} = {};
        continue
    end

    % -------------------------------------------------------------
    % Construct the common left basis.
    %
    % All blocks in the row have the same output variable SET, but
    % their vars.out ordering may differ because sopvar stores
    %
    %   vars.out = [S2,S3].
    %
    % Therefore bases are matched by variable name, not by position.
    % -------------------------------------------------------------
    rowVars = P.vars(P.space_out(i,:));
    Zi = cell(1,numel(rowVars));

    for k = 1:numel(rowVars)
        v = rowVars{k};
        zk = [];

        for j = jj
            B = P.C{i,j};

            pos = find(strcmp(B.vars.out,v),1);
            zk = union(zk,B.ZL{pos}(:),'sorted');
        end

        Zi{k} = zk;
    end

    Zcommon{i} = Zi;

    % Number of monomials in the common tensor-product basis.
    nNew = prod(cellfun(@numel,Zi));

    % -------------------------------------------------------------
    % Rewrite every populated block onto the common basis.
    % -------------------------------------------------------------
    for j = jj
        B = P.C{i,j};

        % The common basis has to be put into this block's vars.out
        % ordering before constructing the coefficient embedding.
        ZLnew = cell(size(B.ZL));

        for k = 1:numel(B.vars.out)
            pos = find(strcmp(rowVars,B.vars.out{k}),1);
            ZLnew{k} = Zi{pos};
        end

        % Map each old tensor-product monomial to its location in ZLnew.
        idx = basisEmbeddingIndex(B.ZL,ZLnew);

        nOld = numel(idx);
        p = B.dims(1);

        % Lift the monomial embedding over the p output components.
        %
        % Coefficient rows are arranged as
        %
        %   component 1: all ZL monomials
        %   component 2: all ZL monomials
        %   ...
        %
        rowOld = (1:p*nOld).';
        comp = floor((rowOld-1)/nOld);
        mono = mod(rowOld-1,nOld) + 1;

        rowNew = comp*nNew + idx(mono);

        params = B.params;

        for a = 1:numel(params)
            A = params{a};

            % Preserve sopvar's zero-block shorthand.
            if isempty(A) || (isscalar(A) && A==0)
                continue
            end

            Anew = sparse(p*nNew,size(A,2));
            Anew(rowNew,:) = A;

            params{a} = Anew;
        end

        % Reconstruct through sopvar so its class invariants remain valid.
        Cnew{i,j} = sopvar(params,B.vars,ZLnew,B.ZR,B.dom,B.dims);
    end
end

% Operator spaces and dimensions have not changed, so reuse P's metadata.
Q = copvar(Cnew,P.metadata());

Q.leftCommonBasis = Zcommon;
end
function Q = commonColumnBasis(P)
% commonColumnBasis
%
% Q = commonColumnBasis(P) rewrites each block column of the copvar P so
% that all populated blocks in column j use the same right monomial basis.
%
% For
%
%   P.C{i,j} = ZL_ij(s)' * C_ij * ZR_ij(theta),
%
% the result has
%
%   Q.C{i,j} = ZL_ij(s)' * Cbar_ij * ZR_j(theta),
%
% where ZR_j is a common basis containing the right basis of every
% populated block in column j.
%
% The represented operator is unchanged. Only zero columns are inserted
% into the coefficient matrices.
%
% The common basis for column j is stored in
%
%   Q.rightCommonBasis{j}.
%
% Column consistency is already guaranteed by the copvar constructor, so
% no checks on vars.in, dom.in, or dim_in are performed here.

[M,N] = size(P.C);

Cnew = P.C;
Zcommon = cell(N,1);

for j = 1:N

    % Find the populated blocks in this column.
    ii = find(~cellfun(@isempty,P.C(:,j)));

    if isempty(ii)
        % Structurally zero column.
        Zcommon{j} = {};
        continue
    end

    % -------------------------------------------------------------
    % Construct the common right basis.
    %
    % As on the left side, use variable names rather than assuming
    % identical positional ordering inside the individual sopvar blocks.
    % -------------------------------------------------------------
    colVars = P.vars(P.space_in(j,:));
    Zj = cell(1,numel(colVars));

    for k = 1:numel(colVars)
        v = colVars{k};
        zk = [];

        for i = ii
            B = P.C{i,j};

            pos = find(strcmp(B.vars.in,v),1);
            zk = union(zk,B.ZR{pos}(:),'sorted');
        end

        Zj{k} = zk;
    end

    Zcommon{j} = Zj;

    nNew = prod(cellfun(@numel,Zj));

    % -------------------------------------------------------------
    % Rewrite every populated block onto the common basis.
    % -------------------------------------------------------------
    for i = ii
        B = P.C{i,j};

        % Put the column-common basis into this block's vars.in order.
        ZRnew = cell(size(B.ZR));

        for k = 1:numel(B.vars.in)
            pos = find(strcmp(colVars,B.vars.in{k}),1);
            ZRnew{k} = Zj{pos};
        end

        % Map old tensor-product monomials into the new basis.
        idx = basisEmbeddingIndex(B.ZR,ZRnew);

        nOld = numel(idx);
        q = B.dims(2);

        % Lift the monomial embedding over the q input components.
        colOld = (1:q*nOld).';
        comp = floor((colOld-1)/nOld);
        mono = mod(colOld-1,nOld) + 1;

        colNew = comp*nNew + idx(mono);

        params = B.params;

        for a = 1:numel(params)
            A = params{a};

            if isempty(A) || (isscalar(A) && A==0)
                continue
            end

            Anew = sparse(size(A,1),q*nNew);
            Anew(:,colNew) = A;

            params{a} = Anew;
        end

        Cnew{i,j} = sopvar(params,B.vars,B.ZL,ZRnew,B.dom,B.dims);
    end
end

Q = copvar(Cnew,P.metadata());

Q.rightCommonBasis = Zcommon;
end

function idx = basisEmbeddingIndex(Zold,Znew)
% basisEmbeddingIndex
%
% idx(k) is the position of tensor-product monomial k of Zold in Znew.
%
% Zold and Znew are cell arrays of one-dimensional exponent lists:
%
%   Z = Z{1} kron Z{2} kron ... kron Z{n}.
%
% The last variable therefore varies fastest.

nv = numel(Zold);

if nv==0
    idx = 1;
    return
end

[tf,idx] = ismember(Zold{1}(:),Znew{1}(:));
if ~all(tf)
    error('copvar:basisEmbedding','Old basis is not contained in new basis.')
end

for k = 2:nv
    [tf,ik] = ismember(Zold{k}(:),Znew{k}(:));

    if ~all(tf)
        error('copvar:basisEmbedding', ...
            'Old basis is not contained in new basis.')
    end

    nk = numel(Znew{k});

    idx = reshape((idx(:)-1)*nk + ik(:).',[],1);
end
end