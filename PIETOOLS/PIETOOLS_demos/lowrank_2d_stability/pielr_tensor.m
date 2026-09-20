function faceN = pielr_tensor(face1,n,part)
% PIELR_TENSOR  Replicate a certified n=1 face to the n-state member of a
% replicated family, with NO solve at the larger n.
%
%   faceN = pielr_tensor(face1,n)        % partition carried by face1
%   faceN = pielr_tensor(face1,n,part)   % partition supplied explicitly
%
% face1 : a certified base face -- either a cert returned by pielr_certify
%         (fields face, S, Ns, part) or a face struct (fields V, S, Ns, part),
%         e.g. the shipped face_n1.mat.  S must be in ORIGINAL-b units (that is
%         what pielr_certify stores).
% n     : state replication factor (the target system is n identical copies of
%         the base system; the target Gram blocks must be exactly n*Ns).
% part  : 1xB cell, rows-per-monomial-group of each base Gram block.  Computed
%         by pielr_certify (cert.part); pass it only to override.
%
% WHY PER GROUP AND NOT A GLOBAL kron(I_n,V_1).  poslpivar_2d stacks the Gram
% rows GROUP BY GROUP (16 monomial groups) and inside each group as
% kron(eye(n2),Z_g), so rows are state-outer/monomial-inner WITHIN a group and
% the layout is block-diagonal over groups.  A global kron(I_n,V_1) interleaves
% the groups wrongly and is a DIFFERENT subspace: MEASURED, both global
% orderings return X = 0 (op rel 1) where the per-group face certifies at the
% n=1 residual.
%
% The partition is CHECKED against the actual Gram side lengths (sum|Z| = Ns1
% here; pielr_certify then checks n*sum|Z| = Ns of the target program when the
% face is used) -- a wrong partition would build the face on a wrong layout and
% must refuse, not proceed.
%
% OUTPUT faceN: struct with V (per-group replicated face), S = kron(I_n,S1)
% (so pielr_certify can try the no-SDP direct lift first), Ns, n, part, meta.

% accept a pielr_certify cert directly
if isfield(face1,'face') && ~isfield(face1,'V')
    f.V = face1.face;  f.S = face1.S;  f.Ns = face1.Ns;
    if isfield(face1,'part'), f.part = face1.part; end
    face1 = f;
end
assert(isfield(face1,'V') && iscell(face1.V),'pielr_tensor: face1 needs .V (cell)');
B = numel(face1.V);
if nargin<3 || isempty(part)
    assert(isfield(face1,'part') && ~isempty(face1.part), ...
        ['pielr_tensor: no group partition available -- the base certificate ' ...
         'was produced on a program whose Gram is not states x bases (see ' ...
         'pielr_certify notes), or face1 predates cert.part']);
    part = face1.part;
end
assert(iscell(part) && numel(part)==B,'pielr_tensor: part must be 1x%d cell',B);
assert(isscalar(n) && n==round(n) && n>=1,'pielr_tensor: n must be a positive integer');
if isfield(face1,'n') && ~isempty(face1.n)
    assert(face1.n==1,'pielr_tensor: only n=1 base faces tensor (got n=%d)',face1.n);
end

Vn = cell(1,B);
for i = 1:B
    p = part{i}(:)';   V1 = face1.V{i};   r = size(V1,2);
    % the partition must tile the base Gram exactly, or the embedding below
    % lands on the wrong rows
    assert(sum(p)==size(V1,1), ...
        'pielr_tensor: block %d partition sums to %d but the face has %d rows', ...
        i,sum(p),size(V1,1));
    if isfield(face1,'Ns') && ~isempty(face1.Ns)
        assert(face1.Ns(i)==size(V1,1), ...
            'pielr_tensor: block %d face rows %d disagree with Ns=%d', ...
            i,size(V1,1),face1.Ns(i));
    end
    % Vn(:,(a-1)*r+k) places the WHOLE base face vector V1(:,k) into the
    % state-a slot of EVERY active group.  Group g occupies n*p_g consecutive
    % rows at n, state-outer inside the group (kron(eye(n2),Z_g)), so the
    % state-a slot of group g is off_n + (a-1)*p_g + (1:p_g) while the base
    % rows of that group are off_1 + (1:p_g).  Columns are orthonormal for
    % free: different a have disjoint support and each column carries V1(:,k)
    % exactly once.
    Vn{i} = zeros(n*sum(p),n*r);
    o1 = 0;  on = 0;
    for g = 1:numel(p)
        blk = V1(o1+(1:p(g)),:);                    % p_g x r
        for a = 1:n
            Vn{i}(on+(a-1)*p(g)+(1:p(g)),(a-1)*r+(1:r)) = blk;
        end
        o1 = o1 + p(g);   on = on + n*p(g);
    end
end

faceN.V = Vn;
if isfield(face1,'S') && ~isempty(face1.S)
    % states decouple in the replicated family, so the base coefficients
    % replicate too: X_n = Vn * kron(I_n,S_1) * Vn' is (a row permutation of)
    % I_n (x) X_1.  MEASURED: this direct lift verifies at n = 1..6 at the
    % base residual, with no SDP solved.
    faceN.S = cellfun(@(S1)kron(eye(n),S1),face1.S,'uni',0);
end
faceN.Ns = n*cellfun(@(v)size(v,1),face1.V(:)');
faceN.n = n;
faceN.part = part;
faceN.meta = sprintf('pielr_tensor: base face replicated per group to n=%d on %s', ...
                     n,char(datetime('now','Format','yyyy-MM-dd')));
end
