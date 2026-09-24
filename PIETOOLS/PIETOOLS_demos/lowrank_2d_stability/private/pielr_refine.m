function [V,R,q,notes] = pielr_refine(prog,H,D,P,A,V,R,q,vb)               % CC, 09/23/2026
% PIELR_REFINE  Shrink an accepted face one direction at a time, per block.
%
% Transcribed from pielr_certify's local refine_blocks, with the gate taken
% from the adapter so it works in either dimension and for any executive.
% Directions are dropped in order of the certified S_i's own spectrum --
% measured weight, not guessed.  Each trial is a tiny solve.
%
% WHY IT MATTERS FOR A BENCHMARK.  The rank the ladder stops at is an upper
% bound set by the ladder's own step, not by the problem: a run that first
% certifies at [3 3 3] may well hold at [2 3 2], and without the shrink the
% reported rank is really "the first rung that worked".  The rank column is
% the headline number for every low-rank claim this package makes, so it
% should not be an artefact of the search order.
%
% A REJECTED SHRINK IS EVIDENCE, NOT A PROOF.  It says this face does not
% contain a certificate at that width; another face might.  Reported as
% "resists", never as a rank floor.

notes = {};
B = numel(P.rows);
improved = true;
while improved
    improved = false;
    for i = 1:B
        if size(V{i},2) <= 1, continue; end
        Vt = V;   Vt{i} = shrink1(V{i},q,P,i);
        try
            [Rt,qt] = restrict_solve(prog,H,P,D.At,D.b,Vt,false,A);
            if Rt.ok
                V = Vt;  R = Rt;  q = qt;  improved = true;
                if vb, fprintf('  refine blk %d -> r=%d  OK    rel %.3g\n', ...
                        i,size(Vt{i},2),Rt.rel); end
            elseif vb
                fprintf('  refine blk %d -> r=%d  FAILS rel %.3g mineig %.2g  <== resists\n', ...
                        i,size(Vt{i},2),Rt.rel,min(Rt.mineig));
            end
        catch ME
            if vb, fprintf('  refine blk %d ERROR %s\n',i,ME.message(1:min(70,end))); end
        end
    end
end
end

% =========================================================================
function Vi = shrink1(Vi,q,P,i)
% Drop the direction the certified Gram block uses least.  S_i = V_i' X_i V_i
% is tiny, so its eigendecomposition costs nothing at any decision count.
N  = round(sqrt(numel(P.rows{i})));
Xi = reshape(q(P.rows{i}),N,N);   Xi = (Xi+Xi')/2;
Si = Vi'*Xi*Vi;   Si = (Si+Si')/2;
[W,Dg] = eig(Si);
[~,ord] = sort(diag(Dg),'descend');
Vi = Vi*W(:,ord(1:end-1));
% keep the columns orthonormal: V_i*W is orthonormal already when W is, but
% orth() guards against a numerically rank-deficient product
Vi = orth(Vi);
end
