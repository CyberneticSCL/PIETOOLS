function A = expand_full(a)
% Every 2/3 substitution of the full-domain entries of a multi-index.
%
% alpha_k = 4 denotes a full-domain integral, which is the SUM of the lower
% and upper integrals over one shared kernel, so a term carrying it stands
% for 2^(number of 4s) ordinary semiseparable terms. Direction 1 varies
% fastest, matching the layout of 'alpha_all'.

f = find(a==4);
if isempty(f)
    A = a;
    return
end
nf = numel(f);
A = repmat(a,2^nf,1);
for t = 1:nf
    blk = 2^(t-1);
    A(:,f(t)) = repmat(reshape(repmat([2,3],blk,1),[],1),2^nf/(2*blk),1);
end

end
