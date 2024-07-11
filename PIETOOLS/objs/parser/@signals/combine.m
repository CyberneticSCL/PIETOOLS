function [out, varargout] = combine(varargin)
out = varargin{1};
for i=2:nargin % find uniques and append into a single vector
    tmp = varargin{i};
    for j=1:length(tmp)
        temp = tmp(j); 
        if ~ismember(temp,out)
            out = [out; temp];
        end
    end
end

varargout = cell(1,nargin);

for i=1:nargin
    varargout{i} = findPermutation(out,varargin{i});
end
end
function P = findPermutation(A,B) % returns P, such that B = P*A
blen = [B.len]; alen = [A.len];
P = zeros(sum(blen),sum(alen));
blen = [0;cumsum(blen)]+1; alen = [0;cumsum(alen)]+1;
[~,idx] = ismember(B,A);
for i=1:length(blen)-1
    tmp = [A.len];
    P(blen(i):blen(i+1)-1,alen(idx(i)):alen(idx(i)+1)-1) = eye(tmp(idx(i)));
end
end