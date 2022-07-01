function [out, varargout] = combine(varargin)
out = varargin{1};
for i=2:nargin
    tmp = varargin{i};
    for j=1:length(tmp.veclength)
        s.type = '()'; s.subs = {j};
        temp = subsref(tmp,s); % some weird error if temp is not used
        if ~ismember(temp.statename,out.statename)
            out = [out; temp];
        end
    end
end

varargout = cell(1,length(out.statename));

for i=1:nargin
    varargout{i} = findPermutation(out,varargin{i});
end
end
function P = findPermutation(A,B)
P = zeros(length(B),length(A));
blen = [0,cumsum(B.veclength)]+1; alen = [0,cumsum(A.veclength)]+1;
[~,idx] = ismember(B.statename, A.statename);
for i=1:length(blen)-1
    P(blen(i):blen(i+1)-1,alen(idx(i)):alen(idx(i)+1)-1) = eye(A.veclength(idx(i)));
end
end