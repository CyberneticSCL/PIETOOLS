function [out, varargout] = combine(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that combines a vector of states, x, to get (x_min,A),
% where x_min is the unique set of vectors in, x, and x = A*x_min
% Input: 
% varargin - states to combine, x{1}...x{N}
% Output:
% out - vector of unique states, x{1}...x{M}
% varargout - Matrices A{1}..A{M}, where x{i} = A{i}*out

out = varargin{1};
for i=2:nargin
    tmp = varargin{i};
    for j=1:length(tmp)
        s.type = '()'; s.subs = {j};
        temp = subsref(tmp,s); 
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
s.type = '.'; s.subs = 'veclength';
P = zeros(sum(subsref(B,s)),sum(subsref(A,s)));
blen = [0;cumsum(subsref(B,s))]+1; alen = [0;cumsum(subsref(A,s))]+1;
[~,idx] = ismember(B,A);
for i=1:length(blen)-1
    tmp = subsref(A,s);
    P(blen(i):blen(i+1)-1,alen(idx(i)):alen(idx(i)+1)-1) = eye(tmp(idx(i)));
end
end