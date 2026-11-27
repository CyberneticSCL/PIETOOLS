function V = polyeval2(X,v1,v2,pts)
% PROVENANCE.  scratchpad/reach1d/polyeval2.m verbatim; op2gal's evaluator.
%
% Evaluate a (matrix-valued) polynomial in at most the two variables v1,v2 on
% a list of points.  pts is npts x 2 giving (v1,v2) at each point.
% Returns nr x nc x npts.  Vectorised over points: the PIETOOLS `subs` route
% is per-point and far too slow for an N^2 grid.
np = size(pts,1);
if isa(X,'double')
    if isempty(X), V = zeros(0,0,np); else, V = repmat(X,1,1,np); end
    return
end
X = polynomial(X);
nr = X.matdim(1); nc = X.matdim(2);
if nr*nc==0 || isempty(X.coefficient), V = zeros(nr,nc,np); return; end
dg = full(X.degmat);  vn = X.varname;  nt = size(dg,1);
M = ones(np,nt);
for k = 1:numel(vn)
    if     strcmp(vn{k},v1), xk = pts(:,1);
    elseif strcmp(vn{k},v2), xk = pts(:,2);
    else,  error('polyeval2: unexpected variable %s',vn{k});
    end
    M = M .* (xk.^(dg(:,k)'));
end
C = full(X.coefficient);          % nt x (nr*nc), column-major over the matrix
V = reshape((M*C).',nr,nc,np);
end
