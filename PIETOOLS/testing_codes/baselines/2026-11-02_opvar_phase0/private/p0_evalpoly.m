function V = p0_evalpoly(p, vnames, pts)
% V = P0_EVALPOLY(P,VNAMES,PTS) evaluates the m x q 'polynomial' (or 'double')
% P at the rows of PTS, the values of the variables named in VNAMES, and
% returns the m x q x size(PTS,1) array. An error is raised when P depends
% on a variable not in VNAMES.
np = size(pts,1);
[m,q] = size(p);
if isempty(vnames)
    V = repmat(double(p),1,1,np);
    return
end
v = p0_vars(vnames);
V = zeros(m,q,np);
for k = 1:np
    V(:,:,k) = double(subs(p, v, pts(k,:)'));
end
end
