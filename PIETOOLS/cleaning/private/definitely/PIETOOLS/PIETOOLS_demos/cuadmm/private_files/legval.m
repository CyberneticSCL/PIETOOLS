function P = legval(x,a,b,D)
% L2[a,b]-ORTHONORMAL shifted Legendre polynomials evaluated at points x, by
% the Bonnet recurrence.  Evaluating by recurrence (rather than constructing
% the polynomials in the monomial basis and substituting) is what keeps this
% conditioned: the symbolic route lost 7 digits by degree 8.
x = x(:);  z = (2*x-a-b)/(b-a);
P = zeros(numel(x),D+1);
p0 = ones(numel(x),1);
P(:,1) = p0/sqrt(b-a);
if D>=1
    p1 = z;
    P(:,2) = p1*sqrt(3/(b-a));
end
for k = 2:D
    p2 = ((2*k-1)*z.*p1 - (k-1)*p0)/k;
    P(:,k+1) = p2*sqrt((2*k+1)/(b-a));
    p0 = p1;  p1 = p2;
end
end
