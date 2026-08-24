function logval = eq(P1,P2,tol)
% Test equality of two sdopvar objects after aligning decision and monomial bases.
if nargin<3
    tol = 1e-14;
end
if ~isa(P1,'sdopvar') && isequal(P1,0)
    logval = eq_zero(P2,tol);
    return
end
if ~isa(P2,'sdopvar') && isequal(P2,0)
    logval = eq_zero(P1,tol);
    return
end
if ~isa(P1,'sdopvar') || ~isa(P2,'sdopvar')
    error('Both inputs must be sdopvar objects, or one input must be zero.');
end

if any(P1.dims~=P2.dims)
    logval = false;
    return
end
if ~isequal(cellstr(string(P1.vars.in(:))),cellstr(string(P2.vars.in(:)))) || ~isequal(cellstr(string(P1.vars.out(:))),cellstr(string(P2.vars.out(:))))
    logval = false;
    return
end
if any(any(P1.dom.in~=P2.dom.in)) || any(any(P1.dom.out~=P2.dom.out))
    logval = false;
    return
end
if numel(P1.params.A)~=numel(P2.params.A)
    logval = false;
    return
end

% Put both objects in the same decision-variable coordinates.
[P1,P2,~] = CombineDecisionBasis(P1,P2);
% Put both objects in the same left and right monomial coordinates.
[~,C1R,C2R] = UnionBasisMonomials(P1.ZR,P2.ZR);
[~,C1L,C2L] = UnionBasisMonomials(P1.ZL,P2.ZL);
C1L = kron(eye(P1.dims(1)),C1L);
C2L = kron(eye(P2.dims(1)),C2L);
C1R = kron(eye(P1.dims(2)),C1R);
C2R = kron(eye(P2.dims(2)),C2R);
T1 = kron(C1R',C1L');
T2 = kron(C2R',C2L');

logval = true;
for ii=1:numel(P1.params.A)
    A1 = T1*P1.params.A{ii};
    A2 = T2*P2.params.A{ii};
    B1 = P1.params.B{ii}*T1.';
    B2 = P2.params.B{ii}*T2.';
    logval = logval && all(abs(A1(:)-A2(:))<tol) && all(abs(B1(:)-B2(:))<tol);
    if ~logval
        return
    end
end
end

function logval = eq_zero(P,tol)
if ~isa(P,'sdopvar')
    logval = isequal(P,0);
    return
end
logval = true;
for ii=1:numel(P.params.A)
    logval = logval && all(abs(P.params.A{ii}(:))<tol) && all(abs(P.params.B{ii}(:))<tol);
    if ~logval
        return
    end
end
end
