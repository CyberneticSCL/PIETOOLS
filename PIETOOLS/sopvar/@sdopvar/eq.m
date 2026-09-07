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

% Put both objects in the same decision-variable and monomial coordinates.  % MMP, 09/07/2026
[ops,T] = sync_basis({P1,P2});                                              % MMP, 09/07/2026
P1 = ops{1};    P2 = ops{2};                                                % MMP, 09/07/2026
%[P1,P2,~] = CombineDecisionBasis(P1,P2);                                  % MMP, 09/07/2026 (was)
%[~,C1R,C2R] = UnionBasisMonomials(P1.ZR,P2.ZR);                           % MMP, 09/07/2026 (was)
%[~,C1L,C2L] = UnionBasisMonomials(P1.ZL,P2.ZL);                           % MMP, 09/07/2026 (was)
%C1L = kron(eye(P1.dims(1)),C1L);                                          % MMP, 09/07/2026 (was)
%C2L = kron(eye(P2.dims(1)),C2L);                                          % MMP, 09/07/2026 (was)
%C1R = kron(eye(P1.dims(2)),C1R);                                          % MMP, 09/07/2026 (was)
%C2R = kron(eye(P2.dims(2)),C2R);                                          % MMP, 09/07/2026 (was)
%T1 = kron(C1R',C1L');                                                     % MMP, 09/07/2026 (was)
%T2 = kron(C2R',C2L');                                                     % MMP, 09/07/2026 (was)

logval = true;
for ii=1:numel(P1.params.A)
    [A1,B1] = apply_basis_map(T{1},P1.params.A{ii},P1.params.B{ii});        % MMP, 09/07/2026
    [A2,B2] = apply_basis_map(T{2},P2.params.A{ii},P2.params.B{ii});        % MMP, 09/07/2026
%   A1 = T1*P1.params.A{ii};                                               % MMP, 09/07/2026 (was)
%   A2 = T2*P2.params.A{ii};                                               % MMP, 09/07/2026 (was)
%   B1 = P1.params.B{ii}*T1.';                                             % MMP, 09/07/2026 (was)
%   B2 = P2.params.B{ii}*T2.';                                             % MMP, 09/07/2026 (was)
    % 'all(abs(D(:))<tol)' on a sparse D materialises a FULL logical of     % MMP, 09/07/2026
    % q*nC entries, since every implicit zero satisfies the comparison:     % MMP, 09/07/2026
    % 1166 MB and 0.21 s measured at q=1e5, nC=1296, to compare a 0.1 MB    % MMP, 09/07/2026
    % object. Testing the stored nonzeros is equivalent for tol>0 and is    % MMP, 09/07/2026
    % O(nnz), and the rows of B are the decision variables.                 % MMP, 09/07/2026
    logval = logval && ~any(abs(nonzeros(A1-A2))>=tol) ...
                    && ~any(abs(nonzeros(B1-B2))>=tol);                     % MMP, 09/07/2026
%   logval = logval && all(abs(A1(:)-A2(:))<tol) && all(abs(B1(:)-B2(:))<tol); % MMP, 09/07/2026 (was)
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
    % Same sparse-comparison trap as above; see the note in eq.           % MMP, 09/07/2026
    logval = logval && ~any(abs(nonzeros(P.params.A{ii}))>=tol) ...
                    && ~any(abs(nonzeros(P.params.B{ii}))>=tol);            % MMP, 09/07/2026
%   logval = logval && all(abs(P.params.A{ii}(:))<tol) && all(abs(P.params.B{ii}(:))<tol); % MMP, 09/07/2026 (was)
    if ~logval
        return
    end
end
end
