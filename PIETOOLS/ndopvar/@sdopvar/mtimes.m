function C = mtimes(A,B)
% Multiplication for sdopvar with numeric matrices/scalars.
if isa(A,'sdopvar') && isnumeric(B)
    if isscalar(B)
        C = A;
        for ii=1:numel(C.params.A)
            C.params.A{ii} = C.params.A{ii}*B;
            C.params.B{ii} = C.params.B{ii}*B;
        end
        return
    end
    if A.dims(2)~=size(B,1)
        error('Inner dimensions do not match.');
    end
    % If B is matrix
    left_size = A.dims(1)*prod(cellfun(@numel,A.ZL));
    R = kron(B,eye(prod(cellfun(@numel,A.ZR))));
    params = MatrixMultiply(A.params,speye(left_size),R);
    C = sdopvar(params,A.vars,A.Zd,A.ZL,A.ZR,A.dom,[A.dims(1),size(B,2)]);
elseif isnumeric(A) && isa(B,'sdopvar')
    if isscalar(A)
        C = B;
        for ii=1:numel(C.params.A)
            C.params.A{ii} = A*C.params.A{ii};
            C.params.B{ii} = A*C.params.B{ii};
        end
        return
    end
    if size(A,2)~=B.dims(1)
        error('Inner dimensions do not match.');
    end
    % If A is matrix
    right_size = B.dims(2)*prod(cellfun(@numel,B.ZR));
    L = kron(A,eye(prod(cellfun(@numel,B.ZL))));
    params = MatrixMultiply(B.params,L,speye(right_size));
    C = sdopvar(params,B.vars,B.Zd,B.ZL,B.ZR,B.dom,[size(A,1),B.dims(2)]);
else
    error('Composition of sdopvar objects is not supported yet.');
end
end
