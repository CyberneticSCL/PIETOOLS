function C = plus(A,B)
% Add two sdopvar operators, preserving affine decision representation.
% Error handling: Checks to ensure A and B are compatible
if any(A.dims~=B.dims)
    error('Dimensions of summands A and B do not match');
end
if any(~strcmp(A.vars.in,B.vars.in)) || any(~strcmp(A.vars.out,B.vars.out))
    error('Summands A and B map between different spaces');
end
if any(any(A.dom.in~=B.dom.in)) || any(any(A.dom.out~=B.dom.out))
    error('Input or output variables in summands A and B have different domains');
end
if numel(A.params.A)~=numel(B.params.A)
    error('number of terms in summands is not equal -- one of them is probably malformed');
end
% First, combine the basis of decision variables Zd of A and B
[A,B,Zd] = CombineDecisionBasis(A,B);
% Then, align monomials in ZL and ZR of A and B, by finding common monomial
% basis, as done for sopvar objects.
[ZR,C1R,C2R] = UnionBasisMonomials(A.ZR,B.ZR);
[ZL,C1L,C2L] = UnionBasisMonomials(A.ZL,B.ZL);
% align sizes, considering the repeated block structures across the 
% input dimension dims(1) and output dimension dims(2) of operators A and B.
C1L = kron(eye(A.dims(1)),C1L); 
C2L = kron(eye(B.dims(1)),C2L);
C1R = kron(eye(A.dims(2)),C1R); 
C2R = kron(eye(B.dims(2)),C2R);
% T_1 = T_1R' \otimes T_1L
T1 = kron(C1R', C1L');
% T_2 = T_2R' \otimes T_2L
T2 = kron(C2R', C2L');
%pa = MatrixMultiply(A.params,C1L',C1R);
%pb = MatrixMultiply(B.params,C2L',C2R);
%params = pa;
% Once C_1, C_2 have compatible sizes and variables,
% C_1(d) + C_2(d) = unvec(A_1 + A_2 + (Bt_1 + Bt_2)*d)
% with the change of monomial basis,
% C_1L C_1(d) C_1R + C_2L C_2(d) C_2R =unvec(T_1 A_1 +T_2 A_2 +(T_1 Bt_1 +T_2 Bt_2)*d)
params_new.A = cell(size(A.params.A));
params_new.B = cell(size(A.params.B));
for i=1:numel(A.params.A)
   % params.A{ii} = params.A{ii}+pb.A{ii};
  %  params.Bt{ii} = params.Bt{ii}+pb.Bt{ii};
    params_new.A{i} = T1*A.params.A{i} + T2*B.params.A{i};
    params_new.B{i} = A.params.B{i}*T1.' + B.params.B{i}*T2.';
end
C = sdopvar(params_new,A.vars,Zd,ZL,ZR,A.dom,A.dims);
end
