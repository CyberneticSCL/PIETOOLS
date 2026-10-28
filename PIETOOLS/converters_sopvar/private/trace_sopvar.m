function Q = trace_sopvar(P,name,bndry)
% Q = TRACE_SOPVAR(P,NAME,BNDRY) evaluates the output of the PI operator P
% at a boundary of the variable NAME. That is, for the operator
%   (Pu)(s,r) = sum int K(s,r,t,...) u(t,...) dt,
% with BNDRY='lower' or 'upper' we return the operator Q such that
%   (Qu)(r) = (Pu)(a,r)    or    (Qu)(r) = (Pu)(b,r),
% for s in [a,b].
%
% INPUTS
% - P:      'sopvar' object;
% - name:   Name of the output variable s at which to evaluate;
% - bndry:  'lower' or 'upper';
%
% OUTPUT
% - Q:      'sopvar' object, not having s as an output variable. If s is an
%           input variable of P as well, then Q integrates over s.
%
% NOTES
% - If s is both an input and an output variable, the operator involves a
%   multiplier and partial integrals in s. At the boundary, the integrals
%   reduce to full integrals, but the multiplier would result in a point
%   evaluation of the input, which is not a PI operator. An error is
%   therefore thrown if P has a non-zero multiplier in s;
% - The variable must be an output variable of P; otherwise, P is returned.

iv = find(strcmp(P.vars.out,name));
if isempty(iv)
    Q = P;
    return
end
if strcmp(bndry,'lower'),   val = P.dom.out(iv,1);
elseif strcmp(bndry,'upper'),   val = P.dom.out(iv,2);
else,   error('Boundary should be specified as ''lower'' or ''upper''.')
end

nz = cellfun(@numel,P.ZL);
w = val.^P.ZL{iv}(:);                   % monomials evaluated at the boundary
q = P.dims(1);
params = P.params;

S3 = intersect(P.vars.out,P.vars.in);   % sorted shared variables
n3 = numel(S3);
if any(strcmp(S3,name))
    % The variable is shared: evaluate partial integrals at the boundary
    d = find(strcmp(S3,name));
    sz = size(params);
    if n3==1,   sz = [3,1];     end
    A = reshape(params,sz);
    perm = [d,setdiff(1:numel(sz),d)];
    A = permute(A,perm);
    sz_p = size(A);
    A = reshape(A,sz_p(1),[]);
    if any(cellfun(@(c) nnz(c)~=0,A(1,:)))
        error(['Cannot evaluate the output of the operator at the boundary, as it involves ',...
                'a multiplier along the evaluated variable.'])
    end
    if strcmp(bndry,'lower'),   keep = 3;
    else,   keep = 2;
    end
    params = A(keep,:);
    if n3-1==0,         params = params(1);
    elseif n3-1==1,     params = reshape(params,[3,1]);
    else,               params = reshape(params,3*ones(1,n3-1));
    end
end
for m=1:numel(params)
    params{m} = contract_rows(full(params{m}),nz,iv,q,w);
end

vars = P.vars;      dom = P.dom;        ZL = P.ZL;
vars.out(iv) = [];  dom.out(iv,:) = [];     ZL(iv) = [];
if isempty(vars.out),   vars.out = {};      end
Q = sopvar(params,vars,ZL,P.ZR,dom,P.dims);
end


function C = contract_rows(C,nz,iv,q,w)
% Contract the monomials associated to the iv-th output variable, in the
% rows of the coefficient matrix C, with the vector w. The rows of C are
% ordered by component, then by monomial with the last variable fastest.
N = numel(nz);
ncol = size(C,2);
Ct = reshape(C,[fliplr(nz),q,ncol]);
if N==1,    Ct = reshape(C,[nz,q,ncol]);    end
dimpos = N-iv+1;
others = setdiff(1:ndims(Ct),dimpos);
Ct = permute(Ct,[dimpos,others]);
Ct = reshape(Ct,nz(iv),[]);
Cn = w(:)'*Ct;
nz_new = nz;    nz_new(iv) = [];
C = reshape(Cn,[prod(nz_new)*q,ncol]);
end
