% Test_apply_intop_N0.m -- Tests the N==0 edge case of apply_op: an 'intop'
% with no spatial variables at all (a plain constant multiplier, no
% integration), whose direction-processing loop never runs. Also exercises
% the single-argument constructor intop(params). Oracle: direct
% multiplication c*v and Kmat*v.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - A multiplier without spatial variables must be constant: a polynomial
%   in a free variable t is rejected by the constructor (checked below).
%   The function it acts on may still depend on t.

tol = 1e-8;
p = 2;

%%%% double input to the single-argument constructor
c = randn;
Kop1 = intop(c);
if ~isequal(Kop1.matdim,[1,1])
    error("intop(double): matdim should be [1,1].")
end
v1 = randn(1,p);
Kv1 = apply_op(Kop1,v1);
if ~isdouble(Kv1) || norm(double(Kv1)-c*v1,'fro')>tol*max(1,norm(c*v1,'fro'))
    error("intop(double) applied via apply_op does not match direct scalar multiplication.")
end

%%%% polynomial input to the single-argument constructor
pvar t
Kmat = polynomial([1, 2; 0, 3]);
try
    intop([1+t, 2; 0, t^2]);
    error("intop(polynomial): a kernel depending on a free variable should not be accepted.")
catch ME
    if ~contains(ME.message,"dummy variables")
        rethrow(ME)
    end
end
Kop2 = intop(Kmat);
if ~isequal(Kop2.matdim,[2,2])
    error("intop(polynomial): matdim should be [2,2].")
end
v2 = rand_poly([2,p],t,2,10);
Kv2 = apply_op(Kop2,v2);
ref2 = Kmat*v2;
err2 = cleanpoly(Kv2 - ref2,1e-10);
nrmref2 = norm(full(ref2.C),'fro');     % reference is non-constant: use coeff-matrix norm
if ~isdouble(err2)
    error("intop(polynomial) applied via apply_op: output is not constant.")
elseif norm(double(err2),'fro')>tol*max(1,nrmref2)
    error("intop(polynomial) applied via apply_op: residual exceeds tolerance.")
end

disp('Test_apply_intop_N0 passed: N==0 (plain multiplier) case handled correctly.')
