% Test_intop_plus_N0.m -- Tests the sum of 'intop' objects without spatial
% variables (N==0: a plain constant/polynomial multiplier, no integration);
% see Test_apply_intop_N0.m for the single-operand case. Sections:
%  (1) sum and difference of such operators (common_intmat with N=0 empty
%      index lists);
%  (2) the 'intop' returned by int(...) of an operator without variables
%      (as in the constant term of piesos_sosvar) must have vars.in and
%      intmat of size 1x0, as isvalid requires of intop(c), so that it can
%      be added to intop(c);
%  (3) the same 1x0 shapes for an 'intop' declared with a vars struct with
%      empty fields {}.
%
% Oracle: the sums are compared with the direct matrix products (K1+K2)*v,
% not with another 'intop' routine.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - Multipliers without spatial variables must be constant (they can still
%   act on a function of the free variable t); kernels are compared as
%   returned by 'kernels'.

tol = 1e-8;
p = 2;
pvar t

%%%% (1) sum and difference of constants
a = intop(2);       b = intop(3);
c = a+b;
if ~isequal(c.matdim,[1,1]) || ~isscalar(c.params) || double(c.params{1})~=5
    error("intop(2)+intop(3) must be the constant 5.")
end
d = a-b;
if double(d.params{1})~=-1
    error("intop(2)-intop(3) must be the constant -1.")
end
v1 = randn(1,p);
Cv1 = apply_op(c,v1);
if ~isdouble(Cv1) || norm(double(Cv1)-5*v1,'fro')>tol*max(1,norm(5*v1,'fro'))
    error("(intop(2)+intop(3)) applied to v must equal 5*v.")
end
% Sum with 0 and with itself
a0 = a+0;       aa = a+a;
if ~isequal(a0.params{1},a.params{1}) || double(aa.params{1})~=4
    error("Sum with 0 or with itself is incorrect.")
end

%%%% (1) matrix-valued multipliers depending on the function variable
K1 = polynomial(randn(2,2));
K2 = polynomial(randn(2,2));
v2 = rand_poly([2,p],t,2,10);
for sgn = [1,-1]
    if sgn==1
        Kop = intop(K1)+intop(K2);      ref = (K1+K2)*v2;
    else
        Kop = intop(K1)-intop(K2);      ref = (K1-K2)*v2;
    end
    Kv = apply_op(Kop,v2);
    err = cleanpoly(Kv-ref,1e-10);
    nrmref = norm(full(ref.C),'fro');
    if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,nrmref)
        error("Sum (sign %d) of matrix multipliers does not match (K1+/-K2)*v.",sgn)
    end
end

%%%% (1) kernels with decision variables
dpvar c1 c2
Kd = intop(c1)+intop(2*c2-1);
Kd_kernels = kernels(Kd);
dif = Kd_kernels{1}-(c1+2*c2-1);
if nnz(dif.C)~=0
    error("Sum of 'dpvar' kernels is incorrect.")
end

%%%% (2) operator without variables obtained by integrating a tensor-PI operator
% The constant function 2 on [0,1], an operator from R to L2[0,1] (as the
% constant blocks of the basis in piesos_sosvar), integrates to the constant 2
pvars = struct('in',{cell(1,0)},'out',{{'s'}});
pdom = struct('in',zeros(0,2),'out',[0,1]);
Kint = int(mat2tensopvar(polynomial(2),pvars,pdom));
if double(Kint.params{1})~=2
    error("int of the constant function 2 on [0,1] must be the constant 2.")
elseif ~isequal(size(Kint.vars.in),[1,0]) || ~isequal(size(Kint.intmat),[1,0])
    error("int of an operator without variables must have 1x0 vars.in and intmat.")
end
Kc = Kint+a;
Cv = apply_op(Kc,v1);
if ~isdouble(Cv) || norm(double(Cv)-4*v1,'fro')>tol*max(1,norm(4*v1,'fro'))
    error("int(2)+intop(2) applied to v must equal 4*v.")
end

%%%% (3) declaration with empty-cell variables {}
K0 = intop({polynomial(7)},{},struct('in',{{}},'out',{cell(1,0)}),zeros(0,2),true);
if ~isequal(size(K0.vars.in),[1,0]) || ~isequal(size(K0.vars.out),[1,0]) || ~isequal(size(K0.intmat),[1,0])
    error("The empty vars and intmat must be stored as 1x0 cells.")
end
K9 = K0+a;
if double(K9.params{1})~=9
    error("intop(7)+intop(2) must be the constant 9.")
end
if ~isequal(K0.vars,a.vars) || ~isequal(K0.intmat,a.intmat)
    error("Constants must have the same vars and intmat.")
end

%%%% Different function spaces must still be rejected
pvar s
Kvar = intop({polynomial(1)},{[0,2]},{'s'},[0,1]);     % int_0^1 v(s) ds
try
    Kvar+a;
    error("Test_intop_plus_N0:noerror","Sum of operators on different function spaces must error.")
catch ME
    if strcmp(ME.identifier,"Test_intop_plus_N0:noerror")
        rethrow(ME)
    end
end

disp('Test_intop_plus_N0 passed: sums of operators without variables are correct.')
