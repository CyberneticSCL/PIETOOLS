% Test_tensopmat_native_blocks.m -- tests blocks of a 'tensopmat' of type
% 'double', 'polynomial' or 'dpvar', acting on R (see @tensopmat/
% tensopmat.m): the block M maps x in R^n to the function M*x in L2^m[S],
% for S the variables that appear in M.
% Parts:
%   (1) construction: dimensions, variables, domains and depmat of the
%       object, from explicit blocks, the domains from the constructor,
%       from other blocks, or from the object itself, and the errors;
%   (2) isvalid, for valid and for invalid objects;
%   (3) apply_op, against the explicit matrix products;
%   (4) uminus, minus, plus, times, mtimes (values compared with the
%       explicit matrix expressions), including sums and products with
%       'tensopvar' operators (4b, 4c);
%   (5) horzcat, vertcat, blkdiag, repmat, repelem, subsref, subsasgn;
%   (6) clean_op, is_id_op;
%   (7) otimes, against the explicit products of functions, including
%       products with 'tensopvar' operators (7b);
%   (8) cost of the constructor, isvalid, clean_op, plus and
%       operator*matrix for a 'dpvar', as a function of the number of
%       decision variables.
% The expected values are computed from the matrices M defining the
% blocks, not from the class itself.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Part (8) asserts a time ratio (q=1e6 within 20x of q=1e3, plus 50 ms)
%   for the constructor, isvalid, clean_op and operator*matrix only; the
%   time of plus is reported, not asserted;
% - Variables of a Kronecker output are named s_1, t_1 by apply_op; the
%   tests strip '_1' before comparing.

tol = 1e-10;
pvar s t
dpvar a b
dom_s = [0,1];      dom_t = [-1,2];
vars_s = struct('in',{cell(1,0)},'out',{{'s'}});        dom_S = struct('in',zeros(0,2),'out',dom_s);
vars_st = struct('in',{cell(1,0)},'out',{{'s','t'}});   dom_ST = struct('in',zeros(0,2),'out',[dom_s;dom_t]);
ncases = 0;

%% (1) Construction
% (1.1) 'double' blocks: no variables
C = tensopmat({[1,2;3,4]});
check_obj(C,{2,2,cell(1,0),zeros(0,2),zeros(1,0),zeros(1,0)},'(1.1)');
C = tensopmat([1,2,3]);         % single matrix
check_obj(C,{1,3,cell(1,0),zeros(0,2),zeros(1,0),zeros(1,0)},'(1.1b)');
% (1.2) 'polynomial' block in one variable
C = tensopmat({[1+s;2*s^2]},vars_s,dom_S);
check_obj(C,{2,1,{'s'},dom_s,1,zeros(1,0)},'(1.2)');
C2 = tensopmat([1+s;2*s^2],vars_s,dom_S);
if ~isequal(C2.vars,C.vars) || ~isequal(C2.dom,C.dom) || ~isequal(size(C2.ops),[1,1])
    error("Test_tensopmat_native_blocks failed (1.2b): single polynomial matrix differs from cell.")
end
% (1.3) variables of different entries: the output space has the union, sorted
C = tensopmat({[s;t^2]},vars_st,dom_ST);
check_obj(C,{2,1,{'s','t'},[dom_s;dom_t],[1,1],zeros(1,0)},'(1.3)');
% (1.4) 'dpvar': decision variables are no spatial variables
C = tensopmat({[a*s+b;3*a*t]},vars_st,dom_ST);
check_obj(C,{2,1,{'s','t'},[dom_s;dom_t],[1,1],zeros(1,0)},'(1.4)');
C = tensopmat({a+2*b},vars_st,dom_ST);
check_obj(C,{1,1,cell(1,0),zeros(0,2),zeros(1,0),zeros(1,0)},'(1.4b)');
% (1.5) only the variables that appear are kept
C = tensopmat({s},vars_st,dom_ST);
check_obj(C,{1,1,{'s'},dom_s,1,zeros(1,0)},'(1.5)');
C = tensopmat({polynomial(3)},vars_st,dom_ST);
check_obj(C,{1,1,cell(1,0),zeros(0,2),zeros(1,0),zeros(1,0)},'(1.5b)');
C = tensopmat({(1+s)-s},vars_st,dom_ST);        % the variable cancels
check_obj(C,{1,1,cell(1,0),zeros(0,2),zeros(1,0),zeros(1,0)},'(1.5c)');
% (1.6) a grid of blocks, each row with its own variables, and different types
ops = {[1+s;s], [2,3+s;4,5]; a*t, [b,a]};
C = tensopmat(ops,vars_st,dom_ST);
if ~isequal(C.dims.out,[2;1]) || ~isequal(C.dims.in,[1;2]) || ~isequal(C.vars.out,{'s','t'}) || ...
        ~isequal(C.depmat.out,[1,0;0,1]) || ~isequal(C.depmat.in,zeros(2,0)) || ~isequal(C.dom.out,[dom_s;dom_t])
    error("Test_tensopmat_native_blocks failed (1.6): grid of blocks.")
end
ncases = ncases+10;
% (1.7) with other blocks: the other blocks set the output space of the row, the variables of the
% blocks M must be among them, and the domains are taken from these blocks
Top = mat2tensopvar([1;2],vars_s,dom_S);
C = tensopmat({Top,[5;6]});                                   % constant in a row that maps to L2[s]
if ~isequal(C.vars.out,{'s'}) || ~isequal(C.depmat.out,1) || ~isequal(C.dom.out,dom_s) || ~isvalid(C)
    error("Test_tensopmat_native_blocks failed (1.7): constant block in a row with other blocks.")
end
C = tensopmat({Top,[5;6]+s});                               % domain of s from the first block
if ~isequal(C.vars.out,{'s'}) || ~isequal(C.depmat.out,1) || ~isequal(C.dom.out,dom_s) || ~isvalid(C)
    error("Test_tensopmat_native_blocks failed (1.7b): domain from the other block.")
end
C = tensopmat({Top;3+s});                                   % domain of s from another row
if ~isequal(C.vars.out,{'s'}) || ~isequal(C.depmat.out,[1;1]) || ~isequal(C.dom.out,dom_s)
    error("Test_tensopmat_native_blocks failed (1.7c): domain from another row.")
end
expect_error(@() tensopmat({Top,[5;6]+t},vars_st,dom_ST),'output variables of the other blocks');
ncases = ncases+4;
% (1.8) the domains from the object itself, as set before calling initialize
C = tensopmat();        C.ops = {[1+s;s*t]};
C.vars = vars_st;       C.dom = dom_ST;
C = initialize(C);
check_obj(C,{2,1,{'s','t'},[dom_s;dom_t],[1,1],zeros(1,0)},'(1.8)');
ncases = ncases+1;
% (1.9) errors
expect_error(@() tensopmat({1+s}),'cannot be determined');
expect_error(@() tensopmat({1+s},vars_t_only(),dom_t_only(dom_t)),'cannot be determined');
expect_error(@() tensopmat({[1+s,2];[3,4,5]},vars_s,dom_S),'dimensions');
expect_error(@() tensopmat({1+s},vars_s),'together');
expect_error(@() tensopmat(Top,vars_s,dom_S),'only be specified');
expect_error(@() tensopmat({'a'}),'must be of type');
% a column with an operator acting on a function space
Top_in = mat2tensopvar(1,struct('in',{{'s'}},'out',{{'s'}}),struct('in',dom_s,'out',dom_s));
expect_error(@() tensopmat({1+s;Top_in},vars_s,dom_S),'act on R');
ncases = ncases+7;

%% (2) isvalid
Cv = tensopmat({[1+s;a*s],[2;3]},vars_s,dom_S);
if ~isvalid(Cv)
    error("Test_tensopmat_native_blocks failed (2.1): valid object not valid.")
end
% a native block in a column acting on an input variable
Cin = tensopmat({Top_in,[1+s]},vars_s,dom_S);           % column 1 acts on L2[s]; column 2 on R
if ~isvalid(Cin) || ~isequal(Cin.depmat.in,[1;0])
    error("Test_tensopmat_native_blocks failed (2.2): object with column acting on L2 and column acting on R.")
end
Cbad = Cin;     Cbad.depmat.in(2,1) = 1;
[lv,msg] = isvalid(Cbad);
if lv || ~contains(msg,'act on R')
    error("Test_tensopmat_native_blocks failed (2.3): depmat.in>0 for a column of matrices not detected: %s",msg)
end
Cbad = Cv;      Cbad.depmat.out(1,1) = 0;               % row 1 has the block with s
[lv,msg] = isvalid(Cbad);
if lv || ~contains(msg,'output variables')
    error("Test_tensopmat_native_blocks failed (2.4): variable not in the output space not detected: %s",msg)
end
Cbad = Cv;      Cbad.ops{1,1} = [1+s;2;3];              % wrong dimension
[lv,msg] = isvalid(Cbad);
if lv
    error("Test_tensopmat_native_blocks failed (2.5): wrong dimension not detected.")
end
ncases = ncases+5;

%% (3) apply_op
P1 = [1+s;s*t];     D1 = [a*s;b+t];    c1 = 3;     D2 = [2,a;b,4];
Cop = tensopmat({P1,D1;c1,5},vars_st,dom_ST);
x1 = 2;     x2 = 7;
Y = apply_op(Cop,{x1,x2});
Yexp = [P1*x1+D1*x2; c1*x1+5*x2];
assert_equal(Y,Yexp,tol,'(3.1)');
% input in the 'dpvar' and 'polynomial'
xd = a+3*b;
Cop2 = tensopmat({P1,[1,2;3,4]},vars_st,dom_ST);
Y = apply_op(Cop2,{xd,[a;b]});
Yexp = P1*xd+[1,2;3,4]*[a;b];
assert_equal(Y,Yexp,tol,'(3.2)');
% the input of a column may also be given as a 1x1 cell, as for the factors of an operator
Y = apply_op(Cop2,{{xd},{[a;b]}});
assert_equal(Y,Yexp,tol,'(3.3)');
ncases = ncases+3;

%% (4) Algebra
Cs = tensopmat({[1+s;s^2]},vars_s,dom_S);
Cd = tensopmat({[a*s;2+b]},vars_s,dom_S);
Cm = uminus(Cs);
assert_equal(Cm.ops{1},-[1+s;s^2],tol,'(4.1)');
Cm = minus(Cs,Cd);
assert_equal(Cm.ops{1},[1+s;s^2]-[a*s;2+b],tol,'(4.2)');
if ~isvalid(Cm) || ~isequal(Cm.vars.out,{'s'})
    error("Test_tensopmat_native_blocks failed (4.2): difference not valid.")
end
% plus: same spaces, different variables (p(s)+2), and different types
Csum = Cs+Cd;
assert_equal(Csum.ops{1},[1+s;s^2]+[a*s;2+b],tol,'(4.3)');
C2 = tensopmat({[2;3]});
Csum = Cs+C2;
assert_equal(Csum.ops{1},[1+s;s^2]+[2;3],tol,'(4.4)');
if ~isequal(Csum.vars.out,{'s'}) || ~isequal(Csum.depmat.out,1) || ~isequal(Csum.dom.out,dom_s) || ~isvalid(Csum)
    error("Test_tensopmat_native_blocks failed (4.4): sum of function and constant.")
end
Csum = C2+Cs;
assert_equal(Csum.ops{1},[1+s;s^2]+[2;3],tol,'(4.4b)');
% two functions of different variables
Ct = tensopmat({[t;t^2]},vars_st,dom_ST);
Csum = Cs+Ct;
assert_equal(Csum.ops{1},[1+s;s^2]+[t;t^2],tol,'(4.5)');
if ~isequal(Csum.vars.out,{'s','t'}) || ~isequal(Csum.depmat.out,[1,1]) || ~isvalid(Csum)
    error("Test_tensopmat_native_blocks failed (4.5): sum of functions in different variables.")
end
% times: scalar and elementwise
Cm = 2*Cs;
assert_equal(Cm.ops{1},2*[1+s;s^2],tol,'(4.6)');
Cm = s*Cs;
assert_equal(Cm.ops{1},s*[1+s;s^2],tol,'(4.7)');
Cm = [2;3].*Cd;
assert_equal(Cm.ops{1},[2;3].*[a*s;2+b],tol,'(4.8)');
expect_error(@() t*Cs,'cannot be determined');                       % domain of t unknown
Cst = vertcat_rows(Cs,Ct);                                          % rows with s, and with t
Cm = t*Cst;
assert_equal(Cm.ops{1},t*[1+s;s^2],tol,'(4.9)');
assert_equal(Cm.ops{2},t*[t;t^2],tol,'(4.9b)');
if ~isequal(Cm.vars.out,{'s','t'}) || ~isequal(Cm.depmat.out,[1,1;0,1]) || ~isvalid(Cm) || ~isequal(Cm.dom.out,[dom_s;dom_t])
    error("Test_tensopmat_native_blocks failed (4.9): variables after multiplication with t.")
end
Cm = t.*Cst;                                                        % scalar, elementwise
assert_equal(Cm.ops{1},t.*[1+s;s^2],tol,'(4.9c)');
if ~isequal(Cm.depmat.out,[1,1;0,1]) || ~isvalid(Cm)
    error("Test_tensopmat_native_blocks failed (4.9c): variables after elementwise multiplication with t.")
end
Cm = [t;1;t^2;3].*Cst;                                              % matrix, elementwise
assert_equal(Cm.ops{1},[t;1].*[1+s;s^2],tol,'(4.9d)');
assert_equal(Cm.ops{2},[t^2;3].*[t;t^2],tol,'(4.9e)');
if ~isequal(Cm.depmat.out,[1,1;0,1]) || ~isvalid(Cm)
    error("Test_tensopmat_native_blocks failed (4.9d): variables after elementwise multiplication with a matrix.")
end
% mtimes: matrix times object, and object times object
Bm = tensopmat({[1+s];[s^2]},vars_s,dom_S);                         % 2 rows of blocks, 1 column
Am = [1,2];
Pm = Am*Bm;
assert_equal(Pm.ops{1},1*(1+s)+2*s^2,tol,'(4.10)');
Ap = tensopmat({[1+s,s]},vars_s,dom_S);                             % R^2 -> L2[s]
Bp = tensopmat({[1;2]});                                            % R^1 -> R^2
Pm = Ap*Bp;
assert_equal(Pm.ops{1},[1+s,s]*[1;2],tol,'(4.11)');
if ~isequal(Pm.vars.out,{'s'}) || ~isvalid(Pm)
    error("Test_tensopmat_native_blocks failed (4.11): product of objects.")
end
Ad = tensopmat({[a*s,b]},vars_s,dom_S);
Pm = Ad*Bp;
assert_equal(Pm.ops{1},[a*s,b]*[1;2],tol,'(4.12)');
Pm = Ap*tensopmat(mat2tensopvar([1;2]));                           % with an operator R^1 -> R^2: the result is an operator
assert_equal(apply_op(Pm,{{3}}),[1+s,s]*[1;2]*3,tol,'(4.12b)');
ncases = ncases+19;

%% (4b) The sum of a block of type 'double'/'polynomial'/'dpvar' with a 'tensopvar' that acts on a constant
% monomial (R): the 'tensopvar' is converted to a 'polynomial', and the matrices added
vars_ss = struct('in',{{'s'}},'out',{{'s'}});       dom_ss = struct('in',dom_s,'out',dom_s);
Top = tensopmat(mat2tensopvar([2+s;3],vars_s,dom_S));       % multiplier from R, to L2[s]
Cm = Cs+Top;
assert_equal(Cm.ops{1},[1+s;s^2]+[2+s;3],tol,'(4.13)');
if ~isequal(Cm.vars.out,{'s'}) || ~isequal(Cm.depmat.out,1) || ~isvalid(Cm) || ~isa(Cm.ops{1},'polynomial')
    error("Test_tensopmat_native_blocks failed (4.13): sum of a matrix and an operator.")
end
Cm = Top+Cs;
assert_equal(Cm.ops{1},[2+s;3]+[1+s;s^2],tol,'(4.14)');
Cm = Cs-Top;
assert_equal(Cm.ops{1},[1+s;s^2]-[2+s;3],tol,'(4.15)');
Cm = Cd+Top;                                                % 'dpvar'
assert_equal(Cm.ops{1},[a*s;2+b]+[2+s;3],tol,'(4.16)');
if ~isa(Cm.ops{1},'dpvar')
    error("Test_tensopmat_native_blocks failed (4.16): the sum is not of type 'dpvar'.")
end
Cm = tensopmat({[1;2]})+tensopmat(mat2tensopvar([2;3]));    % a constant with a constant, no variables
assert_equal(Cm.ops{1},[3;5],tol,'(4.17)');
% in other variables than the matrix: the output space has the union
Ttt = tensopmat(mat2tensopvar([t;1],vars_st,dom_ST));
Cm = Cs+Ttt;
assert_equal(Cm.ops{1},[1+s;s^2]+[t;1],tol,'(4.18)');
if ~isvalid(Cm) || ~isequal(Cm.vars.out,{'s','t'})
    error("Test_tensopmat_native_blocks failed (4.18): output space of the sum.")
end
% an operator with several factors and terms, acting on R: compared with the values of the operator itself
F1 = mat2tensopvar([1+s,2;2,s],vars_s,dom_S);     F2 = mat2tensopvar([s,3;1,s],vars_s,dom_S);
G1 = mat2tensopvar([3,s^2;1,2],vars_s,dom_S);     G2 = mat2tensopvar([2,s;s,1],vars_s,dom_S);
T1 = otimes(F1,F2,[1,2]);                       T2 = otimes(G1,G2,[1,2]);
N2 = tensopmat({[1,2,3,4;5,6,7,8]+s*[0,1,0,1;1,0,1,0]},vars_s,dom_S);       % 2x4, as the operators
x1 = [2;1];     x2 = [1;-3];       xx = kron(x1,x2);
for kk=1:2
    if kk==1,   Tk = tensopmat(T1);     else,   Tk = tensopmat(T1+T2);      end     % one term, two terms
    Cm = N2+Tk;
    if size(Cm,2)~=4 || ~isvalid(Cm) || ~isa(Cm.ops{1},'polynomial')
        error("Test_tensopmat_native_blocks failed (4.19.%d): sum with a multi-factor operator.",kk)
    end
    Yexp = apply_op(N2,{xx})+apply_op(Tk,{{x1,x2}});    % the operator applied to the factors of the monomial
    assert_equal(apply_op(Cm,{xx}),Yexp,tol,sprintf('(4.19.%d)',kk));
end
% a Kronecker product of the rows, in different variables
Kt = otimes(mat2tensopvar([1+s;2],vars_s,dom_S),mat2tensopvar([t;1],vars_st,dom_ST),[2,2]);     % 4x1, in s and t
Nk = tensopmat({[1;2;3;4]*(1+0*s)+s*[0;1;0;1]},vars_s,dom_S);
Cm = Nk+tensopmat(Kt);
Yk = apply_op(tensopmat(Kt),{{1,1}});       Yk.varname = strrep(Yk.varname,'_1','');
assert_equal(Cm.ops{1},[1;2;3;4]+s*[0;1;0;1]+Yk,tol,'(4.20)');
% an operator acting on a function space, a functional, or factors sharing a variable in a Kronecker
% product of the rows, cannot be converted
Tin = tensopmat(mat2tensopvar(1,vars_ss,dom_ss));
Kint = tensopmat(intop({1},{[0,2]},{'s'},[0,1],true));
Cn1 = tensopmat({1+s},vars_s,dom_S);
expect_error(@() Cn1+Tin,'acts on R');
expect_error(@() Tin+Cn1,'acts on R');
expect_error(@() Cn1+Kint,'acts on R');
Tsh = tensopmat(otimes(mat2tensopvar([1+s;2],vars_s,dom_S),mat2tensopvar([s;1],vars_s,dom_S),[2,2]));    % factors share s: no polynomial
expect_error(@() Nk+Tsh,'acts on R');
% as the sum in the product of a matrix with an object (mtimes)
Bmix = [tensopmat({1+s},vars_s,dom_S);tensopmat(mat2tensopvar(s^2,vars_s,dom_S))];
Pm = [1,2]*Bmix;
assert_equal(Pm.ops{1},1*(1+s)+2*s^2,tol,'(4.21)');
Bmix = [tensopmat(mat2tensopvar(s^2,vars_s,dom_S));tensopmat({1+s},vars_s,dom_S)];     % the other order
Pm = [1,2]*Bmix;
assert_equal(Pm.ops{1},1*s^2+2*(1+s),tol,'(4.22)');
ncases = ncases+15;

%% (4c) The product of a block of type 'double'/'polynomial'/'dpvar' with an operator (mtimes)
% (b) operator times matrix: the operator is applied to the columns of the matrix, giving a matrix
q = 2+s;        Aq = tensopmat(mat2tensopvar(q,vars_ss,dom_ss));                % multiplier on L2[s]
Bn = tensopmat({[1+s^2,s]},vars_s,dom_S);                                       % R^2 -> L2[s]
Pm = Aq*Bn;
assert_equal(Pm.ops{1},q*[1+s^2,s],tol,'(4.23)');
if ~isvalid(Pm) || ~isa(Pm.ops{1},'polynomial') || ~isequal(Pm.vars.out,{'s'}) || ~isequal(Pm.dims.in,2)
    error("Test_tensopmat_native_blocks failed (4.23): product of an operator and a matrix.")
end
% a functional: the integral of the functions, as a constant matrix
Fint = tensopmat(mat2tensopvar(3,struct('in',{{'s'}},'out',{cell(1,0)}),struct('in',dom_s,'out',zeros(0,2))));  % x -> 3 int_0^1 x
Pm = Fint*Bn;
assert_equal(Pm.ops{1},[3*(1+1/3),3/2],tol,'(4.24)');
if ~isvalid(Pm) || ~isempty(Pm.vars.out)
    error("Test_tensopmat_native_blocks failed (4.24): functional applied to a matrix.")
end
Kint = tensopmat(intop({1},{[0,2]},{'s'},[0,1],true));                          % x -> int_0^1 x
Pm = Kint*Bn;
assert_equal(Pm.ops{1},[4/3,1/2],tol,'(4.25)');
% vector-valued, an operator with a matrix of multipliers
qq = [2+s,0;1,3-s];
Aqq = tensopmat(mat2tensopvar(qq,vars_ss,dom_ss));
Bn2 = tensopmat({[1+s,s^2;2,s]},vars_s,dom_S);
Pm = Aqq*Bn2;
assert_equal(Pm.ops{1},qq*[1+s,s^2;2,s],tol,'(4.26)');
% a single factor with a Kronecker product of the (single) factor: no renaming of the variables
Aqk = tensopmat(sopvar2tensopvar(mat2sopvar(q,vars_ss,dom_ss),[2,2]));
Pm = Aqk*Bn;
assert_equal(Pm.ops{1},q*[1+s^2,s],tol,'(4.27)');
% a 'dpvar': the operator is applied to the monomials, and the coefficients combined
Bd = tensopmat({[a*s+b,2*a]},vars_s,dom_S);
Pm = Aq*Bd;
assert_equal(Pm.ops{1},q*[a*s+b,2*a],tol,'(4.28)');
if ~isa(Pm.ops{1},'dpvar') || ~isvalid(Pm)
    error("Test_tensopmat_native_blocks failed (4.28): the product with a block of type 'dpvar' is no 'dpvar'.")
end
Pm = Fint*Bd;
assert_equal(Pm.ops{1},[3*(a/2+b),6*a],tol,'(4.29)');
Pm = Kint*Bd;
assert_equal(Pm.ops{1},[a/2+b,2*a],tol,'(4.30)');
% a 'dpvar' depending on two spatial variables, of which the operator acts on one
Bd2 = tensopmat({[a*s*t+b*t^2;s+a]},vars_st,dom_ST);
Aq2 = tensopmat(mat2tensopvar(diag([q,q+1]),struct('in',{{'s','t'}},'out',{{'s','t'}}),struct('in',[dom_s;dom_t],'out',[dom_s;dom_t])));
Pm = Aq2*Bd2;
assert_equal(Pm.ops{1},[q*(a*s*t+b*t^2);(q+1)*(s+a)],tol,'(4.31)');
% (a) matrix times operator that maps to R: the matrix is converted to a multiplier operator
Mp = tensopmat({[1+s;s^2]},vars_s,dom_S);
Pm = Mp*Fint;
x = 2+s^2;      Ix = 3*(2+1/3);                                                 % 3*int_0^1 x
assert_equal(apply_op(Pm,{x}),[1+s;s^2]*Ix,tol,'(4.32)');
if ~isvalid(Pm) || ~isequal(Pm.dims.out,2) || ~isequal(Pm.vars.out,{'s'}) || ~isequal(Pm.vars.in,{'s'})
    error("Test_tensopmat_native_blocks failed (4.32): matrix times functional.")
end
Pm = tensopmat({2})*Fint;                                                       % a constant
assert_equal(apply_op(Pm,{x}),2*Ix,tol,'(4.33)');
% a matrix of products: a block matrix with a native and an operator column, times native and constant blocks;
% sum of the product of the matrix (a) and the operator (b), added after converting to a polynomial
Aml = [tensopmat({1+s},vars_s,dom_S),Aq];               % [p(s), q(s)]: R x L2[s] -> L2[s]
Bml = [tensopmat(mat2tensopvar(2));tensopmat({3+s^2},vars_s,dom_S)];   % R -> R x L2[s]
Pm = Aml*Bml;
assert_equal(Pm.ops{1},(1+s)*2+q*(3+s^2),tol,'(4.34)');
Bml = [tensopmat({3+s^2},vars_s,dom_S);tensopmat(mat2tensopvar(2))];
Aml = [Aq,tensopmat({1+s},vars_s,dom_S)];
Pm = Aml*Bml;                                                                   % the other order
assert_equal(Pm.ops{1},q*(3+s^2)+(1+s)*2,tol,'(4.35)');
% an integral operator, a product of the multiplier and the functional: applied to matrices
Aint = Mp*Fint;                                                                 % x -> [1+s;s^2]*3*int x, mapping L2[s] to L2^2[s]
Pm = Aint*Bn;
assert_equal(Pm.ops{1},[1+s;s^2]*3*[4/3,1/2],tol,'(4.36)');
Pm = Aint*Bd;
assert_equal(Pm.ops{1},[1+s;s^2]*3*[a/2+b,2*a],tol,'(4.37)');
% rows of matrices in different variables: each converted to the multiplier to the output space of its own row
Mrows = [tensopmat({1+s},vars_s,dom_S);tensopmat({t^2+1},vars_st,dom_ST)];
Pm = Mrows*Fint;
Yr = apply_op(Pm,{x});
assert_equal(Yr,[(1+s)*Ix;(t^2+1)*Ix],tol,'(4.38)');
if ~isvalid(Pm) || ~isequal(Pm.depmat.out,[1,0;0,1])
    error("Test_tensopmat_native_blocks failed (4.38): output spaces of the rows.")
end
% not supported: a 'dpvar' times an operator (the operators cannot depend on decision variables), and a
% matrix times a functional of type 'intop'
expect_error(@() tensopmat({a*s},vars_s,dom_S)*Fint,'decision variables');
expect_error(@() Mp*Kint,'functional');
ncases = ncases+19;

%% (5) Concatenation, repetition, indexing, assignment
Cs2 = tensopmat({[3+s;t]},vars_st,dom_ST);
Ch = [Cs,Cs2];
if ~isequal(size(Ch),[2,2]) || ~isvalid(Ch) || ~isequal(Ch.vars.out,{'s','t'}) || ~isequal(Ch.dom.out,[dom_s;dom_t])
    error("Test_tensopmat_native_blocks failed (5.1): horzcat.")
end
assert_equal(Ch.ops{1,2},[3+s;t],tol,'(5.1)');
Ch = [Cs,tensopmat({[7;8]})];
if ~isvalid(Ch) || ~isequal(Ch.depmat.out,1) || ~isequal(Ch.vars.out,{'s'})
    error("Test_tensopmat_native_blocks failed (5.2): horzcat of function and constant.")
end
Ch = [tensopmat({[7;8]}),Cs];
if ~isvalid(Ch) || ~isequal(Ch.depmat.out,1) || ~isequal(Ch.vars.out,{'s'})
    error("Test_tensopmat_native_blocks failed (5.3): horzcat of constant and function.")
end
Ch = [tensopmat(mat2tensopvar([1;2],vars_s,dom_S)),tensopmat({[7;8]})];
if ~isvalid(Ch) || ~isequal(Ch.depmat.out,1)
    error("Test_tensopmat_native_blocks failed (5.4): horzcat of operator and constant.")
end
Cv = [Cs;Cd];
if ~isequal(size(Cv),[4,1]) || ~isvalid(Cv) || ~isequal(Cv.dims.out,[2;2])
    error("Test_tensopmat_native_blocks failed (5.5): vertcat.")
end
Cb = blkdiag(Cs,Ct);
if ~isequal(size(Cb),[4,2]) || ~isvalid(Cb) || ~isequal(Cb.depmat.out,[1,0;0,1]) || ~isempty(Cb.ops{1,2}) || ~isempty(Cb.ops{2,1})
    error("Test_tensopmat_native_blocks failed (5.6): blkdiag.")
end
Cr = repmat(Cs,2,3);
if ~isequal(size(Cr),[4,3]) || ~isvalid(Cr)
    error("Test_tensopmat_native_blocks failed (5.7): repmat.")
end
assert_equal(Cr.ops{2,3},[1+s;s^2],tol,'(5.7)');
Cr = repelem(Cs,2,3);
if ~isequal(size(Cr),[4,3]) || ~isvalid(Cr)
    error("Test_tensopmat_native_blocks failed (5.8): repelem.")
end
Csub = Cs(2,1);
if ~isequal(size(Csub),[1,1]) || ~isvalid(Csub)
    error("Test_tensopmat_native_blocks failed (5.9): indexing.")
end
assert_equal(Csub.ops{1},s^2,tol,'(5.9)');
Ca = Ch;        Ca{1,1} = [4*s;5];
if ~isvalid(Ca)
    error("Test_tensopmat_native_blocks failed (5.10): assignment.")
end
assert_equal(Ca.ops{1,1},[4*s;5],tol,'(5.10)');
Ca = Cs;        Ca{1,1} = [a;2*b*s];
assert_equal(Ca.ops{1,1},[a;2*b*s],tol,'(5.11)');
expect_error(@() assign_bad(Cs),'must be specified');
ncases = ncases+11;

%% (6) clean_op, is_id_op
Cn = tensopmat({[1+1e-14*s;s]},vars_s,dom_S);
[Cc,is_zero] = clean_op(Cn,1e-10);
assert_equal(Cc.ops{1},[1;s],tol,'(6.1)');
if is_zero || ~isvalid(Cc)
    error("Test_tensopmat_native_blocks failed (6.1): clean_op.")
end
[Cc,is_zero] = clean_op(tensopmat({[1e-14*s;0*s+1e-15]},vars_s,dom_S),1e-10);
if ~is_zero || ~isvalid(Cc)
    error("Test_tensopmat_native_blocks failed (6.2): zero block not detected.")
end
Cc = clean_op(tensopmat({[a*(1+1e-14*s);b]},vars_s,dom_S),1e-10);
assert_equal(Cc.ops{1},[a;b],tol,'(6.3)');
if ~is_id_op(tensopmat({eye(2)})) || ~is_id_op(tensopmat({polynomial(eye(2)+1e-14*s)},vars_s,dom_S))
    error("Test_tensopmat_native_blocks failed (6.4): identity not detected.")
end
if is_id_op(tensopmat({2*eye(2)})) || is_id_op(tensopmat({[1,0;0,1]+s*[0,1;0,0]},vars_s,dom_S)) || is_id_op(tensopmat({[1,1;0,1]}))
    error("Test_tensopmat_native_blocks failed (6.5): non-identity found to be identity.")
end
% an identity block matrix, and a block that is the identity matrix in a row that maps to L2[s]
Cid = tensopmat({eye(2),zeros(2,1);zeros(1,2),polynomial(1)});
if ~is_id_op(Cid)
    error("Test_tensopmat_native_blocks failed (6.6): identity block matrix.")
end
Cid = tensopmat({eye(2),mat2tensopvar([0;0],vars_s,dom_S);[],1});
if is_id_op(Cid) || ~isvalid(Cid)
    error("Test_tensopmat_native_blocks failed (6.7): identity matrix in a row mapping to L2[s].")
end
ncases = ncases+8;

%% (7) otimes
% (A*x).*(B*y) = (A (x) B)*kron(x,y), for rows multiplied pointwise
A = tensopmat({[1+s,2;s,3]},vars_s,dom_S);      B = tensopmat({[2-s,t;1,a]},vars_st,dom_ST);
x = [1;-2];     y = [3;5];
Cab = otimes(A,B,[1,2]);
Yexp = ([1+s,2;s,3]*x).*([2-s,t;1,a]*y);
assert_equal(apply_op(Cab,{kron(x,y)}),Yexp,tol,'(7.1)');
if ~isvalid(Cab) || ~isequal(Cab.dims.in,4) || ~isequal(Cab.dims.out,2) || ~isequal(Cab.vars.out,{'s','t'}) || ~isequal(Cab.depmat.out,[1,1])
    error("Test_tensopmat_native_blocks failed (7.1): product object.")
end
% pointwise rows and columns: the elementwise product A.*B
Cab = otimes(A,B,[1,1]);
Yexp = ([1+s,2;s,3].*[2-s,t;1,a])*x;
assert_equal(apply_op(Cab,{x}),Yexp,tol,'(7.2)');
% Kronecker rows, in different variables
A1 = tensopmat({[1+s;2]},vars_s,dom_S);         B1 = tensopmat({[t;a*t]},vars_st,dom_ST);
Cab = otimes(A1,B1,[2,2]);
Yexp = [(1+s)*t; (1+s)*(a*t); 2*t; 2*(a*t)];            % kron([1+s;2],[t;a*t]), 'dpvar' has no kron
assert_equal(apply_op(Cab,{1}),Yexp,tol,'(7.3)');
if ~isequal(Cab.dims.out,4) || ~isvalid(Cab)
    error("Test_tensopmat_native_blocks failed (7.3): Kronecker product.")
end
expect_error(@() otimes(A1,A1,[2,2]),'share a spatial variable');
% constants
Cab = otimes(tensopmat({[2,3]}),tensopmat({[1+s,s]},vars_s,dom_S),[1,2]);
Yexp = ([2,3]*[1;2]).*([1+s,s]*[3;4]);
assert_equal(apply_op(Cab,{kron([1;2],[3;4])}),Yexp,tol,'(7.4)');
ncases = ncases+7;

%% (7b) otimes of a block of type 'double'/'polynomial' with a 'tensopvar': the block is converted to the
% multiplier operator from R to the output space of its row, and the regular product is taken. The
% result has two factors, acting on R and on the function space: (A*x).*(B*y) = (A (x) B)*{x,y}.
dom_ss = struct('in',dom_s,'out',dom_s);        vars_ss = struct('in',{{'s'}},'out',{{'s'}});
q = 2+s;        qq = [2+s,0;1,3-s];
Bq = tensopmat(mat2tensopvar(q,vars_ss,dom_ss));        % multiplier on L2[s]
y = 3+s^2;
% polynomial block times operator, and operator times polynomial block
Ap = tensopmat({1+s^2},vars_s,dom_S);
Cab = otimes(Ap,Bq,[1,2]);
assert_equal(apply_op(Cab,{{2,y}}),(1+s^2)*2*(q*y),tol,'(7.5)');
if ~isvalid(Cab) || ~isequal(Cab.dims.in,1) || ~isequal(Cab.depmat.out,1) || ~isequal(Cab.vars.out,{'s'}) || ~isequal(size(Cab.ops{1}.ops),[1,2])
    error("Test_tensopmat_native_blocks failed (7.5): product object.")
end
Cba = otimes(Bq,Ap,[1,2]);
assert_equal(apply_op(Cba,{{y,2}}),(q*y)*(1+s^2)*2,tol,'(7.6)');
% vector-valued, with a Kronecker product of the columns
Av = tensopmat({[1+s,2;s,3]},vars_s,dom_S);
Bv = tensopmat(mat2tensopvar(qq,vars_ss,dom_ss));
yv = [3+s^2;1-s];       xv = [1;-2];
Cab = otimes(Av,Bv,[1,2]);
assert_equal(apply_op(Cab,{{xv,yv}}),([1+s,2;s,3]*xv).*(qq*yv),tol,'(7.7)');
% a constant (no variables): the row of the other operator sets the output space
Cc = tensopmat({[2;3]});
Cab = otimes(Cc,Bv,[1,2]);
assert_equal(apply_op(Cab,{{1,yv}}),[2;3].*(qq*yv),tol,'(7.8)');
% a constant in a row that also contains an operator
Crow = [tensopmat(mat2tensopvar([1;1],vars_s,dom_S)),tensopmat({[5;6]})];      % row maps to L2[s]
Cab = otimes(Crow,Bv,[1,2]);
assert_equal(apply_op(Cab,{{1,yv},{1,yv}}),([1;1]+[5;6]).*(qq*yv),tol,'(7.9)');
% Kronecker product of the rows: the block of type 'polynomial' in other variables than the operator
At = tensopmat({[t;1]},vars_st,dom_ST);
Cab = otimes(At,Bq,[2,2]);
Yexp = [t*(q*y); 1*(q*y)];
Y = apply_op(Cab,{{1,y}});
Y.varname = strrep(Y.varname,'_1','');          % apply_op names the variables of a Kronecker output s_1, t_1
assert_equal(Y,Yexp,tol,'(7.10)');
if ~isequal(Cab.depmat.out,[1,1])
    error("Test_tensopmat_native_blocks failed (7.10): output space of the product.")
end
% the block of type 'polynomial' as second factor, in other variables than the first: its row, not that of the first
Cba = otimes(Bq,At,[2,2]);
Y = apply_op(Cba,{{y,1}});
Y.varname = strrep(Y.varname,'_1','');
assert_equal(Y,[(q*y)*t; (q*y)*1],tol,'(7.11)');
% a scalar constant scales a functional ('intop'); a function times a functional is not supported
Kint = tensopmat(intop({1},{[0,2]},{'s'},[0,1],true));          % x -> int_0^1 x(s) ds
Cab = otimes(tensopmat({2}),Kint,[1,2]);
assert_equal(apply_op(Cab,{1+s}),polynomial(3),tol,'(7.12)');   % 2*int_0^1 (1+s) ds
Cba = otimes(Kint,tensopmat({2}),[1,2]);
assert_equal(apply_op(Cba,{1+s}),polynomial(3),tol,'(7.13)');
expect_error(@() otimes(Ap,Kint,[1,2]),'functional');
expect_error(@() otimes(Kint,Ap,[1,2]),'functional');
% a block of type 'dpvar' cannot be converted, as the operators cannot depend on decision variables
expect_error(@() otimes(tensopmat({a*s},vars_s,dom_S),Bq,[1,2]),'decision variables');
expect_error(@() otimes(Bq,tensopmat({a*s},vars_s,dom_S),[1,2]),'decision variables');
ncases = ncases+13;

%% (8) Cost for a number of decision variables q
% The blocks have a row for every decision variable in the coefficients of a 'dpvar', so the cost of
% the operations on the blocks must be linear in the number of nonzero coefficients, not in q
qs = [1e3,1e5,1e6];
times_q = zeros(numel(qs),5);
Amul = tensopmat(mat2tensopvar(diag([2+s,3+s]),vars_ss,dom_ss));      % multiplier on L2^2[s], applied to the 2x2 'dpvar'
for iq=1:numel(qs)
    q = qs(iq);
    nz = 20;                                    % nonzero coefficients, independent of q
    idcs = round(linspace(2,q+1,nz))';          % rows: constant term first, then q decision variables
    Cq = sparse(idcs,repmat([1;2],nz/2,1),1+rand(nz,1),(q+1)*2,4);     % 2x2 matrix of 'dpvar', 2 monomials
    dvarnames = cellstr(compose('d%d',(1:q)'));
    Dq = dpvar(Cq,[0;1],{'s'},dvarnames,[2,2]);
    t0 = tic;       Cop = tensopmat({Dq},vars_s,dom_S);         times_q(iq,1) = toc(t0);
    t0 = tic;       lv = isvalid(Cop);                          times_q(iq,2) = toc(t0);
    t0 = tic;       Cc = clean_op(Cop,1e-12);                   times_q(iq,3) = toc(t0);
    t0 = tic;       Cp = Cop+Cop;                               times_q(iq,4) = toc(t0);
    t0 = tic;       Pq = Amul*Cop;                              times_q(iq,5) = toc(t0);
    if ~lv || ~isequal(Cop.vars.out,{'s'}) || ~isequal(Cc.vars.out,{'s'}) || ~isa(Pq.ops{1},'dpvar') || ~isequal(size(Pq),[2,2])
        error("Test_tensopmat_native_blocks failed (8): object for q=%d.",q)
    end
end
% q=1e6 within a factor 20 of q=1e3 (plus 50 ms); linear in q would be 1000x. plus is only reported
if any(times_q(end,[1:3,5])>20*times_q(1,[1:3,5])+0.05)
    error("Test_tensopmat_native_blocks failed (8): cost grows with the number of decision variables: %s s for q=%s.",mat2str(times_q,3),mat2str(qs))
end
ncases = ncases+1;
fprintf('Test_tensopmat_native_blocks: times [constructor, isvalid, clean_op, plus, operator*matrix] for q=%s: %s s\n',mat2str(qs),mat2str(times_q,2))

fprintf('Test_tensopmat_native_blocks passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function check_obj(C,exp,label)
% exp = {m,n,vars.out,dom.out,depmat.out,depmat.in} for a single-block object
[m,n,vout,dout,dpout,dpin] = exp{:};
if ~isequal(C.dims.out,m) || ~isequal(C.dims.in,n) || ~isequal(C.vars.out,vout) || ~isequal(C.vars.in,cell(1,0)) ...
        || ~isequal(C.dom.out,dout) || ~isequal(C.depmat.out,dpout) || ~isequal(C.depmat.in,dpin) || ~isequal(C.dom.in,zeros(0,2))
    error("Test_tensopmat_native_blocks failed %s: dims %s,%s, vars %s, dom %s, depmat.out %s, depmat.in %s.",label, ...
        mat2str(C.dims.out),mat2str(C.dims.in),strjoin(C.vars.out,','),mat2str(C.dom.out),mat2str(C.depmat.out),mat2str(C.depmat.in))
end
[lv,msg] = isvalid(C);
if ~lv
    error("Test_tensopmat_native_blocks failed %s: isvalid: %s",label,msg)
end
end

function assert_equal(Y,Yref,tol,label)
% Y and Yref of type double, polynomial or dpvar of the same size, equal up to tol
if ~isequal(size(Y),size(Yref))
    error("Test_tensopmat_native_blocks failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
D = Y-Yref;
if isa(D,'double'),         c = D;
elseif isa(D,'dpvar'),      c = D.C;
else,                       c = D.coefficient;      end
err = full(max([0;abs(nonzeros(c))]));
if err>tol
    error("Test_tensopmat_native_blocks failed %s: residual %.3e.",label,err)
end
% not vacuous: the reference is nonzero
if isa(Yref,'double'),      cr = Yref;
elseif isa(Yref,'dpvar'),   cr = Yref.C;
else,                       cr = Yref.coefficient;  end
if isempty(nonzeros(cr))
    error("Test_tensopmat_native_blocks %s: reference value is zero; test would be vacuous.",label)
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("Test_tensopmat_native_blocks failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("Test_tensopmat_native_blocks failed: expected an error containing '%s', but none was thrown.",str)
end

function v = vars_t_only()
v = struct('in',{cell(1,0)},'out',{{'t'}});
end

function d = dom_t_only(dom_t)
d = struct('in',zeros(0,2),'out',dom_t);
end

function C = vertcat_rows(A,B)
% vertical concatenation of two single-block objects
C = [A;B];
end

function C = assign_bad(C)
C{1,1} = 'a';
end
