% Test_get_tensorvars.m -- verify the function
% get_tensorvars(vars,degs,opts), which returns the names si_1,...,si_di
% (or si_dum1,...,si_dumdi) of the copies of spatial variables si.
%
% The names are checked against the definition, written out by hand: inputs
% of type char and cellstr (row and column), scalar, vector and zero
% degrees, more than 9 copies, the options 'dum', 'full' and 'bare' (Part
% 1);
% inputs of type polynomial, kept in the order of the elements (not the
% sorted order of varname), and the options 'pvar'/'polynomial' and
% 'cellstr'/'cell' to choose the output type (Part 2); the names used by
% the classes: Kop.varnames.in, Kop.varnames.in_dum and Kop.pvars.in of
% 'intop' objects for N=1-3 variables of different degrees, and the
% variables of the output of 'apply_op' for a tensor-product (otype(1)=2)
% output, also with a single copy of a variable, and the names accepted by
% 'apply_op' of 'intop': the copies, or the bare variable for a single copy
% only (Part 3); invalid inputs (Part 4).
%
% DJ (CC), 11/19/2026: Initial coding, as Test_get_varcopies, renamed with
%   the function to Test_get_tensorvars;
% DJ (CC), 11/19/2026: Part 3: names accepted by 'apply_op' of 'intop';
% DJ (CC), 11/19/2026: Option 'bare'; Part 3: names written out, as the
%   classes now take their names from get_tensorvars;

pvar s1 s2 s3

%% Part 1: inputs of type char and cellstr
chk(isequal(get_tensorvars('s1',3),{{'s1_1','s1_2','s1_3'}}),'char input');
chk(isequal(get_tensorvars('s1',1),{{'s1_1'}}),'char input, a single copy is numbered');
chk(isequal(get_tensorvars({'s1','s2'},[2,1]),{{'s1_1','s1_2'},{'s2_1'}}),'cellstr input, vector degrees');
chk(isequal(get_tensorvars({'s1';'s2'},[2;1]),{{'s1_1','s1_2'},{'s2_1'}}),'column cellstr input, column degrees');
chk(isequal(get_tensorvars({'s1','s2'},2),{{'s1_1','s1_2'},{'s2_1','s2_2'}}),'scalar degree for all variables');
chk(isequal(get_tensorvars({'s2','s1'},[1,2]),{{'s2_1'},{'s1_1','s1_2'}}),'variables kept in the given order');
chk(isequal(get_tensorvars({'x','theta'},[1,2]),{{'x_1'},{'theta_1','theta_2'}}),'names of different lengths');
names10 = get_tensorvars('s',11);
chk(isequal(names10{1}([1,10,11]),{'s_1','s_10','s_11'}) && numel(unique(names10{1}))==11,'more than 9 copies');
% A variable of degree 0 has no copies: a 1 x 0 cellstr
out = get_tensorvars({'s1','s2'},[0,2]);
chk(isequal(size(out),[1,2]) && iscell(out{1}) && isequal(size(out{1}),[1,0]) && isequal(out{2},{'s2_1','s2_2'}),'degree 0');
out = get_tensorvars(cell(1,0),[]);
chk(iscell(out) && isequal(size(out),[1,0]),'no variables');
out = get_tensorvars(cell(1,0),2);
chk(iscell(out) && isequal(size(out),[1,0]),'no variables, scalar degree');
% Option 'dum': the dummy variables
chk(isequal(get_tensorvars({'s1','s2'},[2,1],'dum'),{{'s1_dum1','s1_dum2'},{'s2_dum1'}}),'option dum');
chk(isequal(get_tensorvars('s1',1,{'dum'}),{{'s1_dum1'}}),'option dum, a single copy is numbered');
% Option 'full': one cellstr with all copies, by variable and then by copy
chk(isequal(get_tensorvars({'s1','s2'},[2,1],'full'),{'s1_1','s1_2','s2_1'}),'option full');
chk(isequal(get_tensorvars({'s1','s2','s3'},[1,0,2],'full'),{'s1_1','s3_1','s3_2'}),'option full, degree 0');
chk(isequal(get_tensorvars({'s1','s2'},[2,1],{'full','dum'}),{'s1_dum1','s1_dum2','s2_dum1'}),'options full and dum');
chk(isequal(get_tensorvars({'s1','s2'},[2,1],{'dum','full'}),{'s1_dum1','s1_dum2','s2_dum1'}),'options dum and full, other order');
out = get_tensorvars({'s1','s2'},0,'full');
chk(iscell(out) && isequal(size(out),[1,0]),'option full, no copies');
% Option 'bare': the variable itself (si_dum with 'dum') for a single copy
% only, numbered copies otherwise
chk(isequal(get_tensorvars({'s1','s2','s3'},[1,2,0],'bare'),{{'s1'},{'s2_1','s2_2'},cell(1,0)}),'option bare');
chk(isequal(get_tensorvars({'s1','s2'},[1,2],{'bare','dum'}),{{'s1_dum'},{'s2_dum1','s2_dum2'}}),'options bare and dum');
chk(isequal(get_tensorvars({'s1','s2'},[1,2],{'full','bare'}),{'s1','s2_1','s2_2'}),'options bare and full');
chk(isequal(get_tensorvars({'s1','s2'},1,{'bare','full'}),{'s1','s2'}),'option bare, all of degree 1: the variables');
chk(isequal(get_tensorvars('s',11,'bare'),get_tensorvars('s',11)),'option bare, more than one copy');
% Default options: the same as no options
chk(isequal(get_tensorvars({'s1','s2'},[2,1],{}),get_tensorvars({'s1','s2'},[2,1])),'empty options');
chk(isequal(get_tensorvars({'s1','s2'},[2,1],'cellstr'),{{'s1_1','s1_2'},{'s2_1'}}),'cellstr input, option cellstr');
disp('Test_get_tensorvars Part 1 (char and cellstr) passed.')

%% Part 2: inputs of type polynomial, and the output type
out = get_tensorvars([s1;s2],[2,1]);
chk(iscell(out) && isequal(size(out),[1,2]) && all(cellfun(@(a) isa(a,'polynomial'),out)),'pvar input: cell of polynomials');
chk(isequal(size(out{1}),[2,1]) && isequal(size(out{2}),[1,1]),'pvar input: column vectors');
chk(isequal(pvar2varname(out{1}),{'s1_1';'s1_2'}) && isequal(pvar2varname(out{2}),{'s2_1'}),'pvar input: names');
chk(ispvar(out{1}) && ispvar(out{2}),'pvar input: vectors of pvar objects');
out = get_tensorvars([s1,s2],[2,1]);
chk(isequal(pvar2varname(out{1}),{'s1_1';'s1_2'}) && isequal(pvar2varname(out{2}),{'s2_1'}),'row pvar input');
% The order of the elements, not the sorted order of varname
out = get_tensorvars([s2;s1],[1,2]);
chk(isequal(pvar2varname(out{1}),{'s2_1'}) && isequal(pvar2varname(out{2}),{'s1_1';'s1_2'}),'pvar input, unsorted');
out = get_tensorvars(s1,2);
chk(isequal(size(out),[1,1]) && isequal(pvar2varname(out{1}),{'s1_1';'s1_2'}),'single pvar input');
out = get_tensorvars([s1;s2],[2,1],{'dum','full'});
chk(isa(out,'polynomial') && isequal(size(out),[3,1]) && isequal(pvar2varname(out),{'s1_dum1';'s1_dum2';'s2_dum1'}),'pvar input, options dum and full');
% Degree 0: an empty polynomial
out = get_tensorvars([s1;s2],[0,1]);
chk(isa(out{1},'polynomial') && isempty(out{1}) && isequal(pvar2varname(out{2}),{'s2_1'}),'pvar input, degree 0');
out = get_tensorvars([s1;s2],0,'full');
chk(isa(out,'polynomial') && isempty(out),'pvar input, option full, no copies');
% Options 'cellstr'/'cell' and 'pvar'/'polynomial' choose the output type
ref = {{'s1_1','s1_2'},{'s2_1'}};
chk(isequal(get_tensorvars([s1;s2],[2,1],'cellstr'),ref) && isequal(get_tensorvars([s1;s2],[2,1],'cell'),ref),'pvar input, option cellstr');
chk(isequal(get_tensorvars([s1;s2],[2,1],{'cell','full'}),[ref{:}]),'pvar input, options cell and full');
for opt={'pvar','polynomial'}
    out = get_tensorvars({'s1','s2'},[2,1],opt{1});
    chk(iscell(out) && isa(out{1},'polynomial') && isequal(pvar2varname(out{1}),{'s1_1';'s1_2'}) && isequal(pvar2varname(out{2}),{'s2_1'}),['cellstr input, option ',opt{1}]);
end
out = get_tensorvars('s1',2,{'pvar','full','dum'});
chk(isa(out,'polynomial') && isequal(pvar2varname(out),{'s1_dum1';'s1_dum2'}),'char input, options pvar, full and dum');
out = get_tensorvars([s1;s2],[2,1],'pvar');
chk(iscell(out) && isequal(pvar2varname(out{1}),{'s1_1';'s1_2'}) && isequal(pvar2varname(out{2}),{'s2_1'}),'pvar input, option pvar');
out = get_tensorvars([s1;s2],[1,2],{'bare','dum'});
chk(isequal(pvar2varname(out{1}),{'s1_dum'}) && isequal(pvar2varname(out{2}),{'s2_dum1';'s2_dum2'}),'pvar input, options bare and dum');
disp('Test_get_tensorvars Part 2 (polynomial) passed.')

%% Part 3: the names used by the classes
% intop: Kop.varnames.in, Kop.varnames.in_dum, and Kop.pvars.in, for
% N=1-3 variables, each integrated over [0,1] in each copy
cases = {{'s1'},1; {'s1'},3; {'s1','s2'},[1,1]; {'s1','s2'},[2,1]; {'s1','s2'},[1,3]; {'s1','s2','s3'},[2,1,2]};
for k=1:size(cases,1)
    [vars,degs] = cases{k,:};
    intmat = arrayfun(@(d) repmat([0,d+1],1,d),degs,'UniformOutput',false);
    Kop = intop({1},intmat,vars,[0,1]);
    msg = sprintf('intop, degrees [%s]',num2str(degs));
    chk(isequal(Kop.deps.in,degs),[msg,', degrees of the functional']);
    chk(isequal(Kop.varnames.in,get_tensorvars(Kop.vars.in,Kop.deps.in)),[msg,', varnames.in']);
    chk(isequal(Kop.varnames.in_dum,get_tensorvars(Kop.vars.in,Kop.deps.in,'dum')),[msg,', varnames.in_dum']);
    pv_ref = Kop.pvars.in;
    pv = get_tensorvars(Kop.vars.in,Kop.deps.in,{'dum','pvar'});
    chk(isequal(size(pv),size(pv_ref)) && all(cellfun(@(a,b) isequal(size(a),size(b)) && isequal(pvar2varname(a),pvar2varname(b)),pv,pv_ref)),[msg,', pvars.in']);
end
% apply_op for a tensor-product output (otype(1)=2): x(s1,s2) o x(s1,s2) is
% a function of two copies of s1 and of s2. The polynomials of the
% Kronecker product are compared by their variables
polyopvar x(s1,s2)
xx = otimes(x,x,2);
Y = eval_polyopvar(xx,{s1+2*s2});
chk(isequal(sort(Y.varname(:)),sort(get_tensorvars({'s1','s2'},2,'full')')),'apply_op, tensor product x o x');
% The classes take their names from get_tensorvars, so also compare with
% the names written out, which pin the convention
chk(isequal(sort(Y.varname(:)),{'s1_1';'s1_2';'s2_1';'s2_2'}),'apply_op, tensor product x o x, names written out');
Kop = intop({1},{[0,3,0,3],[0,2]},{'s1','s2'},[0,1]);
chk(isequal(Kop.varnames.in,{{'s1_1','s1_2'},{'s2_1'}}) && isequal(Kop.varnames.in_dum,{{'s1_dum1','s1_dum2'},{'s2_dum1'}}),'intop varnames, written out');
chk(isequal(pvar2varname(Kop.pvars.in{1}),{'s1_dum1';'s1_dum2'}) && isequal(pvar2varname(Kop.pvars.in{2}),{'s2_dum1'}),'intop pvars, written out');
% x(s1) o y(s2): a single copy of s1 and s2, numbered as well (apply_poly
% renames them back to s1, s2 as depmat.out<=1; eval_polyopvar does not)
polyopvar u(s1) w(s2)
uw = otimes(u,w,2);
Y = eval_polyopvar(uw,{2*s1;3*s2});
chk(isequal(sort(Y.varname(:)),sort(get_tensorvars({'s1','s2'},1,'full')')),'apply_op, tensor product u(s1) o w(s2)');
% apply_op of intop takes a function of the copies s_j, or of the bare s
% only if the functional has one copy of s; values are the integrals over
% [0,1] in each copy
K1 = intop({1},{[0,2]},{'s1'},[0,1]);
chk(abs(double(apply_op(K1,s1))-1/2)<1e-12 && abs(double(apply_op(K1,polynomial({'s1_1'})))-1/2)<1e-12,'intop degree 1: s1 and s1_1');
K2 = intop({1},{[0,3,0,3]},{'s1'},[0,1]);
s11 = polynomial({'s1_1'});     s12 = polynomial({'s1_2'});
chk(abs(double(apply_op(K2,s11*s12^2))-1/6)<1e-12,'intop degree 2: copies s1_1, s1_2');
chk(abs(double(apply_op(K2,{s1,s1^2}))-1/6)<1e-12,'intop degree 2: factors in s1');
expect_error(@() apply_op(K2,s1),'several copies','intop degree 2: bare s1');
expect_error(@() apply_op(K2,s1*s12),'several copies','intop degree 2: bare s1 with s1_2');
K12 = intop({1},{[0,2],[0,3,0,3]},{'s1','s2'},[0,1]);
chk(abs(double(apply_op(K12,s1*polynomial({'s2_1'})))-1/4)<1e-12,'intop degrees [1,2]: bare s1, copy s2_1');
expect_error(@() apply_op(K12,s1*s2),'several copies','intop degrees [1,2]: bare s2');
disp('Test_get_tensorvars Part 3 (classes) passed.')

%% Part 4: invalid inputs
expect_error(@() get_tensorvars('s1',-1),'nonnegative integers','negative degree');
expect_error(@() get_tensorvars('s1',1.5),'nonnegative integers','fractional degree');
expect_error(@() get_tensorvars('s1','a'),'nonnegative integers','char degree');
expect_error(@() get_tensorvars({'s1','s2'},[1,2,3]),'Number of degrees','too many degrees');
expect_error(@() get_tensorvars({'s1','s2','s3'},[1,2]),'Number of degrees','too few degrees');
expect_error(@() get_tensorvars(1,1),'Variables must be','numeric variables');
expect_error(@() get_tensorvars({'s1',2},1),'Variables must be','cell with a number');
expect_error(@() get_tensorvars(s1^2,1),'Variables must be','polynomial that is not a pvar');
expect_error(@() get_tensorvars([s1;s1+s2],1),'Variables must be','vector with an element that is not a pvar');
expect_error(@() get_tensorvars('s1',1,'foo'),'Unknown option','unknown option');
expect_error(@() get_tensorvars('s1',1,{'full','dummy'}),'Unknown option','unknown option in a cell');
expect_error(@() get_tensorvars('s1',1,1),'Options must be','numeric options');
expect_error(@() get_tensorvars('s1',1,{'pvar','cellstr'}),'cannot both','options pvar and cellstr');
expect_error(@() get_tensorvars(s1,1,{'polynomial','cell'}),'cannot both','options polynomial and cell');
disp('Test_get_tensorvars Part 4 (invalid inputs) passed.')

disp('Test_get_tensorvars passed.')

function chk(cond,msg)
% The condition must hold
if ~cond
    error("Test_get_tensorvars failed: %s.",msg)
end
end

function expect_error(fun,pat,msg)
% The call must be rejected, with an error message containing pat
try
    fun();
catch ME
    if ~contains(ME.message,pat)
        error("Test_get_tensorvars failed: %s was rejected with unexpected error: %s",msg,ME.message)
    end
    return
end
error("Test_get_tensorvars failed: %s was not rejected.",msg)
end
