% Test_tensopvar_subsasgn.m -- tests @tensopvar/subsasgn.m 'case {}':
% setting a single term/factor Cop{i,j} from a 'sopvar' or 'intop'
% object, and a block Cop{idcs1,idcs2} from a 'tensopvar' object. The
% result must be validated (via @tensopvar/initialize) against the rest
% of the operator: errors on a dimension/variable/domain mismatch,
% out-of-range indices, wrong-class right-hand sides, mismatched block
% shapes and chained assignment (C{i,j}.field=...). The stored objects
% are compared with the assigned ones (isequal for plain-data fields, the
% difference of the sparse coefficients for intop.params).
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Every factor depends on all N input variables and has no output
%   variables (vars.out=cell(1,0)), so any column is eligible for an
%   'intop' (a pure functional) replacement, and a matching replacement
%   operator can be built from the vars/dom/size of an existing factor
%   Cop.ops{1,jc}.

ntrials = 8;
for trial = 1:ntrials

d = randi(2)+1;             % 2 or 3 factors
nterms = 2;                 % >=2 terms, to exercise cross-term consistency
N = randi(3)+1;              % 2--4 input variables
m = randi(3)+1;               % common row dimension of every factor
n_j = randi(3)+1;             % common column dimension of every factor
dnsty = 0.5;

dims = [repmat(m,d,1), repmat(n_j,d,1)];
vars = struct('in',N,'out',0);
dom = [0,1];
depmat = struct('in',true(d,N),'out',false(d,0));

Cop = rand_tensopvar(dims,vars,dom,depmat,nterms,[2,2],2,dnsty);
if ~isvalid(Cop)
    error("Trial %d: test setup produced an invalid tensopvar.",trial)
end

jc = d;             % column (factor position) to replace
i_target = 2;       % term (row) to replace; must be ~=1, see dims check below

vars_jc = Cop.ops{1,jc}.vars.in;
dom_jc = Cop.ops{1,jc}.dom.in;
szj = size(Cop.ops{1,jc});
vars_jc_struct = struct();
vars_jc_struct.in = vars_jc;
vars_jc_struct.out = cell(1,0);
dom_jc_struct = struct();
dom_jc_struct.in = dom_jc;
dom_jc_struct.out = zeros(0,2);

%%%% Positive: Cop{i,j} = <same-structure 'sopvar'>
Cop2 = Cop;
new_op = rand_sopvar(szj,vars_jc_struct,dom_jc_struct,2,dnsty);
Cop2{i_target,jc} = new_op;
if ~isequal(Cop2.ops{i_target,jc},new_op)
    error("Trial %d: Cop{i,j}='sopvar' did not store the new operator.",trial)
end
for i=1:nterms
    for j=1:d
        if i==i_target && j==jc
            continue
        elseif ~isequal(Cop2.ops{i,j},Cop.ops{i,j})
            error("Trial %d: Cop{i,j}='sopvar' assignment altered an unrelated factor.",trial)
        end
    end
end
if ~Cop2.chkval || ~isvalid(Cop2)
    error("Trial %d: Cop{i,j}='sopvar' assignment left the object improperly validated.",trial)
end

%%%% Positive: Cop{i,j} = <same-structure 'intop'>
new_Kop = rand_intop(szj,vars_jc,dom_jc,1,1);
Cop3 = Cop;
Cop3{i_target,jc} = new_Kop;
stored_Kop = Cop3.ops{i_target,jc};
tol = 1e-10;
if ~isa(stored_Kop,'intop') || ~isequal(stored_Kop.intmat,new_Kop.intmat) ...
        || ~isequal(stored_Kop.vars,new_Kop.vars) || ~isequal(stored_Kop.dom,new_Kop.dom)
    error("Trial %d: Cop{i,j}='intop' did not store the new functional.",trial)
end
for k=1:numel(new_Kop.params)
    errk = stored_Kop.params{k}-new_Kop.params{k};      % sparse coefficients
    if nnz(errk)>0 && max(max(abs(errk)))>tol
        error("Trial %d: Cop{i,j}='intop' stored a functional with an altered kernel.",trial)
    end
end
if ~Cop3.chkval || ~isvalid(Cop3)
    error("Trial %d: Cop{i,j}='intop' assignment left the object improperly validated.",trial)
end

%%%% Positive: Cop{idcs1,jc} = <same-structure 'tensopvar' block>
idcs1 = 1:nterms;
RHS_block = tensopvar();
RHS_block.ops = cell(nterms,1);
for i=1:nterms
    RHS_block.ops{i,1} = -Cop.ops{i,jc};
end
Cop4 = Cop;
Cop4{idcs1,jc} = RHS_block;
for i=1:nterms
    if ~isequal(Cop4.ops{i,jc},RHS_block.ops{i,1})
        error("Trial %d: Cop{idcs1,jc}='tensopvar' block assignment failed for term %d.",trial,i)
    end
end
if ~Cop4.chkval || ~isvalid(Cop4)
    error("Trial %d: Cop{idcs1,jc}='tensopvar' block assignment left the object improperly validated.",trial)
end

%%%% Negative: mismatched matrix dimensions should error
bad_op = rand_sopvar([szj(1)+1,szj(2)],vars_jc_struct,dom_jc_struct,2,dnsty);
threw = false;
try
    Cop5 = Cop;
    Cop5{i_target,jc} = bad_op; %#ok<NASGU>
catch
    threw = true;
end
if ~threw
    error("Trial %d: mismatched-dimension replacement should have errored.",trial)
end

%%%% Negative: wrong-class right-hand side for a single element
threw = false;
try
    Cop6 = Cop;
    Cop6{i_target,jc} = 5; %#ok<NASGU>
catch ME
    threw = true;
    if ~contains(ME.message,"'sopvar', 'intop', or 'tensopvar'")
        error("Trial %d: wrong-class scalar assignment errored with unexpected message: %s",trial,ME.message)
    end
end
if ~threw
    error("Trial %d: wrong-class scalar assignment should have errored.",trial)
end

%%%% Negative: wrong-class right-hand side for a block ('sopvar' instead of 'tensopvar')
threw = false;
try
    Cop7 = Cop;
    Cop7{idcs1,jc} = Cop.ops{1,jc}; %#ok<NASGU>
catch ME
    threw = true;
    if ~contains(ME.message,"must be specified as 'tensopvar'")
        error("Trial %d: wrong-class block assignment errored with unexpected message: %s",trial,ME.message)
    end
end
if ~threw
    error("Trial %d: wrong-class block assignment should have errored.",trial)
end

%%%% Negative: mismatched block shape
threw = false;
try
    Cop8 = Cop;
    RHS_bad_shape = tensopvar();
    RHS_bad_shape.ops = RHS_block.ops(1,1);    % 1x1, but idcs1 selects nterms rows
    Cop8{idcs1,jc} = RHS_bad_shape; %#ok<NASGU>
catch ME
    threw = true;
    if ~contains(ME.message,'must match the number of indices')
        error("Trial %d: mismatched block shape errored with unexpected message: %s",trial,ME.message)
    end
end
if ~threw
    error("Trial %d: mismatched block shape should have errored.",trial)
end

%%%% Negative: out-of-range index
threw = false;
try
    Cop9 = Cop;
    Cop9{nterms+1,jc} = Cop.ops{1,jc}; %#ok<NASGU>
catch ME
    threw = true;
    if ~contains(ME.message,'cannot exceed number of components')
        error("Trial %d: out-of-range row index errored with unexpected message: %s",trial,ME.message)
    end
end
if ~threw
    error("Trial %d: out-of-range row index should have errored.",trial)
end

%%%% Negative: chained assignment is not supported
threw = false;
try
    Cop10 = Cop;
    Cop10{1,jc}.dims = szj; %#ok<NASGU>
catch ME
    threw = true;
    if ~contains(ME.message,'Chained assignment')
        error("Trial %d: chained assignment errored with unexpected message: %s",trial,ME.message)
    end
end
if ~threw
    error("Trial %d: chained assignment should have errored.",trial)
end

end

fprintf('Test_tensopvar_subsasgn passed (%d trials).\n',ntrials)
