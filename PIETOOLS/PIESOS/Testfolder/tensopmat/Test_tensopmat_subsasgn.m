% Test_tensopmat_subsasgn.m -- tests @tensopmat/subsasgn.m 'case {}':
% setting a single block Cop{i,j} directly from a 'tensopvar' or 'intop'
% object, setting a block range Cop{idcs1,idcs2} from a 'tensopmat'
% object, and checking that the result is validated (via
% @tensopmat/initialize) against the rest of the row(s)/column(s) the new
% block(s) belong to, erroring on a dimension/variable/domain mismatch;
% also rejection of out-of-range indices, wrong-class right-hand sides,
% mismatched block shapes, and chained assignment (C{i,j}.field=...).

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - A mismatched-dimension replacement is only guaranteed to be caught if
%   the target block (i0,j0) shares its row with another populated block
%   (initialize.m takes a row's dimension and output variable pool from
%   the first block it finds; a sole occupant of its row and column may
%   have any shape);
% - The 'intop' sub-test needs (i0,j0) NOT to be the first populated block
%   of its row: 'intop.vars.out' is always empty, so it would redefine the
%   row's output-variable pool to empty and orphan its siblings'
%   variables. 'find_sibling_block' returns a row with >=2 populated
%   blocks and its LAST populated column, and regenerates (20 attempts)
%   otherwise.

ntrials = 8;
for trial = 1:ntrials

p = randi([3,4]);  d = randi([3,4]);
N = randi([2,3]);  M = randi([2,3]);
dom = [0,1];

[Cop,i0,j0] = find_sibling_block(p,d,N,M,dom);
if ~isvalid(Cop)
    error("Trial %d: test setup produced an invalid tensopmat.",trial)
end

blk = Cop.ops{i0,j0};
vars_blk = blk.vars;    dom_blk = blk.dom;    szblk = size(blk);

%%%% Positive: Cop{i,j} = <same-structure 'tensopvar'>
Cop2 = Cop;
new_block = rand_tensopvar(szblk,...
    struct('in',{vars_blk.in},'out',{vars_blk.out}),...
    struct('in',dom_blk.in,'out',dom_blk.out),...
    struct('in',true(1,numel(vars_blk.in)),'out',true(1,numel(vars_blk.out))),...
    1,[1,2],1,[]);
Cop2{i0,j0} = new_block;
if ~isequal(Cop2.ops{i0,j0},new_block)
    error("Trial %d: Cop{i,j}='tensopvar' did not store the new block.",trial)
end
for i=1:p
    for j=1:d
        if i==i0 && j==j0
            continue
        elseif ~isequal(Cop2.ops{i,j},Cop.ops{i,j})
            error("Trial %d: Cop{i,j}='tensopvar' assignment altered an unrelated block.",trial)
        end
    end
end
if ~Cop2.chkval || ~isvalid(Cop2)
    error("Trial %d: Cop{i,j}='tensopvar' assignment left the object improperly validated.",trial)
end

%%%% Positive: Cop{i,j} = <same-structure 'intop'>
% degs=1 matches the degree-1 input dependency of every rand_tensopmat block
% (initialize.m hardcodes deps_in_j=ones(...)), consistent with siblings in column j0
new_Kop = rand_intop(szblk,vars_blk.in,dom_blk.in,1,1);
Cop3 = Cop;
Cop3{i0,j0} = new_Kop;
stored_Kop = Cop3.ops{i0,j0};
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

%%%% Positive: Cop{idcs1,idcs2} = <same-structure 'tensopmat' block>
% Negate a random 2x2 sub-grid of (possibly empty) blocks: same structure,
% changed content
idcs1 = sort(randperm(p,min(2,p)));
idcs2 = sort(randperm(d,min(2,d)));
RHS_block = tensopmat();
RHS_block.ops = cell(numel(idcs1),numel(idcs2));
for ii=1:numel(idcs1)
    for jj=1:numel(idcs2)
        blk_ij = Cop.ops{idcs1(ii),idcs2(jj)};
        if ~isempty(blk_ij)
            blk_ij = -blk_ij;
        end
        RHS_block.ops{ii,jj} = blk_ij;
    end
end
Cop4 = Cop;
Cop4{idcs1,idcs2} = RHS_block;
for ii=1:numel(idcs1)
    for jj=1:numel(idcs2)
        if ~isequal(Cop4.ops{idcs1(ii),idcs2(jj)},RHS_block.ops{ii,jj})
            error("Trial %d: Cop{idcs1,idcs2}='tensopmat' block assignment failed for (%d,%d).",trial,ii,jj)
        end
    end
end
if ~Cop4.chkval || ~isvalid(Cop4)
    error("Trial %d: Cop{idcs1,idcs2}='tensopmat' block assignment left the object improperly validated.",trial)
end

%%%% Negative: mismatched matrix dimensions should error
bad_block = rand_tensopvar([szblk(1)+1,szblk(2)+1],...
    struct('in',{vars_blk.in},'out',{vars_blk.out}),...
    struct('in',dom_blk.in,'out',dom_blk.out),...
    struct('in',true(1,numel(vars_blk.in)),'out',true(1,numel(vars_blk.out))),...
    1,[1,2],1,[]);
threw = false;
try
    Cop5 = Cop;
    Cop5{i0,j0} = bad_block; %#ok<NASGU>
catch
    threw = true;
end
if ~threw
    error("Trial %d: mismatched-dimension replacement should have errored.",trial)
end

%%%% Negative: wrong-class right-hand side for a single block
threw = false;
try
    Cop6 = Cop;
    Cop6{i0,j0} = 'a'; %#ok<NASGU>
catch ME
    threw = true;
    if ~contains(ME.message,"'tensopvar' or 'intop'")
        error("Trial %d: wrong-class single-block assignment errored with unexpected message: %s",trial,ME.message)
    end
end
if ~threw
    error("Trial %d: wrong-class single-block assignment should have errored.",trial)
end

%%%% Negative: wrong-class right-hand side for a block ('tensopvar' instead of 'tensopmat')
threw = false;
try
    Cop7 = Cop;
    Cop7{idcs1,idcs2} = Cop.ops{i0,j0}; %#ok<NASGU>
catch ME
    threw = true;
    if ~contains(ME.message,"must be specified as 'tensopmat'")
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
    RHS_bad_shape = tensopmat();
    RHS_bad_shape.ops = RHS_block.ops(1,1);    % 1x1, but idcs1/idcs2 select more
    Cop8{idcs1,idcs2} = RHS_bad_shape; %#ok<NASGU>
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
    Cop9{p+1,j0} = Cop.ops{i0,j0}; %#ok<NASGU>
catch ME
    threw = true;
    if ~contains(ME.message,'cannot exceed number of row-blocks')
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
    Cop10{i0,j0}.dims = szblk; %#ok<NASGU>
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

fprintf('Test_tensopmat_subsasgn passed (%d trials).\n',ntrials)


function [Cop,i0,j0] = find_sibling_block(p,d,N,M,dom)
% Draw random 'tensopmat' objects (high block density) until a row has >=2
% populated blocks; i0 = that row, j0 = its LAST populated column (see NOTES)
for attempt=1:20
    Cop = rand_tensopmat(p,d,N,M,dom,0.9);
    for ii=1:p
        row_nonempty = find(~cellfun(@isempty,Cop.ops(ii,:)));
        if numel(row_nonempty)>=2
            i0 = ii;  j0 = row_nonempty(end);
            return
        end
    end
end
error("Could not generate a 'tensopmat' with a row containing >=2 populated blocks after 20 attempts.")
end
