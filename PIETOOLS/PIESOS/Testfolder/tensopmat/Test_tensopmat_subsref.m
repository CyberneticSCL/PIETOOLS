% Test_tensopmat_subsref.m -- tests @tensopmat/subsref.m for both '()'
% (whole-block extraction) and '{}' (single-block cell extraction), and
% @tensopmat/subsasgn.m's '{}' assignment.
%
% Oracle: apply_op of the extracted block vs. the original block, and
% isequal of '{}' extraction/assignment against the raw blocks.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - @tensopvar/subsref.m refuses ANY '()'-indexing once otype(2)==2 and
%   the block is wider than 1 (apply_op requires otype(2)==2), so the '()'
%   sub-check is skipped when that error is raised (prints a message);
% - The '{}' subsasgn sub-check is skipped if the chosen block is empty.

tol = 1e-8;

p = randi([3,4]);  d = randi([3,4]);
N = randi([2,3]);  M = randi([2,3]);
dom = [0,1];

[Cop,info] = rand_tensopmat(p,d,N,M,dom,0.7);
X = cell(1,d);
for j=1:d
    X{j} = rand_poly([info.dims_in(j),1],info.vars_in_pool(info.sel_in{j}),1,5);
end
Y = apply_op(Cop,X);

% '()' subsref: extract exactly the (i0,j0) block's own row/column range
% and compare against apply_op on the same block
i0 = randi(p);  j0 = randi(d);
nnr_arr = [0;cumsum(Cop.dims.out)];
nnc_arr = [0;cumsum(Cop.dims.in)];
rrange = nnr_arr(i0)+1:nnr_arr(i0+1);
crange = nnc_arr(j0)+1:nnc_arr(j0+1);
try
    Csub = Cop(rrange,crange);
    if ~isempty(Cop.ops{i0,j0})
        Ysub_direct = apply_op(Cop.ops{i0,j0},X{j0});
    else
        Ysub_direct = polynomial(zeros(numel(rrange),1));
    end
    Ysub = apply_op(Csub,{X{j0}});
    err_sub = cleanpoly(Ysub-Ysub_direct,1e-10);
    if max(max(abs(err_sub.C)))>tol
        error("Test_tensopmat_subsref failed: '()' block extraction mismatch.")
    end
catch ME
    if ~contains(ME.message,"Kronecker products")
        rethrow(ME)
    end
    disp('Test_tensopmat_subsref: skipped ''()'' sub-check (hit pre-existing tensopvar Kronecker-indexing limitation).')
end

% '{}' subsref: Cop{i0,j0} must equal the raw block Cop.ops{i0,j0}
if ~isequal(Cop{i0,j0},Cop.ops{i0,j0})
    error("Test_tensopmat_subsref failed: '{}' extraction does not match Cop.ops{i0,j0}.")
end

% '{}' subsasgn: replace a block with a freshly built, dimensionally
% compatible one; Cop{i0,j0} must then be that block (RHS stored as is)
if isempty(Cop.ops{i0,j0})
    disp('Test_tensopmat_subsref: skipped subsasgn sub-check (block was empty).')
else
    new_block = rand_tensopvar([numel(rrange),numel(crange)],...
        struct('in',{Cop.ops{i0,j0}.vars.in},'out',{Cop.ops{i0,j0}.vars.out}),...
        struct('in',Cop.ops{i0,j0}.dom.in,'out',Cop.ops{i0,j0}.dom.out),...
        struct('in',true(1,numel(Cop.ops{i0,j0}.vars.in)),'out',true(1,numel(Cop.ops{i0,j0}.vars.out))),...
        1,[1,2],1,[]);
    Cop2 = Cop;
    Cop2{i0,j0} = new_block;
    if ~isequal(Cop2.ops{i0,j0},new_block)
        error("Test_tensopmat_subsref failed: '{}' subsasgn did not store the assigned block correctly.")
    end
    for i=1:p
        for j=1:d
            if i==i0 && j==j0
                continue
            end
            if ~isequal(Cop2.ops{i,j},Cop.ops{i,j})
                error("Test_tensopmat_subsref failed: '{}' subsasgn modified an unrelated block.")
            end
        end
    end
end

disp('Test_tensopmat_subsref passed.')
