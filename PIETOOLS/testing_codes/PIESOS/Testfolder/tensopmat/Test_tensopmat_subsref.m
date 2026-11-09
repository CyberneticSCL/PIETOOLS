% Test_tensopmat_subsref.m -- tests @tensopmat/subsref.m for both '()'
% (whole-block extraction) and '{}' (single block, or ':' and index vectors).
%
% Oracle: apply_op of the extracted block vs. the original block, and
% isequal of '{}' extraction against the raw blocks.

% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Add '{}' with ':' and index vectors (the 'tensopmat'
%   of the selected blocks), empty rows, errors; subsasgn tests moved to
%   Test_tensopmat_subsasgn;

% NOTES
% - @tensopvar/subsref.m refuses ANY '()'-indexing once otype(2)==2 and
%   the block is wider than 1 (apply_op requires otype(2)==2), so the '()'
%   sub-check is skipped when that error is raised (prints a message).

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

%%%% '{}' with ':' and index vectors: the 'tensopmat' of the selected blocks.
% Oracle: apply_op of the raw blocks, summed over the selected columns
ncases = 0;
for trial=1:6
    p = randi([3,4]);  d = randi([3,4]);
    [Cop,info] = rand_tensopmat(p,d,randi([2,3]),randi([2,3]),dom,0.7);
    X = cell(1,d);
    for j=1:d
        X{j} = rand_poly([info.dims_in(j),1],info.vars_in_pool(info.sel_in{j}),1,5);
    end
    sels = {':',':'; randi(p),':'; ':',randi(d); randperm(p,2),randperm(d,2); randperm(p),1:d; 1:p,randperm(d,1); [1,1],[d,1]};
    for k=1:size(sels,1)
        I = sels{k,1};      J = sels{k,2};
        rows = expand_idx(I,p);     cols = expand_idx(J,d);
        lbl = sprintf('trial %d, C{%s,%s}',trial,idx_str(I),idx_str(J));
        Csub = Cop{I,J};
        if ~isa(Csub,'tensopmat') || ~isequal(size(Csub.ops),[numel(rows),numel(cols)])
            error("Test_tensopmat_subsref failed (%s): not a 'tensopmat' of the selected blocks.",lbl)
        elseif ~isequal(Csub.ops,Cop.ops(rows,cols))
            error("Test_tensopmat_subsref failed (%s): blocks differ from Cop.ops(rows,cols), or their order.",lbl)
        elseif ~isequal(Csub.dims.out(:),Cop.dims.out(rows(:))) || ~isequal(Csub.dims.in(:),Cop.dims.in(cols(:)))
            error("Test_tensopmat_subsref failed (%s): dims of the selected blocks.",lbl)
        elseif ~isequal(Csub.depmat.out,Cop.depmat.out(rows,:)) || ~isequal(Csub.depmat.in,Cop.depmat.in(cols,:))
            error("Test_tensopmat_subsref failed (%s): depmat of the selected blocks.",lbl)
        elseif ~isequal(Csub.matdim,[sum(Cop.dims.out(rows)),sum(Cop.dims.in(cols))])
            error("Test_tensopmat_subsref failed (%s): matrix dimensions.",lbl)
        elseif ~Csub.chkval
            error("Test_tensopmat_subsref failed (%s): extracted object is not marked as checked, so it would be re-initialized from its blocks.",lbl)
        end
        % apply to the inputs of the selected columns, if every selected row has a block
        pop = ~cellfun(@isempty,Cop.ops(rows,cols));
        if all(any(pop,2))
            Ysub = apply_op(Csub,X(cols));
            Yexp = [];
            for i=rows(:)'
                yi = polynomial(zeros(Cop.dims.out(i),1));
                for j=cols(:)'
                    if ~isempty(Cop.ops{i,j})
                        yi = yi + apply_op(Cop.ops{i,j},X{j});
                    end
                end
                Yexp = [Yexp;yi];       %#ok<AGROW>
            end
            err = cleanpoly(Ysub-Yexp,1e-10);
            if max(max(abs(err.C)))>tol*max(1,max(max(abs(Yexp.C))))
                error("Test_tensopmat_subsref failed (%s): apply_op of the selected blocks does not match the raw blocks.",lbl)
            end
        end
        ncases = ncases+1;
    end
    % a row without any block keeps its dimensions (taken from Cop, not inferred from the blocks)
    r0 = randi(p);
    D = Cop;    D.ops(r0,:) = {[]};
    Dsub = D{r0,':'};
    if ~isa(Dsub,'tensopmat') || any(~cellfun(@isempty,Dsub.ops(:))) || ~isequal(Dsub.dims.out,Cop.dims.out(r0)) || ~isequal(Dsub.dims.in(:),Cop.dims.in(:))
        error("Test_tensopmat_subsref failed (trial %d): empty row D{r0,:}.",trial)
    end
    % a single block selected with ':' is still a 'tensopmat', C{i,j} the block itself
    [i1,j1] = find(~cellfun(@isempty,Cop.ops),1,'first');
    C1 = tensopmat(Cop.ops(i1,j1));
    if ~isa(C1{':',':'},'tensopmat') || ~isequal(C1{':',':'}.ops,C1.ops) || ~isequal(C1{1,1},Cop.ops{i1,j1})
        error("Test_tensopmat_subsref failed (trial %d): ':' on a 1x1 'tensopmat' must return a 'tensopmat'.",trial)
    end
    % chained references apply to the extracted object
    if ~isequal(Cop{':',1}.dims.in,Cop.dims.in(1)) || ~isequal(Cop{':',':'}.dims.out,Cop.dims.out)
        error("Test_tensopmat_subsref failed (trial %d): chained reference after C{':',j}.",trial)
    end
    % errors
    assert_error(@() Cop{':'},'Linear indexing',"C{':'}");
    assert_error(@() Cop{1,1,1},'At most two indices',"C{1,1,1}");
    assert_error(@() Cop{p+1,':'},'Row indices cannot exceed',"row beyond the blocks");
    assert_error(@() Cop{':',d+1},'Column indices cannot exceed',"column beyond the blocks");
    assert_error(@() Cop{[1,p+1],1},'Row indices cannot exceed',"index vector beyond the blocks");
    % 'end' evaluates to 1 for a 'tensopmat', so C{p,end} is C{p,1}: not supported (see header of subsref)
    if ~isequal(Cop{p,end},Cop.ops{p,1}) || ~isequal(Cop{end,1},Cop.ops{1,1})
        error("Test_tensopmat_subsref failed (trial %d): 'end' no longer evaluates to 1, update the header of subsref.",trial)
    end
    ncases = ncases+5;
end

fprintf('Test_tensopmat_subsref passed (%d cases of ''{}'' with '':'' and index vectors).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function idx = expand_idx(I,n)
% Indices selected by I in a dimension of n blocks
if ischar(I),   idx = 1:n;      else,   idx = reshape(I,1,[]);    end
end

function s = idx_str(I)
if ischar(I),   s = ':';        else,   s = mat2str(I);     end
end

function assert_error(f,needle,label)
try
    f();
catch err
    if ~contains(err.message,needle)
        error('Test_tensopmat_subsref failed (%s): wrong error "%s".',label,err.message)
    end
    return
end
error('Test_tensopmat_subsref failed (%s): expected an error containing "%s".',label,needle)
end
