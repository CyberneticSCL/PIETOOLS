function [Pop,info] = rand_tensopmat(p,d,N,M,dom,blk_dnsty,dnsty,info_in,prefix_in,prefix_out,otype1_opt,nterms_max,degs_max)
% [POP,INFO] = RAND_TENSOPMAT(P,D,N,M,DOM,BLK_DNSTY,DNSTY,INFO_IN,...
%   PREFIX_IN,PREFIX_OUT,OTYPE1_OPT,NTERMS_MAX,DEGS_MAX)
% returns a random p x d 'tensopmat' object, with each populated block a
% single-factor 'tensopvar' object built by rand_tensopvar, with the
% factors declared as random 'sopvar' objects (rand_sopvar).
%
% INPUTS
% - p:          number of block rows. Ignored if 'info_in' has the field
%               'sel_out';
% - d:          number of block columns. Ignored if 'info_in' has the field
%               'sel_in';
% - N:          number of input variables in the pool shared by all blocks,
%               named [prefix_in,'1'],...,[prefix_in,num2str(N)]. Each
%               block column depends on a random nonempty subset of these.
%               Ignored if 'info_in' has the field 'vars_in_pool';
% - M:          number of output variables in the pool shared by all
%               blocks, named with prefix_out as the input variables. Each
%               block row depends on a random nonempty subset. Ignored if
%               'info_in' has the field 'vars_out_pool';
% - dom:        1 x 2 array specifying the domain of all variables in the
%               pools generated here;
% - blk_dnsty:  (optional) scalar in [0,1], probability that block (i,j) is
%               populated. Defaults to 0.7. At least one block per row and
%               per column is always populated;
% - dnsty:      (optional) scalar in (0,1), density of the sparse
%               coefficient matrices of the 'sopvar' objects underlying
%               each block. Defaults to [], for which rand_sopvar chooses
%               the density;
% - info_in:    (optional) struct, to reuse the variable pools, variable
%               selection, domains and block sizes of an earlier call (e.g.
%               to build a second operator with a matching input or output
%               function space for horzcat/vertcat/mtimes tests). Defaults
%               to struct(). Each of the following fields is used if
%               present:
%                 'vars_in_pool', 'dom_in_pool':  names and domains of the
%                       input variable pool;
%                 'sel_in', 'dims_in':  for each block column, the indices
%                       in the pool of the input variables on which it
%                       depends, and its column dimension. Then d is taken
%                       from numel(dims_in);
%                 'vars_out_pool', 'dom_out_pool', 'sel_out', 'dims_out':
%                       the same for the output side. Then p is taken from
%                       numel(dims_out);
% - prefix_in:  (optional) char, prefix of the names of the generated input
%               variables. Defaults to 's';
% - prefix_out: (optional) char, prefix of the names of the generated
%               output variables. Defaults to 't';
% - otype1_opt: (optional) 1, 2 or 'rand', specifying otype(1) of the
%               'tensopvar' in each block (otype(2) is always 2). Defaults
%               to 1:
%                 1       - every block has otype=[1,2];
%                 2       - every block has otype=[2,2];
%                 'rand'  - each block independently draws otype(1) in
%                           {1,2} uniformly;
% - nterms_max: (optional) positive integer. Each block independently draws
%               nterms=randi(nterms_max) terms (nterms x 1 cell of
%               factors, summed) for its 'tensopvar'. Defaults to 1;
% - degs_max:   (optional) positive integer. Each block independently draws
%               degs=randi(degs_max) as the maximal degree of the monomials
%               in each variable on which its 'sopvar' kernels depend.
%               Defaults to 1;
%
% OUTPUTS
% - Pop:        random 'tensopmat' object with p x d blocks, each block
%               empty or a 'tensopvar' of size dims_out(i) x dims_in(j);
% - info:       struct with the fields 'vars_in_pool', 'dom_in_pool',
%               'sel_in', 'dims_in', 'vars_out_pool', 'dom_out_pool',
%               'sel_out', 'dims_out' (see INFO_IN), describing the input
%               and output side used to build Pop. Passing it (or parts of
%               it) as 'info_in' to a later call forces a compatible
%               function space;
%
% NOTES
% - Blocks have a single factor (size(ops,2)==1);
% - Blocks with otype(2)==1 are not built, as apply_op for 'tensopvar'
%   requires otype(2)==2 (no pointwise products of input monomial
%   factors);
% - A block with otype(1)==2 as right operand of a product with a
%   'tensopmat' or 'tensopvar' may give an error in @tensopvar/mtimes.m;
%   pass otype1_opt=1 for such operands;
% - For otype(1)==2, apply_op appends "_<k>" to the names of the output
%   variables (the coefficients are unchanged); cross-comparisons must
%   normalize these, see 'strip_kron_suffix' in Test_tensopmat_mtimes.m;
% - '()'-indexing is not possible for blocks with otype(2)==2 and more
%   than one column (or otype(1)==2 and more than one row), due to the
%   "'()'-type subsref ... not supported for Kronecker products" error of
%   @tensopvar/subsref.m.

% DJ (CC), 10/15/2026: Initial coding;

if nargin<6 || isempty(blk_dnsty)
    blk_dnsty = 0.7;
end
if nargin<7
    dnsty = [];
end
if nargin<8
    info_in = struct();
end
if nargin<9 || isempty(prefix_in)
    prefix_in = 's';
end
if nargin<10 || isempty(prefix_out)
    prefix_out = 't';
end
if nargin<11 || isempty(otype1_opt)
    otype1_opt = 1;
elseif ~isequal(otype1_opt,1) && ~isequal(otype1_opt,2) && ~isequal(otype1_opt,'rand')
    error("otype1_opt must be 1, 2, or 'rand'.")
end
if nargin<12 || isempty(nterms_max)
    nterms_max = 1;
end
if nargin<13 || isempty(degs_max)
    degs_max = 1;
end

if isfield(info_in,'vars_in_pool')
    vars_in_pool = info_in.vars_in_pool;
    dom_in_pool = info_in.dom_in_pool;
else
    vars_in_pool = arrayfun(@(i) [prefix_in,num2str(i)],1:N,'UniformOutput',false);
    dom_in_pool = repmat(dom,N,1);
end
if isfield(info_in,'sel_in')
    sel_in = info_in.sel_in;
    dims_in = info_in.dims_in;
    d = numel(dims_in);
else
    Npool = numel(vars_in_pool);
    sel_in = cell(d,1);
    for j=1:d
        nsel = randi(Npool);
        sel_in{j} = sort(randperm(Npool,nsel));
    end
    dims_in = randi(2,d,1);
end
if isfield(info_in,'vars_out_pool')
    vars_out_pool = info_in.vars_out_pool;
    dom_out_pool = info_in.dom_out_pool;
else
    vars_out_pool = arrayfun(@(i) [prefix_out,num2str(i)],1:M,'UniformOutput',false);
    dom_out_pool = repmat(dom,M,1);
end
if isfield(info_in,'sel_out')
    sel_out = info_in.sel_out;
    dims_out = info_in.dims_out;
    p = numel(dims_out);
else
    Mpool = numel(vars_out_pool);
    sel_out = cell(p,1);
    for i=1:p
        nsel = randi(Mpool);
        sel_out{i} = sort(randperm(Mpool,nsel));
    end
    dims_out = randi(2,p,1);
end

% Randomly populate blocks, guaranteeing at least one nonempty block per
% row and per column
ops = cell(p,d);
for i=1:p
    for j=1:d
        if rand()<=blk_dnsty
            ops{i,j} = build_block(i,j);
        end
    end
end
for j=1:d
    if all(cellfun(@isempty,ops(:,j)))
        i = randi(p);
        ops{i,j} = build_block(i,j);
    end
end
for i=1:p
    if all(cellfun(@isempty,ops(i,:)))
        j = randi(d);
        ops{i,j} = build_block(i,j);
    end
end

Pop = tensopmat(ops);

info = struct('vars_in_pool',{vars_in_pool},'dom_in_pool',dom_in_pool,'sel_in',{sel_in},'dims_in',dims_in,...
               'vars_out_pool',{vars_out_pool},'dom_out_pool',dom_out_pool,'sel_out',{sel_out},'dims_out',dims_out);

    function op_ij = build_block(i,j)
        vars_ij = struct('in',{vars_in_pool(sel_in{j})},'out',{vars_out_pool(sel_out{i})});
        dom_ij = struct('in',dom_in_pool(sel_in{j},:),'out',dom_out_pool(sel_out{i},:));
        depmat_ij = struct('in',true(1,numel(sel_in{j})),'out',true(1,numel(sel_out{i})));
        dims_ij = [dims_out(i),dims_in(j)];
        % otype(2) is always 2 (see NOTES)
        if isequal(otype1_opt,'rand')
            otype1_ij = randi(2);
        else
            otype1_ij = otype1_opt;
        end
        nterms_ij = randi(max(nterms_max,1));
        degs_ij = randi(max(degs_max,1));
        op_ij = rand_tensopvar(dims_ij,vars_ij,dom_ij,depmat_ij,nterms_ij,[otype1_ij,2],degs_ij,dnsty);
    end

end
