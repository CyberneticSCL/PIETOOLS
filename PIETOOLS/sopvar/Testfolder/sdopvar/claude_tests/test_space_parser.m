function test_space_parser()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_SPACE_PARSER checks 'parse_copvar_spaces' and 'check_reserved_names'
% (sopvar/misc/conventions) against the parsers they replace:
%
% (A) SEMANTICS on random inputs in 1..4 spatial directions: the registry is
%     the sorted union of the names given, each mask row marks exactly the
%     names of its space, each domain row is the interval given for that
%     NAME (whatever form and order the domain was given in), and the
%     dimensions are the counts given, expanded from a scalar;
%
% (B) the local parsers, copied VERBATIM below from the files at b3808f0f:
%     'lpivar_cdopvar' (parse_spaces, norm_spaces, parse_dom and the
%     registry lines; 'spaces2meta_sop' holds byte-identical copies) and
%     the inline parse of 'copquadvar'. Same outputs on every input both
%     accept, same message on every input both refuse, and each documented
%     difference asserted as such;
%
% (C) the LIVE routines: 'lpivar_cdopvar', 'copquadvar' and
%     'zeros_copvar_sop' (which calls 'spaces2meta_sop') on the same inputs,
%     reading the container metadata they return, and the messages they
%     raise.
%
% Initial coding MMP, 09/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(20260930);
nchk = 0;
pvar x
prog0 = lpiprogram(x,[0 1]);
pool = {'u','s','y','t','x2'};                 % unsorted on purpose
cov = struct('rect',0,'plain',0,'empty_numeric',0,'copquadvar',0, ...
             'row',0,'array',0,'struct',0,'nv',zeros(1,6));

% ========================================== accepted inputs, random forms
for trial = 1:120
    nvt = randi(4);                             % 1..4 directions
    names = pool(randperm(numel(pool),nvt));
    M = randi(3);
    sp = cell(1,M);
    for i = 1:M
        pick = names(rand(1,nvt)<0.6);
        sp{i} = pick(randperm(numel(pick)));    % any order within a space
        if isempty(sp{i}) && rand<0.5,  sp{i} = [];     end     % R^q as []
    end
    rect = rand<0.5;
    if rect                                     % rectangular: struct forms
        N = randi(3);
        si = cell(1,N);
        for j = 1:N
            pick = names(rand(1,nvt)<0.6);
            si{j} = pick;
        end
        spaces = struct();  spaces.out = sp;    spaces.in = si;
        dims = struct('out',randi(3,M,1),'in',randi(3,N,1));
        if rand<0.3,    dims.out = 2;   end    % scalar expanded
    else
        spaces = sp;
        if M==1 && rand<0.5,    spaces = sp{1};     end     % plain cellstr
        if M==1 && iscell(spaces) && numel(spaces)==1 && rand<0.5
            spaces = spaces{1};                 % a char
        end
        dims = randi(3,M,1);
        if rand<0.3,    dims = 2;   end
    end
    % The registry as given: every name placed in some space.
    used = sp(~cellfun(@isnumeric,sp));     used = [used{:}];
    if rect,    used = [used, si{:}];   end
    [dom,dtruth,dform] = rand_domain(unique(used),names);
    [meta,so,si] = parse_copvar_spaces(dims,spaces,dom);
    cov.(dform) = cov.(dform)+1;
    cov.rect = cov.rect+rect;
    cov.plain = cov.plain+(~rect && ~iscell(spaces) || (iscellstr(spaces) && ~isempty(spaces)));
    cov.empty_numeric = cov.empty_numeric+any(cellfun(@isnumeric,sp));
    cov.nv(numel(meta.vars)+1) = cov.nv(numel(meta.vars)+1)+1;
    if isstruct(dims),  dgo = dims.out;     dgi = dims.in;
    else,               dgo = dims;         dgi = dims;
    end
    if isscalar(dgo),   dgo = repmat(dgo,numel(so),1);  end
    if isscalar(dgi),   dgi = repmat(dgi,numel(si),1);  end
    assert(isequal(meta.dim_out,dgo(:)) && isequal(meta.dim_in,dgi(:)),'dims as given');
    % (A) semantics
    allnames = [so{:},si{:}];
    assert(isequal(meta.vars,reshape(sort(unique(allnames)),1,[])),'registry');
    assert(issorted(meta.vars) && isrow(meta.vars),'registry sorted row');
    for i = 1:numel(so)
        assert(isequal(meta.space_out(i,:),ismember(meta.vars,so{i})),'mask out');
        assert(isequal(rowv(sort(so{i})),rowv(meta.vars(meta.space_out(i,:)))),'mask out names');
    end
    for j = 1:numel(si)
        assert(isequal(rowv(sort(si{j})),rowv(meta.vars(meta.space_in(j,:)))),'mask in names');
    end
    for d = 1:numel(meta.vars)
        assert(isequal(meta.dom(d,:),dtruth.(meta.vars{d})),'domain by name');
    end
    assert(iscolumn(meta.dim_out) && numel(meta.dim_out)==numel(so),'dim_out');
    assert(iscolumn(meta.dim_in) && numel(meta.dim_in)==numel(si),'dim_in');
    assert(isequal(fieldnames(meta)',{'vars','dom','space_out','space_in','dim_out','dim_in'}),'field order');
    % (B) lpivar_cdopvar's parsers
    [meta_l,so_l,si_l] = lpivar_cdopvar_parse(dims,spaces,dom);
    assert(isequal(meta,meta_l) && isequal(so,so_l) && isequal(si,si_l),'vs lpivar_cdopvar copy');
    % (C) live: zeros_copvar_sop (spaces2meta_sop) and lpivar_cdopvar
    Z = zeros_copvar_sop(dims,spaces,dom);
    assert(same_meta(Z,meta),'vs zeros_copvar_sop');
    [~,Q] = lpivar_cdopvar(prog0,dims,spaces,dom,0);
    assert(same_meta(Q,meta),'vs lpivar_cdopvar');
    nchk = nchk+12;
    % copquadvar: square forms only
    if ~rect && ~any(cellfun(@(c) isnumeric(c),sp))  % [] is not a copquadvar form
        [sp_c,vars_c,mask_c,dims_c,dom_c] = copquadvar_parse(dims,spaces,dom);
        assert(isequal(sp_c,so) && isequal(vars_c,meta.vars) && isequal(mask_c,meta.space_out) ...
            && isequal(dims_c,meta.dim_out) && isequal(dom_c,meta.dom),'vs copquadvar copy');
        opts = struct('include',{num2cell(ones(1,numel(so)))});     % one basis operator per space
        [~,P] = copquadvar(prog0,dims,spaces,dom,0,opts);
        assert(same_meta(P,meta),'vs copquadvar');
        cov.copquadvar = cov.copquadvar+1;
        nchk = nchk+2;
    end
end
fprintf(['test_space_parser coverage: %d rectangular, %d plain cellstr/char, %d with [] spaces, '...
         '%d copquadvar; domain forms row %d, array %d, struct %d; registry sizes 0..4: %s\n'], ...
        cov.rect,cov.plain,cov.empty_numeric,cov.copquadvar,cov.row,cov.array,cov.struct,mat2str(cov.nv));
assert(all([cov.rect,cov.plain,cov.empty_numeric,cov.copquadvar,cov.row,cov.array,cov.struct]>0) ...
    && all(cov.nv(2:5)>0),'a form was not exercised');

% A polynomial space and polynomial domain names.
pvar s t
dvp = [s,t];    D = [0 1;-1 2];     vn = dvp.varname;   % rows pair with vn
dstr = struct('vars',dvp,'dom',D);
[meta,so] = parse_copvar_spaces([1;2],{[],[t;s]},dstr);
assert(isequal(meta.vars,{'s','t'}) && isequal(sort(so{2}),{'s','t'}),'polynomial forms');
for d = 1:2
    assert(isequal(meta.dom(d,:),D(strcmp(vn,meta.vars{d}),:)),'polynomial domain names');
end
[meta_l,so_l] = lpivar_cdopvar_parse([1;2],{[],[t;s]},dstr);
assert(isequal(meta,meta_l) && isequal(so,so_l),'polynomial vs copy');
% R^q only: empty registry.
[meta,so] = parse_copvar_spaces(3,{},[0 1]);
assert(isempty(meta.vars) && isequal(size(meta.dom),[0 2]) && isequal(so,{cell(1,0)}) ...
    && isequal(size(meta.space_out),[1 0]),'R^q only');
nchk = nchk+3;

% ========================================================= refused inputs
sarr = struct('out',{{'s'},{'t'}});
dA = struct('vars',{{'s'},{'t'}},'dom',[0 1]);          % a struct ARRAY domain
cases = {
 % label,               dims,   spaces,                     dom,                                   copquadvar
 'struct array spaces', 1,      sarr,                       [0 1],                                 'accepts'
 'no field out',        1,      struct('in',{{{'s'}}}),     [0 1],                                 'differs'
 'space not cellstr',   1,      {{1}},                      [0 1],                                 'differs'
 'repeated name',       1,      {{'s','s'}},                [0 1],                                 'same'
 'dims count',          [1;1;1],{{},{'s'}},                 [0 1],                                 'differs'
 'dims zero',           0,      {{'s'}},                    [0 1],                                 'differs'
 'dims fraction',       1.5,    {{'s'}},                    [0 1],                                 'differs'
 'reserved _int',       1,      {{'s_int'}},                [0 1],                                 'same'
 'reserved _dum',       1,      {{},{'a','b_dum'}},         [0 1],                                 'same'
 'dom struct fields',   1,      {{'s'}},                    struct('vars',{{'s'}}),                'same'
 'dom struct array',    1,      {{'s','t'}},                dA,                                    'differs'
 'dom struct pairing',  1,      {{'s'}},                    struct('vars',{{'s','t'}},'dom',[0 1]),'same'
 'dom missing name',    1,      {{'s','t'}},                struct('vars',{{'s'}},'dom',[0 1]),    'same'
 'dom columns',         1,      {{'s'}},                    [0 1 2],                               'same'
 'dom rows',            1,      {{'s','t','u'}},            [0 1;0 1],                             'differs'
 'dom order',           1,      {{'s'}},                    [1 0],                                 'same'
 'spaces before dims',  0,      {{'s','s'}},                [0 1],                                 'same'
 'dims before names',   0,      {{'s_int'}},                [0 1],                                 'differs'
 'names before dom',    1,      {{'s_int'}},                [1 0],                                 'same'
 };
cq_live_parser = 0;     % live copquadvar messages that are the parser's
for c = 1:size(cases,1)
    [lbl,dm,spc,dmn,cq] = cases{c,:};
    m_new = errmsg(@() parse_copvar_spaces(dm,spc,dmn));
    m_old = errmsg(@() lpivar_cdopvar_parse(dm,spc,dmn));
    assert(~isempty(m_new) && strcmp(m_new,m_old),'%s: parser "%s" vs lpivar copy "%s"',lbl,m_new,m_old);
    m_live = errmsg(@() lpivar_cdopvar(prog0,dm,spc,dmn,0));
    assert(strcmp(m_new,m_live),'%s: vs live lpivar_cdopvar "%s"',lbl,m_live);
    m_zero = errmsg(@() zeros_copvar_sop(dm,spc,dmn));
    assert(strcmp(m_new,m_zero),'%s: vs live zeros_copvar_sop "%s"',lbl,m_zero);
    m_cq = errmsg(@() copquadvar_parse(dm,spc,dmn));
    switch cq
        case 'same',     assert(strcmp(m_new,m_cq),'%s: vs copquadvar "%s"',lbl,m_cq);
        case 'differs',  assert(~isempty(m_cq) && ~strcmp(m_new,m_cq),'%s: copquadvar should differ',lbl);
        case 'accepts',  assert(isempty(m_cq),'%s: copquadvar accepted this',lbl);
    end
    if ~strcmp(cq,'accepts')
        % Live copquadvar: its inline message before it adopts the parser,
        % the parser's after.
        m_cql = errmsg(@() copquadvar(prog0,dm,spc,dmn,0));
        assert(strcmp(m_cq,m_cql) || strcmp(m_new,m_cql),'%s: copquadvar live "%s"',lbl,m_cql);
        cq_live_parser = cq_live_parser + (strcmp(m_new,m_cql) && ~strcmp(m_cq,m_cql));
    end
    nchk = nchk+5;
end
fprintf('test_space_parser: live copquadvar raised the parser''s message on %d of the %d differing cases.\n', ...
        cq_live_parser,sum(strcmp(cases(:,5),'differs')));
% The copquadvar struct-array trap: it read the first element only.
[sp_c,vars_c] = copquadvar_parse(1,sarr,[0 1]);
assert(isequal(sp_c,{{'s'}}) && isequal(vars_c,{'s'}),'copquadvar struct-array trap');

% Documented differences on inputs one side accepts.
% (1) A non-cell space list: MATLAB's brace error before, copquadvar's text now.
m_new = errmsg(@() parse_copvar_spaces(1,5,[0 1]));
m_old = errmsg(@() lpivar_cdopvar_parse(1,5,[0 1]));
m_cq  = errmsg(@() copquadvar_parse(1,5,[0 1]));
assert(strcmp(m_new,m_cq) && ~strcmp(m_new,m_old) && contains(m_old,'Brace indexing'),'non-cell spaces');
% (2) [] is one R^q space for the parser and lpivar_cdopvar; copquadvar refused it.
meta = parse_copvar_spaces(2,[],[0 1]);
assert(isempty(meta.vars) && isequal(meta.dim_out,2),'[] spaces');
assert(~isempty(errmsg(@() copquadvar_parse(2,[],[0 1]))),'copquadvar refused []');
% (3) '' is one R^q space; copquadvar declared a variable named ''.
meta = parse_copvar_spaces(1,'',[0 1]);
[~,vars_c] = copquadvar_parse(1,'',[0 1]);
assert(isempty(meta.vars) && isequal(vars_c,{''}),'empty char');
% (4) struct dims: accepted; copquadvar failed inside MATLAB.
meta = parse_copvar_spaces(struct('out',2,'in',2),{{'s'}},[0 1]);
assert(isequal(meta.dim_out,2) && isequal(meta.dim_in,2),'struct dims');
assert(~isempty(errmsg(@() copquadvar_parse(struct('out',2,'in',2),{{'s'}},[0 1]))),'copquadvar struct dims');
nchk = nchk+9;

% ======================== the self-adjoint check copquadvar keeps (API doc)
sq = {
 % spaces,                                                 dims,                        copquadvar today, recommended
 struct('out',{{{'s'}}},'in',{{{'s'}}}),                    1,                           true,  true
 struct('out',{{{'s'}}}),                                   1,                           true,  true
 struct('out',{{{'s'}}},'in',{{{'t'}}}),                    1,                           false, false
 struct('out',{{{'s'}}},'in',{{'s'}}),                      1,                           false, true   % equal once normalized
 struct('out',{{{'s'}}},'in',{{{'s'}}}),                    struct('out',1,'in',2),      false, false
 };
for c = 1:size(sq,1)
    [spc,dm,ok_today,ok_rec] = sq{c,:};
    today = isempty(errmsg(@() copquadvar_parse(dm,spc,[0 1])));
    [meta,so,si] = parse_copvar_spaces(dm,spc,[0 1]);
    rec = isempty(errmsg(@() copquadvar_square_check(so,si,meta)));
    assert(today==ok_today && rec==ok_rec,'square check case %d',c);
    % Live copquadvar agrees with one of the two, before or after adoption.
    live = isempty(errmsg(@() copquadvar(prog0,dm,spc,[0 1],0)));
    assert(live==today || live==rec,'live copquadvar square case %d',c);
    if ~rec
        assert(strcmp(errmsg(@() copquadvar_square_check(so,si,meta)), ...
            'A self-adjoint operator has equal input and output spaces.'),'square message');
    end
    nchk = nchk+1;
end

% ============================================== check_reserved_names
for nm = {'s_int','a_dum','x_int'}
    m1 = errmsg(@() check_reserved_names({'q',nm{1}}));
    m2 = errmsg(@() copquadvar_parse(1,{{'q',nm{1}}},[0 1]));
    assert(~isempty(m1) && strcmp(m1,m2),'reserved message');
    nchk = nchk+1;
end
check_reserved_names({'s_int2','dum','int_s','s_dumb'});      % not a suffix
check_reserved_names(cell(1,0));
nchk = nchk+2;

fprintf('test_space_parser: %d checks passed.\n',nchk);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dom,truth,form] = rand_domain(reg,names)
% A random interval per name, and the domain argument in one of the three
% accepted forms: one row for every variable, an array in the order of the
% sorted registry REG, or a struct over NAMES (a superset of REG, plus an
% unused name), in random order. TRUTH maps each name to its interval.
nv = numel(names);
iv = sort(rand(nv,2)*4-2,2);    iv(:,2) = iv(:,2)+0.1;
truth = struct();
r = rand;
if r<0.3                                        % one row for every variable
    for d = 1:nv,   truth.(names{d}) = iv(1,:);     end
    dom = iv(1,:);      form = 'row';
    return
end
for d = 1:nv,   truth.(names{d}) = iv(d,:);     end
if r<0.6                                        % array, SORTED registry order
    reg = sort(reg);
    dom = zeros(numel(reg),2);
    for d = 1:numel(reg),   dom(d,:) = truth.(reg{d});  end
    form = 'array';
else                                            % struct, any order, extra names
    ord = randperm(nv);
    dv = [names(ord),{'zz'}];   dd = [iv(ord,:);0 1];
    dom = struct();     dom.vars = dv;      dom.dom = dd;
    form = 'struct';
end
end

function c = rowv(c)
% A name list as a row; sort() of an empty cell returns 0 x 0.
c = reshape(c,1,[]);
end

function tf = same_meta(P,meta)
tf = isequal(P.vars,meta.vars) && isequal(P.dom,meta.dom) ...
  && isequal(P.space_out,meta.space_out) && isequal(P.space_in,meta.space_in) ...
  && isequal(P.dim_out,meta.dim_out) && isequal(P.dim_in,meta.dim_in);
end

function copquadvar_square_check(sp_out,sp_in,meta)
% The check copquadvar keeps after adopting the parser (helpers_API.md).
if ~isequal(sp_in,sp_out) || ~isequal(meta.dim_in,meta.dim_out)
    error("A self-adjoint operator has equal input and output spaces.")
end
end

function m = errmsg(f)
m = '';
try
    f();
catch ME
    m = ME.message;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% lpivar_cdopvar, VERBATIM at b3808f0f (lines 126-138 and 211-215 as one
%%% function; the local functions renamed lpivar_cdopvar_<name>)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [meta,sp_out,sp_in] = lpivar_cdopvar_parse(dims,spaces,dom)
[sp_out,sp_in,d_out,d_in] = lpivar_cdopvar_parse_spaces(dims,spaces);
M = numel(sp_out);      N = numel(sp_in);
vars = reshape(unique([sp_out{:}, sp_in{:}]),1,[]);     nv = numel(vars);
is_reserved = ~cellfun(@isempty,regexp(vars,'_(int|dum)$','once'));
if any(is_reserved)
    error("Spatial variable names may not end in '_int' or '_dum'; "...
          +"'"+string(vars{find(is_reserved,1)})+"' does.")
end
dom = lpivar_cdopvar_parse_dom(dom,vars);
so = false(M,nv);       si = false(N,nv);
for i = 1:M,    so(i,:) = ismember(vars,sp_out{i});    end
for j = 1:N,    si(j,:) = ismember(vars,sp_in{j});     end
meta = struct('vars',{vars},'dom',dom,'space_out',so,'space_in',si,...
    'dim_out',d_out(:),'dim_in',d_in(:));
end

function [so,si,dout,din] = lpivar_cdopvar_parse_spaces(dims,spaces)
% Output and input space lists, and component counts.
if isstruct(spaces)
    if ~isscalar(spaces)
        error("'spaces' should be a single struct, not an array; struct() with "...
              +"a cell value returns an array, so assign the fields instead.")
    end
    if ~isfield(spaces,'out')
        error("A 'struct' spaces argument needs a field 'out'.")
    end
    so = lpivar_cdopvar_norm_spaces(spaces.out);
    if isfield(spaces,'in'),    si = lpivar_cdopvar_norm_spaces(spaces.in);
    else,                       si = so;
    end
else
    so = lpivar_cdopvar_norm_spaces(spaces);   si = so;
end
M = numel(so);      N = numel(si);
if isstruct(dims)
    dout = reshape(dims.out,[],1);
    if isfield(dims,'in'),  din = reshape(dims.in,[],1);
    else,                   din = dout;
    end
else
    dout = reshape(dims,[],1);  din = dout;
end
if isscalar(dout),  dout = repmat(dout,M,1);    end
if isscalar(din),   din = repmat(din,N,1);      end
if numel(dout)~=M || numel(din)~=N
    error("Dimensions should have one entry per space: %d output and %d input.",M,N)
end
if any([dout;din]<1) || any([dout;din]~=round([dout;din]))
    error("Component counts should be positive integers.")
end
end

function sp = lpivar_cdopvar_norm_spaces(sp)
% A plain cellstr is ONE space; otherwise a cell of cellstr, one per space.
if isempty(sp),     sp = {cell(1,0)};   return,     end
if iscellstr(sp) || ischar(sp) || isa(sp,'polynomial'),  sp = {sp};  end
sp = reshape(sp,1,[]);
for k = 1:numel(sp)
    sk = sp{k};
    if isa(sk,'polynomial'),            sk = sk.varname(:)';    end
    if ischar(sk),                      sk = {sk};              end
    if isnumeric(sk) && isempty(sk),    sk = cell(1,0);         end
    if ~iscellstr(sk)
        error("Space %d should be a 'cellstr' of variable names; an empty "...
              +"one is the finite-dimensional space R^q.",k)
    end
    sk = reshape(sk,1,[]);
    if numel(unique(sk))~=numel(sk)
        error("Space %d repeats a variable name.",k)
    end
    sp{k} = sk;
end
end

function dom = lpivar_cdopvar_parse_dom(dom,vars)
% As 'copquadvar': nv x 2 in registry order, one 1 x 2 row for all, or a
% struct pairing names with intervals.
nv = numel(vars);
if isa(dom,'struct')
    if ~isfield(dom,'vars') || ~isfield(dom,'dom')
        error("A 'struct' domain should have fields 'vars' and 'dom'.")
    end
    if ~isscalar(dom)
        error("A 'struct' domain should be a single struct, not an array; "...
              +"struct() with a cell value returns an array.")
    end
    dvars = dom.vars;
    if isa(dvars,'polynomial'),  dvars = dvars.varname(:)';   end
    if ischar(dvars),            dvars = {dvars};             end
    dvars = reshape(dvars,1,[]);
    if size(dom.dom,1)~=numel(dvars) || size(dom.dom,2)~=2
        error("A 'struct' domain should pair each name in 'vars' with a row of 'dom'.")
    end
    [tf,loc] = ismember(vars,dvars);
    if ~all(tf)
        error("No domain was given for the variable '"+string(vars{find(~tf,1)})+"'.")
    end
    dom = dom.dom(loc,:);
end
if nv==0
    dom = zeros(0,2);
    return
end
if size(dom,2)~=2
    error("Domains should be specified as an nv x 2 array.")
end
if size(dom,1)==1 && nv~=1
    dom = repmat(dom,nv,1);
elseif size(dom,1)~=nv
    error("Domains should be an nv x 2 array for the %d registry variables.",nv)
end
if any(dom(:,2)<=dom(:,1))
    error("Each domain should satisfy dom(d,1) < dom(d,2).")
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% copquadvar, VERBATIM at b3808f0f (lines 264-390, the input processing
%%% of spaces, registry, dimensions, reserved names and domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spaces,vars,mask,dims,dom] = copquadvar_parse(dims,spaces,dom)
if isa(spaces,'struct')
    if ~isfield(spaces,'out')
        error("A 'struct' space specification should have a field 'out'.")
    end
    if isfield(spaces,'in') && ~isequal(spaces.in,spaces.out)
        error("A self-adjoint operator has equal input and output spaces.")
    end
    spaces = spaces.out;
end
if ischar(spaces) || isa(spaces,'polynomial')
    spaces = {spaces};
end
if iscellstr(spaces)
    spaces = {spaces};
end
if ~iscell(spaces)
    error("Spaces should be specified as a cell of 'cellstr' objects.")
end
spaces = reshape(spaces,1,[]);      M = numel(spaces);
if M<1
    error("At least one space must be specified.")
end
for k = 1:M
    sk = spaces{k};
    if isa(sk,'polynomial'),            sk = sk.varname(:)';    end
    if ischar(sk),                      sk = {sk};              end
    if isnumeric(sk) && isempty(sk),    sk = cell(1,0);         end
    if ~iscellstr(sk)
        error("Space "+num2str(k)+" should be a 'cellstr' of variable names; "...
              +"an empty one is the finite-dimensional space R^m.")
    end
    sk = reshape(sk,1,[]);
    if numel(unique(sk))~=numel(sk)
        error("Space "+num2str(k)+" repeats a variable name.")
    end
    spaces{k} = sk;
end

vars = reshape(unique([spaces{:}]),1,[]);       nv = numel(vars);
mask = false(M,nv);
for k = 1:M
    mask(k,:) = ismember(vars,spaces{k});
end
own = cell(1,M);
for k = 1:M
    own{k} = reshape(find(mask(k,:)),1,[]);
end

dims = reshape(dims,[],1);
if isscalar(dims)
    dims = repmat(dims,M,1);
elseif numel(dims)~=M
    error("Dimensions should be a scalar or have one entry per space.")
end
if any(dims<1) || any(dims~=round(dims))
    error("Dimensions of the operator should be positive integers.")
end

vars_int = strcat(vars,'_int');
vars_dum = strcat(vars,'_dum');
is_reserved = ~cellfun(@isempty,regexp(vars,'_(int|dum)$','once'));
if any(is_reserved)
    error("Spatial variable names may not end in '_int' or '_dum'; "...
          +"'"+string(vars{find(is_reserved,1)})+"' does.")
end

if isa(dom,'struct')
    if ~isfield(dom,'vars') || ~isfield(dom,'dom')
        error("A 'struct' domain should have fields 'vars' and 'dom'.")
    end
    if ~isscalar(dom)
        error("A 'struct' domain should be a single struct, not a "...
              +num2str(numel(dom))+"-element array; build it by assigning "...
              +"the fields, since struct() with a cell value returns an array.")
    end
    dvars = dom.vars;
    if isa(dvars,'polynomial'),  dvars = dvars.varname(:)';   end
    if ischar(dvars),            dvars = {dvars};             end
    dvars = reshape(dvars,1,[]);
    if size(dom.dom,1)~=numel(dvars) || size(dom.dom,2)~=2
        error("A 'struct' domain should pair each name in 'vars' with a row "...
              +"of 'dom'.")
    end
    [tf_d,loc] = ismember(vars,dvars);
    if ~all(tf_d)
        error("No domain was given for the variable '"...
              +string(vars{find(~tf_d,1)})+"'.")
    end
    dom = dom.dom(loc,:);
end
if nv==0
    dom = zeros(0,2);
else
    if size(dom,2)~=2
        error("Domains should be specified as an nv x 2 array.")
    end
    if size(dom,1)==1 && nv~=1
        dom = repmat(dom,nv,1);
    elseif size(dom,1)~=nv
        error("Domains should be specified as an nv x 2 array for nv "...
              +"registry variables.")
    end
    if any(dom(:,2)<=dom(:,1))
        error("Each domain should satisfy dom(d,1) < dom(d,2).")
    end
end
end
