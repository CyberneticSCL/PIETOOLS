function Test_tensopvar_clean_op()
% Test_tensopvar_clean_op.m -- one local function for each case:
%   case_basic              verify @tensopvar/clean_op.m correctly removes
%                           terms (rows of Cop.ops) that are zero up to
%                           tolerance, leaves nonzero terms and the
%                           represented operator untouched, and always
%                           retains
%   case_nested             that @tensopvar/clean_op.m recurses into nested
%                           'tensopvar' and 'intop' factors, dropping a
%                           term whose nested factor is entirely zero and
%                           leaving the action of the
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_tensopvar_clean_op_nested into this test,
%   as cases;

cases = {@case_basic, @case_nested};
nfail = 0;
for i=1:numel(cases)
    tstart = tic;
    try
        cases{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(cases{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(cases{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_tensopvar_clean_op: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_tensopvar_clean_op: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- verify @tensopvar/clean_op.m correctly
% removes terms (rows of Cop.ops) that are zero up to tolerance, leaves
% nonzero terms and the represented operator untouched, and always retains
% at least one term so matdim/vars/dom/depmat stay well-defined.
%
% The oracle is the action on random inputs (apply_op) of the operator
% before and after cleaning.
%

    tol = 1e-8;

    for trial = 1:8
        rng(6400+trial);

        dims = [randi(2),randi(2); randi(2),randi(2)];
        vars = struct('in',{{'s1','s2'}},'out',{{'s3','s4'}});
        dom = struct('in',repmat([0,1],2,1),'out',repmat([0,1],2,1));
        depmat = struct('in',logical(eye(2)),'out',logical(eye(2)));
        Bop = rand_tensopvar(dims,vars,dom,depmat,3,[2,2],1,0.5);

        % Case 1: zero out one term's first factor: that term is dropped, the action is unchanged
        Bop_bad = Bop;
        Bop_bad.ops{2,1} = 0*Bop_bad.ops{2,1};
        Bop_clean = clean_op(Bop_bad);
        if size(Bop_clean.ops,1)~=2
            error("Test_tensopvar_clean_op failed: expected 2 surviving terms, got %d (trial=%d).",size(Bop_clean.ops,1),trial)
        end
        if ~isequal(Bop_clean.matdim,Bop_bad.matdim) || ~isequal(Bop_clean.vars,Bop_bad.vars) || ...
                ~isequal(Bop_clean.dom,Bop_bad.dom) || ~isequal(Bop_clean.depmat,Bop_bad.depmat)
            error("Test_tensopvar_clean_op failed: matdim/vars/dom/depmat changed after cleaning (trial=%d).",trial)
        end
        pvars = polynomial(Bop.vars.in(:));
        v = {rand_poly([dims(1,2),1],pvars(1),3,10), rand_poly([dims(2,2),1],pvars(2),3,10)};
        Y_bad = apply_op(Bop_bad,v);
        Y_clean = apply_op(Bop_clean,v);
        err = cleanpoly(Y_bad-Y_clean,1e-10);
        if max(max(abs(err.C)))>tol
            error("Test_tensopvar_clean_op failed: cleaned operator does not act identically (trial=%d).",trial)
        end

        % Case 2: every term zero: exactly 1 term retained, acts as zero map
        Bop_allzero = Bop;
        for i=1:3
            Bop_allzero.ops{i,1} = 0*Bop_allzero.ops{i,1};
        end
        Bop_allzero_clean = clean_op(Bop_allzero);
        if size(Bop_allzero_clean.ops,1)~=1
            error("Test_tensopvar_clean_op failed: all-zero case should retain exactly 1 term, got %d (trial=%d).",size(Bop_allzero_clean.ops,1),trial)
        end
        if ~isequal(Bop_allzero_clean.matdim,Bop_allzero.matdim)
            error("Test_tensopvar_clean_op failed: matdim not preserved in all-zero case (trial=%d).",trial)
        end
        Y_allzero = apply_op(Bop_allzero_clean,v);
        if max(max(abs(double(cleanpoly(Y_allzero,1e-10).C))))>tol
            error("Test_tensopvar_clean_op failed: all-zero case did not act as the zero operator (trial=%d).",trial)
        end

        % Case 3: nothing zero: no term removed
        Bop_clean2 = clean_op(Bop);
        if size(Bop_clean2.ops,1)~=3
            error("Test_tensopvar_clean_op failed: no-op case incorrectly removed a term (trial=%d).",trial)
        end
    end

    disp('Test_tensopvar_clean_op passed.')
end

%% --------------------------------------------------------------------------
function case_nested()
% case_nested -- tests that @tensopvar/clean_op.m
% recurses into nested 'tensopvar' and 'intop' factors, dropping a term
% whose nested factor is entirely zero and leaving the action of the
% represented operator unchanged. Swept over 1-3 spatial variables for
% both the nested factor and the plain 'sopvar' factor alongside it. The
% oracle is the action on random inputs (apply_op) before and after
% cleaning.
%

    tol = 1e-8;

    for trial = 1:6
        rng(6500+trial);

        Nw = randi(3);   % number of spatial variables in the nested factor
        Nu = randi(3);   % number of spatial variables in the plain sopvar factor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Part 1: nested 'tensopvar' factor in column 1.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        inner_dims = [randi(2),randi(2)];
        wvars_in = arrayfun(@(k)sprintf('w%d',k),1:Nw,'UniformOutput',false);
        inner_vars = struct('in',{wvars_in},'out',{{'wout'}});
        inner_dom = struct('in',repmat([0,1],Nw,1),'out',[0,1]);
        inner_depmat = struct('in',true(1,Nw),'out',true(1,1));
        Bin1 = rand_tensopvar(inner_dims,inner_vars,inner_dom,inner_depmat,2,[2,2],2,0.8);
        Bin2 = rand_tensopvar(inner_dims,inner_vars,inner_dom,inner_depmat,2,[2,2],2,0.8);
        for ii=1:size(Bin2.ops,1)
            Bin2.ops{ii,1} = 0*Bin2.ops{ii,1};
        end

        uvars_in = arrayfun(@(k)sprintf('u%d',k),1:Nu,'UniformOutput',false);
        svars = struct('in',{uvars_in},'out',{{'sout'}});
        sdom = struct('in',repmat([0,1],Nu,1),'out',[0,1]);
        sdegs = struct('in',2*ones(1,Nu),'out',2);
        msop = [randi(2),randi(2)];
        Sop1 = rand_sopvar(msop,svars,sdom,sdegs,0.8);
        Sop2 = rand_sopvar(msop,svars,sdom,sdegs,0.8);

        Cop = tensopvar();
        Cop.otype = [2,2];
        Cop.ops = cell(2,2);
        Cop.ops{1,1} = Bin1;    Cop.ops{1,2} = Sop1;
        Cop.ops{2,1} = Bin2;    Cop.ops{2,2} = Sop2;
        Cop = initialize(Cop);

        [Cop_clean,is_zero] = clean_op(Cop);
        if is_zero
            error("case_nested failed (Part 1): operator incorrectly flagged entirely zero (trial=%d).",trial)
        elseif size(Cop_clean.ops,1)~=1
            error("case_nested failed (Part 1): expected 1 surviving term, got %d (trial=%d).",size(Cop_clean.ops,1),trial)
        end

        v = {rand_poly([inner_dims(1,2),1],polynomial(Bin1.vars.in(:)),2,20), ...
             rand_poly([msop(2),1],polynomial(Sop1.vars.in(:)),2,20)};
        Y_orig = apply_op(Cop,v);
        Y_clean = apply_op(Cop_clean,v);
        err = cleanpoly(Y_orig-Y_clean,1e-10);
        if max(max(abs(full(err.C))))>tol
            error("case_nested failed (Part 1): cleaned operator does not act identically (trial=%d).",trial)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Part 2: 'intop' factor in column 1.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_int = randi(2);   n_int = randi(2);
        Kop1 = rand_intop([p_int,n_int],{'t1'},[0,1],1,2);
        Kop2 = rand_intop([p_int,n_int],{'t1'},[0,1],1,2);
        for kk=1:numel(Kop2.params)
            Kop2.params{kk} = 0*Kop2.params{kk};
        end

        Sop3 = rand_sopvar(msop,svars,sdom,sdegs,0.8);
        Sop4 = rand_sopvar(msop,svars,sdom,sdegs,0.8);

        Dop = tensopvar();
        Dop.otype = [2,2];
        Dop.ops = cell(2,2);
        Dop.ops{1,1} = Kop1;    Dop.ops{1,2} = Sop3;
        Dop.ops{2,1} = Kop2;    Dop.ops{2,2} = Sop4;
        Dop = initialize(Dop);

        [Dop_clean,is_zero2] = clean_op(Dop);
        if is_zero2
            error("case_nested failed (Part 2): operator incorrectly flagged entirely zero (trial=%d).",trial)
        elseif size(Dop_clean.ops,1)~=1
            error("case_nested failed (Part 2): expected 1 surviving term, got %d (trial=%d).",size(Dop_clean.ops,1),trial)
        end

        v2 = {rand_poly([n_int,1],polynomial({'t1'}),1,20), ...
              rand_poly([msop(2),1],polynomial(Sop3.vars.in(:)),2,20)};
        Y2_orig = apply_op(Dop,v2);
        Y2_clean = apply_op(Dop_clean,v2);
        err2 = cleanpoly(Y2_orig-Y2_clean,1e-10);
        if max(max(abs(full(err2.C))))>tol
            error("case_nested failed (Part 2): cleaned operator does not act identically (trial=%d).",trial)
        end
    end

    disp('case_nested passed.')
end
