% Test_tensopvar_clean_op_nested.m -- tests that @tensopvar/clean_op.m
% recurses into nested 'tensopvar' and 'intop' factors, dropping a term
% whose nested factor is entirely zero and leaving the action of the
% represented operator unchanged. Swept over 1-3 spatial variables for
% both the nested factor and the plain 'sopvar' factor alongside it. The
% oracle is the action on random inputs (apply_op) before and after
% cleaning.
%
% DJ (CC), 10/15/2026: Initial coding;

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
        error("Test_tensopvar_clean_op_nested failed (Part 1): operator incorrectly flagged entirely zero (trial=%d).",trial)
    elseif size(Cop_clean.ops,1)~=1
        error("Test_tensopvar_clean_op_nested failed (Part 1): expected 1 surviving term, got %d (trial=%d).",size(Cop_clean.ops,1),trial)
    end

    v = {rand_poly([inner_dims(1,2),1],polynomial(Bin1.vars.in(:)),2,20), ...
         rand_poly([msop(2),1],polynomial(Sop1.vars.in(:)),2,20)};
    Y_orig = apply_op(Cop,v);
    Y_clean = apply_op(Cop_clean,v);
    err = cleanpoly(Y_orig-Y_clean,1e-10);
    if max(max(abs(full(err.C))))>tol
        error("Test_tensopvar_clean_op_nested failed (Part 1): cleaned operator does not act identically (trial=%d).",trial)
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
        error("Test_tensopvar_clean_op_nested failed (Part 2): operator incorrectly flagged entirely zero (trial=%d).",trial)
    elseif size(Dop_clean.ops,1)~=1
        error("Test_tensopvar_clean_op_nested failed (Part 2): expected 1 surviving term, got %d (trial=%d).",size(Dop_clean.ops,1),trial)
    end

    v2 = {rand_poly([n_int,1],polynomial({'t1'}),1,20), ...
          rand_poly([msop(2),1],polynomial(Sop3.vars.in(:)),2,20)};
    Y2_orig = apply_op(Dop,v2);
    Y2_clean = apply_op(Dop_clean,v2);
    err2 = cleanpoly(Y2_orig-Y2_clean,1e-10);
    if max(max(abs(full(err2.C))))>tol
        error("Test_tensopvar_clean_op_nested failed (Part 2): cleaned operator does not act identically (trial=%d).",trial)
    end
end

disp('Test_tensopvar_clean_op_nested passed.')
