function Test_polyopvar_subs()
% Test_polyopvar_subs.m -- one local function for each case:
%   case_class              of @polyopvar/subs.m (all routes: subs_var,
%                           subs_factor, subs_const), complementing the
%                           Test_polyopvar_subs_* tests with a randomized
%                           sweep and the error and
%   case_const              substitution by a constant (subs_const in
%                           @polyopvar/subs.m), i.e
%   case_const_partial      substitution of only some of the states in a
%                           monomial by constants ('double' or 'polynomial'
%                           values, subs_const in @polyopvar/subs.m)
%   case_fctr               substitution fnew = subs_factor(f,fctr,val),
%                           replacing factor number fctr by val, for a
%                           polynomial functional f (intop coefficient, no
%                           output variables)
%   case_fctr_function      substitution of a factor in a polynomial
%                           function (tensopvar coefficient with output
%                           variables), as opposed to
%                           Test_polyopvar_subs_fctr.m (functionals,
%                           structural check
%   case_multirow           substitution (all routes of @polyopvar/subs.m)
%                           in polynomials with multiple rows of
%                           coefficient operators, built by vertcat of two
%                           polynomials (f = [f1;f2]), possibly
%   case_name               substitution of states by name, and of several
%                           factors at once, in polynomial functionals f =
%                           sum_k Aop_k*Z_k(x) (intop coefficients, as
%                           V(x,u) of piesos_LF),
%   case_nested             substitution of some, or all, of the states in
%                           a monomial acted on by a nested operator (a
%                           factor of the 'tensopvar' is itself a
%                           'tensopvar' acting on several of the factors of
%   case_ode_rows           substitution of the PDE states u in the
%                           monomials Z(u) by the rows of Tx=T*x, for PIEs
%                           with ODE and PDE states (rows of different
%                           spatial dependence: the ODE rows depend on no
%                           variable
%   case_var                exercise the subs_var (variable-name
%                           substitution) path of @polyopvar/subs.m, as
%                           distinct from Test_polyopvar_subs_fctr.m, which
%                           only exercises subs_factor
%   case_var_function       substitution of a state variable by name
%                           (subs_var in @polyopvar/subs.m) in (1)
%                           polynomial functions f(a,b)(s) = Kop*(a o b),
%                           tensopvar
%   case_vector_states      substitution (all routes of @polyopvar/subs.m)
%                           in polynomials of vector-valued states: (1)
%                           subs_var by name, by a polynomial function
%                           val(y) of a vector
%
% DJ (CC), 10/23/2026: Merge Test_polyopvar_subs_class,
%   Test_polyopvar_subs_const, Test_polyopvar_subs_const_partial,
%   Test_polyopvar_subs_fctr, Test_polyopvar_subs_fctr_function,
%   Test_polyopvar_subs_multirow, Test_polyopvar_subs_name,
%   Test_polyopvar_subs_nested, Test_polyopvar_subs_ode_rows,
%   Test_polyopvar_subs_var, Test_polyopvar_subs_var_function,
%   Test_polyopvar_subs_vector_states into this test, as cases;
% DJ (CC), 10/28/2026: In case_class, substitute constants in functionals of
%   degree d>1 (compared with the tensor product evaluation) instead of
%   expecting an error;

cases = {@case_class, @case_const, @case_const_partial, @case_fctr, @case_fctr_function, @case_multirow, @case_name, @case_nested, @case_ode_rows, @case_var, @case_var_function, @case_vector_states};
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
    error('Test_polyopvar_subs: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_polyopvar_subs: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_class()
% case_class -- tests of @polyopvar/subs.m (all routes:
% subs_var, subs_factor, subs_const), complementing the
% Test_polyopvar_subs_* tests with a randomized sweep and the error and
% edge cases.
%
% Expected values are computed independently of subs: the original
% polynomial is evaluated (eval_polyopvar, or apply_op of the coefficient
% operators on the factors) at the constant or at the polynomial val(y)
% substituted, and compared with the result of the substitution evaluated
% at the remaining states.
% Parts:
%   (1) dispatch and input validation;
%   (2) renaming of states (subs_var, var and val cellstr);
%   (3) errors and edge cases of subs_var, subs_factor, subs_const;
%   (4) result is a valid polynomial with the output space of f
%       unchanged;
%   (5) sweep over 1-3 spatial variables, pointwise/Kronecker output,
%       vector-valued states, repeated factors: the first state and all
%       states, in random order, by name by a constant ('double',
%       'polynomial', 'dpvar') and by a function of a new state (sorted
%       between the states, so that the states are renumbered), the
%       first and last factor by position, and (5b) polynomials with
%       several monomials and several rows;
%   (6) functionals (intop coefficients).
% Cost: about 1.5 minutes.
%
% DJ, 10/15/2026: Initial coding;

    tol = 1e-8;
    ncases = 0;
    pvar s1 s2 s3
    svar_all = {s1,s2,s3};
    allnames = {'a','b','c'};
    sz1 = [1;1;1];

    %% (1) dispatch and input validation
    rng(1);
    f = make_f(1,[1,2,3],sz1,'pw');
    xv = struct('a',rand_pol(svar_all(1),1),'b',rand_pol(svar_all(1),1),'c',rand_pol(svar_all(1),1));
    yy = polyopvar('y',{'s1'},[0,1]);
    expect_error(@() subs(5,{'a'},yy),"must be of type 'polyopvar'");        % dispatches to polyopvar/subs through val
    expect_error(@() subs(f,'a',1),"must be expressed as 'cellstr'");
    expect_error(@() subs(f,'a',yy),"must be expressed as 'polyopvar' or numeric");
    expect_error(@() subs(f,{'a'},'y'),'');
    ncases = ncases+4;
    % f that has not been validated is initialized first
    f2 = f;     f2.chkval = false;
    c = randn;
    assert_close(eval_new(subs(f2,{'a'},c),xv,sz1,allnames),eval_polyopvar(f,{polynomial(c);xv.b;xv.c}),tol,'(1) chkval=false');
    % ... also if its fields are as set by the user, not yet validated by initialize
    rng(1);      g = make_f(1,[1,2,3],sz1,'pw',true);       % same operator as f
    assert_close(eval_new(subs(g,{'a'},c),xv,sz1,allnames),eval_polyopvar(f,{polynomial(c);xv.b;xv.c}),tol,'(1) not initialized');
    ncases = ncases+2;

    %% (2) renaming of states
    rng(2);
    f = make_f(1,[1,2,3],sz1,'pw') + make_f(1,[1,3,2],sz1,'pw');
    Yref = eval_polyopvar(f,{xv.a;xv.b;xv.c});
    % a --> z moves the state after c; names given as column or row, states not in f are ignored
    fnew = subs(f,{'a';'q'},{'z';'w'});
    assert(isequal(fnew.varname,{'b';'c';'z'}),'(2) rename: states not sorted');
    assert_close(eval_polyopvar(fnew,{xv.b;xv.c;xv.a}),Yref,tol,'(2) a --> z');
    fnew = subs(f,{'a','b'},{'b','a'});
    assert_close(eval_polyopvar(fnew,{xv.b;xv.a;xv.c}),Yref,tol,'(2) swap a,b');
    % names listed in another order than the states
    fnew = subs(f,{'c','a'},{'x','y'});
    assert(isequal(fnew.varname,{'b';'x';'y'}),'(2) rename, var not in the order of the states');
    assert_close(eval_polyopvar(fnew,{xv.b;xv.c;xv.a}),Yref,tol,'(2) c,a --> x,y');
    fnew = subs(f,{'q'},{'w'});
    assert(isequal(fnew.varname,f.varname) && isequal(fnew.Z,f.Z),'(2) rename of a state that is not in f');
    % c --> b merges the states, and the monomials a*b*b that now coincide
    fnew = subs(f,{'c'},{'b'});
    assert(isequal(fnew.varname,{'a';'b'}) && size(fnew.Z,1)==1,'(2) c --> b: states or monomials not merged');
    assert_close(eval_polyopvar(fnew,{xv.a;xv.b}),eval_polyopvar(f,{xv.a;xv.b;xv.b}),tol,'(2) c --> b');
    % states of different size cannot be merged
    fv = make_f(1,[1,2],[2;1;1],'pw');
    expect_error(@() subs(fv,{'a'},{'b'}),'same function space');
    ncases = ncases+8;

    %% (3) errors and edge cases
    rng(3);
    svars = {'s1'};     dom = struct('in',[0,1],'out',[0,1]);
    f = make_f(1,[1,2,3],sz1,'pw');
    val1 = make_val(svars,dom,'y',1,1,1);
    val2 = [val1;val1];
    % subs_var
    expect_error(@() subs(f,yy+yy,val1),'individial states');
    expect_error(@() subs(f,[yy;yy],val2),'unique');
    expect_error(@() subs(f,{'a','b'},val1),'Number of rows');
    expect_error(@() subs(f,{'a'},make_val(svars,dom,'y',2,1,1)),'Output space');
    expect_error(@() subs(f,{'a'},make_val({'s2'},struct('in',[0,1],'out',[0,1]),'y',1,1,1)),'Output space');
    % names that are not states of f are ignored, together with their rows of val: nothing is substituted
    fnew = subs(f,{'q'},val1);
    assert(isequal(fnew.Z,f.Z) && isequal(fnew.varname,f.varname),'(3) unknown state should leave f unchanged');
    vq = make_val(svars,dom,'y',1,1,1);     va = make_val(svars,dom,'y',1,1,2);
    yv = struct('b',xv.b,'c',xv.c,'y',rand_pol(svar_all(1),1));
    Yref = eval_new(subs(f,{'a'},va),yv,[sz1;1],[allnames,{'y'}]);
    assert_close(eval_new(subs(f,{'q','a'},[vq;va]),yv,[sz1;1],[allnames,{'y'}]),Yref,tol,'(3) unknown and known state');
    expect_error(@() subs(f,{'a','q','r'},val2),'Number of rows');ncases = ncases+6;
    % subs_factor
    expect_error(@() subs(f,1.5,val1),'unique positive integers');
    expect_error(@() subs(f,0,val1),'unique positive integers');
    expect_error(@() subs(f,-1,val1),'unique positive integers');
    expect_error(@() subs(f,[1,1],val2),'unique positive integers');
    expect_error(@() subs(f,4,val1),'cannot exceed');
    expect_error(@() subs(f,1,'y'),'non-''polyopvar''');
    expect_error(@() subs(f,[1,2],val1),'Number of rows');
    expect_error(@() subs(f,1,make_val(svars,dom,'y',2,1,1)),'Input function space');       % size of the state
    expect_error(@() subs(f,1,make_val({'s2'},dom,'y',1,1,1)),'Input function space');            % variable of f
    expect_error(@() subs(make_f(2,[1,2],sz1,'pw'),1,make_val({'s1'},dom,'y',1,1,1)),'Input function space');    % fewer variables than the state
    expect_error(@() subs(f+make_f(1,[1,2],sz1,'pw'),1,val1),'multiple monomials');
    fk = subs(f,{'a'},randn);               % Z=[0,2,3]: the factor 1 is a multiplier from R
    expect_error(@() subs(fk,1,val1),'not multipliers');
    ncases = ncases+9;
    % subs_const
    % dd declared with an output: the anonymous functions below cannot capture a
    % variable that dpvar creates in the workspace
    dd = dpvar('dd');
    expect_error(@() subs(f,{'a'},dd),'decision variables');
    expect_error(@() subs(f,{'a','b'},[1,2]),'column vector');
    expect_error(@() subs(f,yy+yy,1),'individial states');
    expect_error(@() subs(f,{'q'},1),'must be states');
    expect_error(@() subs(f,{'a','a'},[1;2]),'unique');
    expect_error(@() subs(f,{'a'},[1;2]),'Number of rows of constant');
    expect_error(@() subs(f,{'a'},s2),'spatial variables');
    expect_error(@() subs(f,{'a','b','c'},[1;2;s2]),'spatial variables');
    expect_error(@() subs(f+make_f(1,[1,2],sz1,'pw'),1,1),'single monomial');
    expect_error(@() subs(f,[1,4],[1;2]),'cannot exceed|not exceeding');
    expect_error(@() subs(f,[1,1],[1;2]),'unique');
    expect_error(@() subs(fk,1,1),'not multipliers');
    expect_error(@() subs(f,[1,2],1),'Number of rows of constant');
    expect_error(@() subs(make_f(1,[1,2],sz1,'kron'),{'a','b'},[1;2]),'tensor product');
    % a state that is not in any monomial leaves f unchanged
    fnew = subs(make_f(1,[1,2],sz1,'pw'),{'c'},1);
    assert(isequal(size(fnew.Z),[1,2]),'(3) unused state should leave f unchanged');
    % dpvar constants without decision variables are accepted
    cc = randn;
    assert_close(eval_new(subs(f,{'a'},dpvar(cc)),xv,sz1,allnames),eval_polyopvar(f,{polynomial(cc);xv.b;xv.c}),tol,'(3) dpvar constant');
    ncases = ncases+16;

    %% (4) the result is a valid polynomial, and keeps the output space of f
    rng(4);
    for N = 1:2
        [svars,dom] = spatial_space(N);
        for kind = {'pw','kron'}
            f = make_f(N,[1,2,3],sz1,kind{1});
            val = make_val(svars,dom,'y',1,1,2);
            cases = {subs(f,{'a'},val), subs(f,2,val), subs(f,{'a'},randn), subs(f,2,randn), subs(f,{'a'},{'z'})};
            for k=1:numel(cases)
                fnew = cases{k};
                g = fnew;   g.chkval = false;   g = initialize(g);
                assert(isequal(fnew.vars.out,f.vars.out) && isequal(fnew.dom.out,f.dom.out) && isequal(fnew.dims.out,f.dims.out) ...
                    && isequal(fnew.depmat.out,f.depmat.out) && fnew.chkval && isequal(g.Z,fnew.Z),sprintf('(4) N=%d,%s,%d',N,kind{1},k));
                ncases = ncases+1;
            end
        end
    end

    %% (5) randomized sweep: functions
    Zlist = {[1,2],[2,3],[1,2,3],[1,2,1],[3,1,1]};
    for N = 1:3
        [svars,dom] = spatial_space(N);     svs = svar_all(1:N);
        for kind = {'pw','kron'}
            for iz = 1:numel(Zlist)
                Z = Zlist{iz};
                % cost: skip 3 factors in 2D with Kronecker output, and 3D except 2 factors with pointwise output
                if N==2 && ~(iz==1 || iz==3 && strcmp(kind{1},'pw')) || N==3 && (numel(Z)>2 || strcmp(kind{1},'kron')),     continue,   end
                szlist = {sz1};     if iz==1 || iz==3,      szlist = {sz1,[2;1;2]};     end     % vector-valued states, cost
                for szc = szlist
                    sz = szc{1};
                    if N>1 && sz(1)>1,      continue,   end
                    seed = 10000*N + 1000*numel(kind{1}) + 100*iz + 10*sz(1);
                    rng(seed);
                    label = sprintf('N=%d,%s,Z=[%s],sz=[%s]',N,kind{1},num2str(Z),num2str(sz'));
                    f = make_f(N,Z,sz,kind{1});
                    x = struct('a',rand_pol(svs,sz(1)),'b',rand_pol(svs,sz(2)),'c',rand_pol(svs,sz(3)));
                    names = allnames(unique(Z));
                    % the first state, and all of them, each by a constant of a different type
                    nsub = numel(names);
                    masks = unique([1,2^nsub-1]);       if N==3,    masks = 1;      end     % cost
                    for mask = masks
                        sub = names(logical(bitget(mask,1:numel(names))));
                        lbl = sprintf('%s,sub={%s}',label,strjoin(sub,','));
                        is_all = numel(sub)==numel(names);
                        forms = {'double','polynomial','dpvar'};
                        for form = forms(1+mod(mask+iz,3))
                            cvals = cellfun(@(nm) rand_const(svs,sz(strcmp(allnames,nm)),form{1}),sub,'UniformOutput',false);
                            if is_all && strcmp(kind{1},'kron')
                                expect_error(@() subs(f,sub,vertcat(cvals{:})),'tensor product');
                                ncases = ncases+1;
                            else
                                ncases = ncases+check_const(f,sub,cvals,x,sz,tol,[lbl,',',form{1}]);
                            end
                        end
                        if N<3,     ncases = ncases+check_poly(f,sub,x,sz,N,1+mod(mask,2),tol,lbl);     end
                    end
                    % each factor by position, by a function of a new state, and by a constant
                    poslist = unique([1,numel(Z)]);     if N>1,     poslist = 1;    end         % cost
                    if N==3,    poslist = [];   end
                    for pos = poslist
                        ncases = ncases+check_pos(f,Z,pos,x,sz,N,tol,[label,sprintf(',pos=%d',pos)]);
                    end
                end
            end
        end
    end

    %% (5b) several monomials and several rows
    for N = 1:2
        [svars,dom] = spatial_space(N);     svs = svar_all(1:N);
        rng(20000+N);
        sz = sz1;
        x = struct('a',rand_pol(svs,1),'b',rand_pol(svs,1),'c',rand_pol(svs,1));
        fs = {  make_f(N,[1,2],sz,'pw') + make_f(N,1,sz,'pw') + make_f(N,[2,3],sz,'pw') + make_f(N,[1,2,3],sz,'pw'),   'sum of monomials'
                [make_f(N,[1,2],sz,'pw'); make_f(N,[1,3],sz,'pw')],                                                  'rows'
                [make_f(N,[1,2],sz,'pw') + make_f(N,3,sz,'pw'); make_f(N,[1,1,3],sz,'pw')],                         'rows and monomials' };
        for k=1:size(fs,1)
            f = fs{k,1};
            for sub = {{'a'},{'b'},{'c'},{'a','c'},{'a','b','c'}}
                lbl = sprintf('(5b) N=%d,%s,sub={%s}',N,fs{k,2},strjoin(sub{1},','));
                for form = {'double','polynomial'}
                    cvals = cellfun(@(nm) rand_const(svs,1,form{1}),sub{1},'UniformOutput',false);
                    ncases = ncases+check_const(f,sub{1},cvals,x,sz,tol,[lbl,',',form{1}]);
                end
                ncases = ncases+check_poly(f,sub{1},x,sz,N,2,tol,lbl);
            end
        end
    end

    %% (6) functionals
    for N = 1:2
        [svars,dom] = spatial_space(N);     svs = svar_all(1:N);
        Zlist_V = {1,[1;2],[1,2],[1,2;3,0],[1,2,0;3,1,2]};
        for iz = 1:numel(Zlist_V)
            Z = Zlist_V{iz};
            if size(Z,2)==3 && N==2,    continue,   end
            rng(30000+10*N+iz);
            V = make_V(svars,dom,Z,0);
            lbl = sprintf('(6) N=%d,Z=%d',N,iz);
            % all states by functions of y
            ncases = ncases+check_poly_V(V,N,tol,lbl);
            % all states by constants, also for monomials of degree d>1, on which the functional acts
            % by the tensor product of the d values: compared with eval_polyopvar for degree one, and with
            % the functional applied to the tensor product built here for the other monomials
            cvals = {rand_const(svs,1,'polynomial');rand_const(svs,1,'polynomial');rand_const(svs,1,'polynomial')};
            fnew = subs(V,allnames,vertcat(cvals{:}));
            if size(Z,2)==1
                Yref = eval_polyopvar(V,cvals);
            else
                Yref = eval_V_tensor(V,cvals,svars);
            end
            assert_close(eval_new(fnew,struct(),sz1,allnames),Yref,tol,[lbl,' (all states by constants)']);
            ncases = ncases+1;
        end
    end

    fprintf('case_class passed (%d checks).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_const()
% case_const -- substitution by a constant (subs_const
% in @polyopvar/subs.m), i.e. 'double' or 'polynomial' values for
% states/factors, by name (cellstr) and by factor index. If only some of
% the factors in a monomial are substituted, the others remain and the
% constant becomes a known factor (part (6)).
%
% Oracle: fnew is evaluated at the remaining states (eval_polyopvar) and
% compared with the original f evaluated at the constants and the same
% remaining states.
%
%
% NOTES:
% - Evaluation uses eval_polyopvar (Testfolder/polyopvar), not apply_poly:
%   apply_poly evaluates by subs, which this test checks.

    pvar s1 s2
    svar_all = {s1,s2};
    tol = 1e-8;
    ncases = 0;

    for Ntot = 1:2
        svs = svar_all(1:Ntot);
        for cform = {'double','polynomial'}
            rng(300+10*Ntot+numel(cform{1}));
            label = string(sprintf('Ntot=%d,%s',Ntot,cform{1}));
            ca = rand_val(svs,cform{1});       cb = rand_val(svs,cform{1});
            b = rand_val(svs,'polynomial');

            %% (1) All factors of a two-factor function monomial, by name
            f = build_f(Ntot,[1,2],'ab');
            fnew = subs(f,{'a';'b'},[ca;cb]);
            assert_deg0(fnew,label+" (1)");
            assert_close_const(eval_new_const(fnew,{}),eval_polyopvar(f,{ca;cb}),tol,label+" (1)");
            % ... the order of the names must be respected
            fnew = subs(f,{'b';'a'},[cb;ca]);
            assert_close_const(eval_new_const(fnew,{}),eval_polyopvar(f,{ca;cb}),tol,label+" (1b)");
            ncases = ncases+2;

            %% (2) Same state twice, a o a, by name
            f = build_f(Ntot,[1,2],'aa');
            fnew = subs(f,{'a'},ca);
            assert_close_const(eval_new_const(fnew,{}),eval_polyopvar(f,{ca;b}),tol,label+" (2)");
            ncases = ncases+1;

            %% (3) Single-factor monomials, function and functional
            for kind = {'function','functional'}
                f = build_f(Ntot,[1,2],'a',kind{1});
                fnew = subs(f,{'a'},ca);
                assert_close_const(eval_new_const(fnew,{}),eval_polyopvar(f,{ca;b}),tol,label+" (3) "+kind{1});
                ncases = ncases+1;
            end

            % ... also with a Kronecker-type (single-factor) output, for which the
            % evaluation renames the output variables s -> s_1
            f = build_f(Ntot,[2,2],'a');
            fnew = subs(f,{'a'},ca);
            Yref = eval_polyopvar(f,{ca;b});
            Yref.varname = regexprep(Yref.varname,"_1$","");
            assert_close_const(eval_new_const(fnew,{}),Yref,tol,label+" (3b)");
            ncases = ncases+1;
            %% (4) Several monomials: those with substituted states change, the rest is kept
            f = build_f(Ntot,[1,2],'a+b');
            fnew = subs(f,{'a'},ca);
            assert_close_const(eval_new_const(fnew,{'b',b}),eval_polyopvar(f,{ca;b}),tol,label+" (4)");
            fnew = subs(f,{'a';'b'},[ca;cb]);
            assert_close_const(eval_new_const(fnew,{}),eval_polyopvar(f,{ca;cb}),tol,label+" (4b)");
            ncases = ncases+2;

            %% (5) Factors by position
            f = build_f(Ntot,[1,2],'ab');
            fnew = subs(f,[1,2],[ca;cb]);
            assert_close_const(eval_new_const(fnew,{}),eval_polyopvar(f,{ca;cb}),tol,label+" (5)");
            fnew = subs(f,[2,1],[cb;ca]);
            assert_close_const(eval_new_const(fnew,{}),eval_polyopvar(f,{ca;cb}),tol,label+" (5b)");
            ncases = ncases+2;

            %% (6) Unsupported cases must error
            f = build_f(Ntot,[1,2],'ab');
            % Partial monomials are supported: the substituted factor remains as a known factor
            % (see also Test_polyopvar_subs/case_const_partial)
            aa = rand_val(svs,'polynomial');
            for kind = {'pw','kron'}
                if strcmp(kind{1},'pw'),    fp = f;
                else,                       fp = build_f(Ntot,[2,2],'ab');      end
                fnew = subs(fp,{'a'},ca);
                assert_close_const(eval_new_const(fnew,{'b',b}),eval_polyopvar(fp,{ca;b}),tol,label+" (6) partial a, "+kind{1});
                fnew = subs(fp,1,ca);
                assert_close_const(eval_new_const(fnew,{'b',b}),eval_polyopvar(fp,{ca;b}),tol,label+" (6) partial factor 1, "+kind{1});
                fnew = subs(fp,{'b'},cb);
                assert_close_const(eval_new_const(fnew,{'a',aa}),eval_polyopvar(fp,{aa;cb}),tol,label+" (6) partial b, "+kind{1});
                fnew = subs(fp,2,cb);
                assert_close_const(eval_new_const(fnew,{'a',aa}),eval_polyopvar(fp,{aa;cb}),tol,label+" (6) partial factor 2, "+kind{1});
                ncases = ncases+4;
            end
            expect_error_const(@() subs(f,{'a';'b'},ca),'Number of rows');            % wrong number of rows
            expect_error_const(@() subs(f,{'c'},ca),'states');                        % unknown state
            expect_error_const(@() subs(f,{'a';'a'},[ca;ca]),'unique');               % repeated state
            expect_error_const(@() subs(f,{'a';'b'},[ca,ca]),'column');               % not a column
            if Ntot==1
                expect_error_const(@() subs(f,{'a';'b'},[s2;s2]),'spatial variables');   % wrong variable
            end
            fk = build_f(Ntot,[2,2],'ab');
            expect_error_const(@() subs(fk,{'a';'b'},[ca;cb]),'tensor product');      % Kronecker output
            ncases = ncases+1;
        end
    end

    fprintf('case_const passed (%d substitution cases).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_const_partial()
% case_const_partial -- substitution of only some of the
% states in a monomial by constants ('double' or 'polynomial' values,
% subs_const in @polyopvar/subs.m). The factors replaced remain in the
% operator as multipliers from R (zeros in Z, see polyopvar) with the
% constant as value, so that e.g. x*y with x=c gives Z=[0,y]. Parts:
%   (A) two factors, pointwise and Kronecker output, by name and by
%       position;
%   (B) three factors, one or two of them replaced;
%   (C) a state repeated in the monomial, a o b o a;
%   (D) states depending on different spatial variables (and on none);
%   (E) vector-valued states;
%   (F) several monomials, with unaffected, partial and full ones, and
%       (F2) combination of consecutive multipliers from R after partial
%       substitution (pointwise output);
%   (G) functionals (intop), and further operations on the result;
%   (I) subs_var, by name, with a polyopvar for states depending on some
%       of the variables only;
%   (H) errors.
%
% Oracle: the substituted polynomial is evaluated at the remaining states
% (eval_polyopvar) and compared with the original polynomial evaluated
% at the constants and the same remaining states.
%
%
% NOTES:
% - Evaluation uses eval_polyopvar (Testfolder/polyopvar), not
%   apply_poly: apply_poly evaluates by subs, which this test checks.
% - Expected Z for a pointwise ('pw') output has consecutive multipliers
%   from R combined (compress_zeros); for a Kronecker output they are
%   left as is.

    pvar s1 s2
    svar_all = {s1,s2};
    tol = 1e-8;
    ncases = 0;
    allnames = {'a','b','c'};

    %% (A)-(C) states depending on all variables, scalar
    for N = 1:2
        svs = svar_all(1:N);
        dep = true(3,N);        sz = [1;1;1];
        for kind = {'pw','kron'}
            for form = {'double','polynomial'}
                rng(100*N+10*numel(kind{1})+numel(form{1}));
                lbl = sprintf('N=%d,%s,%s',N,kind{1},form{1});
                xv = struct('a',rand_const_const_partial(svs,1,'polynomial'),'b',rand_const_const_partial(svs,1,'polynomial'),'c',rand_const_const_partial(svs,1,'polynomial'));
                ca = rand_const_const_partial(svs,1,form{1});     cb = rand_const_const_partial(svs,1,form{1});     cc = rand_const_const_partial(svs,1,form{1});
                % (A) a o b
                f = make_f_const_partial(N,[1,2],dep,sz,kind{1});
                ncases = ncases+check_partial(f,{'a'},{ca},xv,sz,[0,2],tol,lbl+" (A) a");
                ncases = ncases+check_partial(f,{'b'},{cb},xv,sz,[1,0],tol,lbl+" (A) b");
                % (B) a o b o c
                f = make_f_const_partial(N,[1,2,3],dep,sz,kind{1});
                ncases = ncases+check_partial(f,{'a'},{ca},xv,sz,[0,2,3],tol,lbl+" (B) a");
                ncases = ncases+check_partial(f,{'b'},{cb},xv,sz,[1,0,3],tol,lbl+" (B) b");
                ncases = ncases+check_partial(f,{'c'},{cc},xv,sz,[1,2,0],tol,lbl+" (B) c");
                ncases = ncases+check_partial(f,{'a';'c'},{ca;cc},xv,sz,[0,2,0],tol,lbl+" (B) a,c");
                ncases = ncases+check_partial(f,{'c';'a'},{cc;ca},xv,sz,[0,2,0],tol,lbl+" (B) c,a");
                % (C) a o b o a: the state in two factors
                f = make_f_const_partial(N,[1,2,1],dep,sz,kind{1});
                ncases = ncases+check_partial(f,{'a'},{ca},xv,sz,[0,2,0],tol,lbl+" (C) a");
                ncases = ncases+check_partial(f,{'b'},{cb},xv,sz,[1,0,1],tol,lbl+" (C) b");
                % one of the two factors a, by position
                % (evaluation by state gives all factors of a state the same value: apply the operator to the factors)
                Yref = apply_op(f.C,{{xv.a,xv.b,ca}});
                fnew = subs(f,3,ca);
                assert_close_const_partial(eval_new_const_partial(fnew,xv,sz,allnames),Yref,tol,lbl+" (C) factor 3");
                if ~isequal(pad(full(fnew.Z),3),[1,2,0]),   error("case_const_partial failed (%s): Z = [%s].",lbl,num2str(full(fnew.Z))),   end
                ncases = ncases+1;
            end
        end
    end

    %% (D) states depending on different variables: a(s1,s2), b(s1), c(no variables)
    dep = [true,true; true,false; false,false];     sz = [1;1;1];
    N = 2;      svs = svar_all;
    for kind = {'pw','kron'}
        for form = {'double','polynomial'}
            rng(500+10*numel(kind{1})+numel(form{1}));
            lbl = sprintf('(D) %s,%s',kind{1},form{1});
            xv = struct('a',rand_const_const_partial(svs,1,'polynomial'),'b',rand_const_const_partial({s1},1,'polynomial'),'c',rand_const_const_partial({},1,'polynomial'));
            ca = rand_const_const_partial(svs,1,form{1});     cb = rand_const_const_partial({s1},1,form{1});    cc = rand_const_const_partial({},1,form{1});
            f = make_f_const_partial(N,[1,2,3],dep,sz,kind{1});
            ncases = ncases+check_partial(f,{'a'},{ca},xv,sz,[0,2,3],tol,lbl+" a");
            ncases = ncases+check_partial(f,{'b'},{cb},xv,sz,[1,0,3],tol,lbl+" b");
            ncases = ncases+check_partial(f,{'c'},{cc},xv,sz,[1,2,0],tol,lbl+" c");
            ncases = ncases+check_partial(f,{'b';'c'},{cb;cc},xv,sz,[1,0,0],tol,lbl+" b,c");
            ncases = ncases+check_partial(f,{'a';'b'},{ca;cb},xv,sz,[0,0,3],tol,lbl+" a,b");
        end
    end

    %% (E) vector-valued states, sizes [2;3;1]
    sz = [2;3;1];
    for N = 1:2
        svs = svar_all(1:N);    dep = true(3,N);
        for form = {'double','polynomial'}
            rng(700+N+numel(form{1}));
            lbl = sprintf('(E) N=%d,%s',N,form{1});
            xv = struct('a',rand_const_const_partial(svs,2,'polynomial'),'b',rand_const_const_partial(svs,3,'polynomial'),'c',rand_const_const_partial(svs,1,'polynomial'));
            ca = rand_const_const_partial(svs,2,form{1});     cb = rand_const_const_partial(svs,3,form{1});     cc = rand_const_const_partial(svs,1,form{1});
            f = make_f_const_partial(N,[1,2,3],dep,sz,'pw');
            ncases = ncases+check_partial(f,{'a'},{ca},xv,sz,[0,2,3],tol,lbl+" a");
            ncases = ncases+check_partial(f,{'b'},{cb},xv,sz,[1,0,3],tol,lbl+" b");
            ncases = ncases+check_partial(f,{'c';'a'},{cc;ca},xv,sz,[0,2,0],tol,lbl+" c,a");
        end
    end

    %% (F) several monomials: a o b (partial), a (full), b (unaffected), plus a o b o c
    N = 1;      svs = svar_all(1);      dep = true(3,N);    sz = [1;1;1];
    rng(900);
    xv = struct('a',rand_const_const_partial(svs,1,'polynomial'),'b',rand_const_const_partial(svs,1,'polynomial'),'c',rand_const_const_partial(svs,1,'polynomial'));
    fm = make_f_const_partial(N,[1,2],dep,sz,'pw') + make_f_const_partial(N,1,dep,sz,'pw') + make_f_const_partial(N,2,dep,sz,'pw') + make_f_const_partial(N,[1,2,3],dep,sz,'pw');
    for form = {'double','polynomial'}
        ca = rand_const_const_partial(svs,1,form{1});     cb = rand_const_const_partial(svs,1,form{1});
        lbl = "(F) "+string(form{1});
        ncases = ncases+check_partial(fm,{'a'},{ca},xv,sz,[],tol,lbl+" a");
        ncases = ncases+check_partial(fm,{'b'},{cb},xv,sz,[],tol,lbl+" b");
        ncases = ncases+check_partial(fm,{'a';'b'},{ca;cb},xv,sz,[],tol,lbl+" a,b");
        % a state that appears in no monomial leaves the polynomial unchanged
        cz = rand_const_const_partial(svs,1,form{1});
        fnew = subs(fm,{'c'},cz);
        xs = xv;    xs.c = cz;
        assert_close_const_partial(eval_new_const_partial(fnew,rmfield(xv,'c'),sz,allnames),eval_polyopvar(fm,{xs.a;xs.b;xs.c}),tol,lbl+" c");
        ncases = ncases+1;
    end

    %% (F2) After partial substitution, consecutive multipliers from R are combined (pointwise output)
    % and the result is cleaned: monomials that coincide are added, states that no longer appear removed.
    % a o b o c and b o a o c, with a,b replaced: both Z=[0,0,c] --> [0,c] and then one monomial
    N = 1;      svs = svar_all(1);      dep = true(3,N);    sz = [1;1;1];
    for kind = {'pw','kron'}
        rng(950+numel(kind{1}));
        lbl = "(F2) "+string(kind{1});
        xv = struct('a',rand_const_const_partial(svs,1,'polynomial'),'b',rand_const_const_partial(svs,1,'polynomial'),'c',rand_const_const_partial(svs,1,'polynomial'));
        ca = rand_const_const_partial(svs,1,'double');        cb = rand_const_const_partial(svs,1,'polynomial');
        fpair = make_f_const_partial(N,[1,2,3],dep,sz,kind{1}) + make_f_const_partial(N,[2,1,3],dep,sz,kind{1});
        xs = xv;    xs.a = ca;      xs.b = cb;
        Yref = eval_polyopvar(fpair,{xs.a;xs.b;xs.c});
        fnew = subs(fpair,{'a';'b'},[ca;cb]);
        assert_close_const_partial(eval_new_const_partial(fnew,rmfield(xv,{'a','b'}),sz,allnames),Yref,tol,lbl);
        Zexp = [0,0,1];     if strcmp(kind{1},'pw'),    Zexp = [0,1];   end
        n = max(size(fnew.Z,2),numel(Zexp));
        if size(fnew.Z,1)~=1 || ~isequal(pad(full(fnew.Z),n),pad(Zexp,n)) || any(ismember({'a','b'},fnew.varname))
            error("case_const_partial failed (%s): Z = [%s] in states {%s}, expected the single monomial [%s] in state c.",lbl,num2str(full(fnew.Z)),strjoin(fnew.varname,','),num2str(Zexp))
        end
        ncases = ncases+2;
    end

    %% (G) functionals, and further operations on the result
    polyopvar a(s1) b(s1)
    V = innerprod(a,b);         % int a*b ds, Z=[1,2]
    rng(1000);
    xa = rand_const_const_partial({s1},1,'polynomial');       xb = rand_const_const_partial({s1},1,'polynomial');
    vals = struct('a',xa,'b',xb);
    for form = {'double','polynomial'}
        cc = rand_const_const_partial({s1},1,form{1});
        tests = {   subs(V,{'b'},cc),       int(xa*cc,s1,0,1)
                    subs(V,2,cc),           int(xa*cc,s1,0,1)
                    subs(V,{'a'},cc),       int(cc*xb,s1,0,1)       };
        for k=1:size(tests,1)
            assert_close_const_partial(eval_V(tests{k,1},vals),tests{k,2},tol,sprintf('(G) functional %d, %s',k,form{1}));
            ncases = ncases+1;
        end
    end
    % further operations: product and sum with other polynomials of the same output space
    N = 1;      svs = svar_all(1);      dep = true(3,N);    sz = [1;1;1];
    rng(1100);
    xv = struct('a',rand_const_const_partial(svs,1,'polynomial'),'b',rand_const_const_partial(svs,1,'polynomial'),'c',rand_const_const_partial(svs,1,'polynomial'));
    ca = rand_const_const_partial(svs,1,'double');
    f = make_f_const_partial(N,[1,2],dep,sz,'pw');        fz = make_f_const_partial(N,3,dep,sz,'pw');
    g = subs(f,{'a'},ca);
    xs = xv;    xs.a = ca;
    Yf = eval_polyopvar(f,{xs.a;xs.b;xs.c});       Yz = eval_polyopvar(fz,{xs.a;xs.b;xs.c});
    assert_close_const_partial(eval_new_const_partial(g.*fz,rmfield(xv,'a'),sz,allnames),Yf.*Yz,tol,'(G) times');
    assert_close_const_partial(eval_new_const_partial(fz.*g,rmfield(xv,'a'),sz,allnames),Yz.*Yf,tol,'(G) times, reversed');
    assert_close_const_partial(eval_new_const_partial(g+fz,rmfield(xv,'a'),sz,allnames),Yf+Yz,tol,'(G) plus');
    ncases = ncases+3;

    %% (I) subs_var, by name, with a polyopvar for states depending on only some of the variables:
    % a(s1,s2), b(s1), c(no variables). The output variables of the substitute are those of the
    % state, a subset of those of f, so that the function spaces must be compared by variable name
    N = 2;      svs = svar_all;     dep = [true,true; true,false; false,false];     sz = [1;1;1];
    polyopvar w(s1) v(s1)
    for kind = {'pw','kron'}
        rng(1200+numel(kind{1}));
        lbl = "(I) "+string(kind{1});
        xv = struct('a',rand_const_const_partial(svs,1,'polynomial'),'b',rand_const_const_partial({s1},1,'polynomial'),'c',rand_const_const_partial({},1,'polynomial'));
        ca = rand_const_const_partial(svs,1,'polynomial');    cb = rand_const_const_partial({s1},1,'polynomial');   cc = rand_const_const_partial({},1,'polynomial');
        ga = known_fun(ca,{'s1','s2'});         gb = known_fun(cb,{'s1'});              gc = known_fun(cc,cell(1,0));
        f = make_f_const_partial(N,[1,2,3],dep,sz,kind{1});
        % known functions
        ncases = ncases+check_known(f,{'a'},{ga},{ca},xv,sz,tol,lbl+" a");
        ncases = ncases+check_known(f,{'b'},{gb},{cb},xv,sz,tol,lbl+" b");
        ncases = ncases+check_known(f,{'c'},{gc},{cc},xv,sz,tol,lbl+" c");
        ncases = ncases+check_known(f,{'b';'c'},{gb;gc},{cb;cc},xv,sz,tol,lbl+" b,c");
        ncases = ncases+check_known(f,{'a';'b'},{ga;gb},{ca;cb},xv,sz,tol,lbl+" a,b");
        % functions of new states w,v on the same variables as b
        xw = rand_const_const_partial({s1},1,'polynomial');   xvv = rand_const_const_partial({s1},1,'polynomial');
        for kv=1:2
            if kv==1,   val = 2*w;          Yb = 2*xw;
            else,       val = w+2*v;        Yb = xw+2*xvv;      end
            fnew = subs(f,{'b'},val);
            Yref = eval_polyopvar(f,{xv.a;Yb;xv.c});
            xs = struct('a',xv.a,'c',xv.c,'w',xw,'v',xvv);
            assert_close_const_partial(eval_new_const_partial(fnew,xs,sz,allnames),Yref,tol,lbl+sprintf(" b=val%d",kv));
            ncases = ncases+1;
        end
    end
    % a substitute on other variables than the state: b(s1) by a function of s2, a(s1,s2) by a function of s1
    % only, or by a function of a variable that f does not have
    f = make_f_const_partial(N,[1,2,3],dep,sz,'pw');
    expect_error_const_partial(@() subs(f,{'b'},known_fun(1+s2,{'s2'})),'Output space');
    expect_error_const_partial(@() subs(f,{'a'},known_fun(1+s1,{'s1'})),'Output space');
    expect_error_const_partial(@() subs(f,{'b'},known_fun(1+s1,{'s3'})),'Output space');
    ncases = ncases+3;

    %% (H) errors
    N = 2;      svs = svar_all;     dep = [true,true; true,false; false,false];     sz = [1;1;1];
    f = make_f_const_partial(N,[1,2,3],dep,sz,'pw');
    % polynomial in a variable that the state does not depend on
    expect_error_const_partial(@() subs(f,{'b'},s2),'spatial variables');
    expect_error_const_partial(@() subs(f,{'c'},s1),'spatial variables');
    expect_error_const_partial(@() subs(f,2,s2),'spatial variables');
    % wrong number of rows, unknown state, factor that is no state, factor beyond the last one
    expect_error_const_partial(@() subs(f,{'a'},[1;2]),'Number of rows');
    expect_error_const_partial(@() subs(f,{'d'},1),'states');
    expect_error_const_partial(@() subs(f,5,1),'degree');
    % all states replaced in a function with a Kronecker output is still not supported
    fk = make_f_const_partial(N,[1,2],dep,sz,'kron');
    expect_error_const_partial(@() subs(fk,{'a';'b'},[1;2]),'tensor product');
    ncases = ncases+7;

    fprintf('case_const_partial passed (%d cases).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_fctr()
% case_fctr -- substitution
% fnew = subs_factor(f,fctr,val), replacing factor number fctr by val,
% for a polynomial functional f (intop coefficient, no output variables).
%
%
% NOTES:
% - Structural check only: fnew must be a validly declared 'polyopvar';
%   its value is not compared with an independent evaluation (see
%   Test_polyopvar_subs/case_const* for value checks).

% % Set parameters defining f and val

    nterms = 1;     % number of terms in val
    nfctrs = 1;     % number of factors in polynomial to substitute
    d1 = 2;         % degree of monomial defining f
    d2 = 1;         % maximal degree of polynomial val
    nx1 = 1;        % number of state variables in f
    nx2 = 1;        % number of state variables in val
    nx_share = min(nx1,nx2);  % number of state varaibles in val that also appear in f
    m = 1;          % maximal vector size of the states
    Ntot = 1;       % number of spatial variables on which f depends

    dnsty = 0.25;

    % Set parameters defining the functional
    otype = [2,2];
    nr = 1;         % number of rows in f
    use_full_int = false;
    pdegv = 3;

    % % Declare the state variable names
    xidcs_full = randperm(nx1+nx2,nx1+nx2);
    xidcs1 = sort(xidcs_full(1:nx1));
    xidcs2 = sort(xidcs_full(nx1+1:end));
    xidcs_full = [xidcs1; xidcs2];
    xshare_idcs1 = sort(randperm(nx1,nx_share));
    xshare_idcs2 = sort(randperm(nx2,nx_share));
    xidcs2(xshare_idcs2) = xidcs1(xshare_idcs1);
    xvarname_full = cell(nx1+nx2,1);
    % All nx1+nx2 entries must be filled: xidcs1/xidcs2 index into them by values up to nx1+nx2.
    for j=1:nx1+nx2
        xvarname_full{j} = ['x',num2str(xidcs_full(j))];
    end
    xvarname1 = xvarname_full(xidcs1);
    xvarname2 = xvarname_full(xidcs2);

    % % Set the size of the state variables
    dims_full = randi(m,nx1+nx2,1);
    dims1 = dims_full(xidcs1);
    dims2 = dims_full(xidcs2);

    % % Declare spatial variable names and domains
    svarname_full = cell(1,Ntot);
    for j=1:Ntot
        svarname_full{j} = ['s',num2str(j)];
    end
    dom_full = randi([0,1],[Ntot,1]) +[zeros(Ntot,1),rand([Ntot,1])];

    % % Set dependency arrays for the variables
    depmat_full = logical(randi([0,1],nx1+nx2,Ntot));
    depmat1 = depmat_full(xidcs1,:);
    depmat2 = depmat_full(xidcs2,:);

    % % Set the variables and domains actually appearing in f and val
    isvar1 = any(depmat1,1);                isvar2 = any(depmat2,1);
    svars1 = svarname_full(isvar1);         svars2 = svarname_full(isvar2);
    doms1 = dom_full(isvar1,:);             doms2 = dom_full(isvar2,:);
    depmat1 = depmat1(:,isvar1);            depmat2 = depmat2(:,isvar2);
    %degs1 = sum(depmat1,1);
    pdegs = randi(2,[1,size(depmat2,2)]);
    pdegsB = struct('in',randi(2,[1,size(depmat1,2)]),'out',randi(2,[1,size(depmat2,2)]));



    % % Declare the polynomial function f
    % Set the monomial defining f
    fZ = randi([1,nx1],1,d1);
    depmat_mid = depmat1(fZ,:);
    % degs1 (degree of each spatial variable in the monomial Aop acts on, see
    % rand_intop) must count repeated occurrences of a state in fZ; sum(depmat1,1)
    % would build Aop for too low a degree and trip the input/output space check
    % in @intop/mtimes_tensopvar.
    degs1 = sum(depmat_mid,1);
    dims_mid = reshape(dims1(fZ),[],1);
    m1 = prod(dims1(fZ));
    % Generate the functional acting on the monomial
    Aop = rand_intop([nr,m1],svars1,doms1,degs1,pdegs);
    % Generate f
    f = polyopvar();
    f.Z = fZ;
    f.C = tensopmat(Aop);
    f.varname = xvarname1;
    f.dims.in = dims1;          f.dims.out = nr;
    f.vars.in = svars1;         f.vars.out = cell(1,0);
    f.dom.in = doms1;           f.dom.out = zeros(0,2);
    f.depmat.in = depmat1;      f.depmat.out = zeros(1,0);



    % % Declare the polynomial val
    fctrs = randperm(d1,nfctrs);
    val = [];
    Bop_cell = cell(nfctrs,1);
    vars = struct('in',{svars2},'out',{svars1});
    dom = struct('in',doms2,'out',doms1);
    depmat = struct('in',depmat2,'out',zeros(1,0));
    for k=1:numel(fctrs)
        % Determine the function space of factor k
        idx = fZ(fctrs(k));
        depmat_k = depmat;
        depmat_k.out = depmat1(idx,:);
        m_k = dims1(idx);
        % Declare a new monomial to substitute
        Z_k = randi([1,nx2],[1,d2]);
        dims_k_2 = dims2(Z_k);
        depmat_k.in = depmat2(Z_k,:);
        % Declare a tensopvar object mapping monomials to factor
        dims_k_1 = round(m_k/d2)*ones(d2,1);
        if otype(1)==2 && prod(dims_k_1,1)~=m_k
            dims_k_1(1) = m_k;
            dims_k_1(2:end) = 1;
        elseif otype(1)==1
            dims_k_1 = m_k*ones(d2,1);
        end
        dims_k = [dims_k_1,dims_k_2];
        Bop_k = rand_tensopvar(dims_k,vars,dom,depmat_k,nterms,otype,pdegsB);
        Bop_cell{k} = Bop_k;
        % Declare the polynomial for substitution
        val_k = polyopvar();
        val_k.Z = Z_k;
        val_k.C = tensopmat(Bop_k);
        val_k.varname = xvarname2;
        val_k.dims.in = dims2;          val_k.dims.out = m_k;
        val_k.vars.in = svars2;         val_k.vars.out = svars1;
        val_k.dom.in = doms2;           val_k.dom.out = doms1;
        val_k.depmat.in = depmat1;      val_k.depmat.out = depmat_k.out;
        val_k.chkval = true;
        % Concatenate
        val = [val; val_k];
    end


    % % Perform the substitution
    fnew = subs(f,fctrs,val);
    if ~isa(fnew,'polyopvar')
        error("case_fctr failed: substitution did not return a 'polyopvar' object.")
    elseif ~fnew.chkval
        error("case_fctr failed: substitution result is not marked as validly-declared.")
    end
    disp('case_fctr passed (structural check only -- see note above).')
end

%% --------------------------------------------------------------------------
function case_fctr_function()
% case_fctr_function -- substitution of a factor in a
% polynomial function (tensopvar coefficient with output variables), as
% opposed to Test_polyopvar_subs/case_fctr (functionals, structural check
% only).
%
% Supported cases (see subs_factor in @polyopvar/subs.m):
%   (a) f and the substituted polynomial are both single-factor;
%   (b) f has d factors (pointwise [1,2] or Kronecker [2,2] output), and
%       each monomial of val is degree one, so the factors correspond
%       one-to-one.
% Everything else must raise an informative error.
%
% Verification is semantic and independent of subs: fnew(x,y) is compared
% with Kop applied to the factor list in which factor "fctr" has been
% replaced by the value val(y), evaluated by eval_polyopvar.
%
%
% NOTES:
% - Evaluation uses eval_polyopvar (test-only, does not use subs), not
%   apply_poly, which evaluates by subs and would make the check circular.

    pvar s1 s2
    svar_all = {s1,s2};
    tol = 1e-8;
    ncases = 0;

    for Ntot = 1:2
        for d = 1:3
            otypes = {[1,2]};
            if d>1,     otypes = {[1,2],[2,2]};     end
            for io = 1:numel(otypes)
                ot = otypes{io};
                for fctr = 1:d
                    for dval = 1:2          % number of monomials in val
                        rng(100*Ntot+10*d+fctr+io);
                        [f,val,Kop] = build_case(Ntot,d,ot,dval);
                        fnew = subs(f,fctr,val);
                        x = rand_pol_fctr_function(svar_all(1:Ntot));
                        y = rand_pol_fctr_function(svar_all(1:Ntot));
                        vy = eval_polyopvar(val,{y});
                        X = repmat({x},1,d);       X{fctr} = vy;
                        Yref = apply_op(Kop,X);
                        xv = cell(numel(fnew.varname),1);
                        for k=1:numel(xv)
                            if strcmp(fnew.varname{k},'x'),     xv{k} = x;
                            else,                               xv{k} = y;      end
                        end
                        Ynew = eval_polyopvar(fnew,xv);
                        % Guard against a vacuous comparison of two zero polynomials
                        Yc = polynomial(Yref).coeff;
                        if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
                            error("case_fctr_function: reference value is zero (Ntot=%d,d=%d,fctr=%d); test would be vacuous.",Ntot,d,fctr)
                        end
                        err = polynomial(Ynew-Yref);
                        if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
                            error("case_fctr_function failed (Ntot=%d,d=%d,otype=[%d,%d],fctr=%d,dval=%d): residual %.3e.",...
                                Ntot,d,ot(1),ot(2),fctr,dval,max(abs(full(err.coeff(:)))))
                        end
                        ncases = ncases+1;
                    end
                end
            end
        end
    end

    % % Unsupported cases must error informatively
    % (1) single-factor f, val containing a degree-2 monomial
    rng(7);
    [f,val] = build_case(1,1,[1,2],1);
    val2 = build_val_deg2(1);
    expect_error_fctr_function(@() subs(f,1,val2),'supported');
    % (2) two-factor pointwise f, val with a degree-2 monomial
    [f,~] = build_case(1,2,[1,2],1);
    expect_error_fctr_function(@() subs(f,1,val2),'supported');
    % (3) multi-row f is supported, see Test_polyopvar_subs/case_multirow

    fprintf('case_fctr_function passed (%d substitution cases).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_multirow()
% case_multirow -- substitution (all routes of
% @polyopvar/subs.m) in polynomials with multiple rows of coefficient
% operators, built by vertcat of two polynomials (f = [f1;f2]), possibly
% with vector-valued states and rows depending on different output
% variables:
%   (1) subs_var by name, by a polynomial function val(y);
%   (2) subs_factor by position, same;
%   (3) subs_const by name and by position, 'double'/'polynomial' vectors.
%
% Verification is independent of subs: the substituted polynomial is
% evaluated at the remaining states and compared with the original f
% evaluated (eval_polyopvar) at a = val(y) or a = constant.
%
%
% NOTES:
% - Evaluation uses eval_polyopvar (test-only, does not use subs), not
%   apply_poly, which evaluates by subs and would make the check circular.

    pvar s1 s2
    svar_all = {s1,s2};
    tol = 1e-8;
    ncases = 0;
    sizes = [1,1,1; 2,3,2];             % [ma,mb,my]
    seed = 0;

    for Ntot = 1:2
        svs = svar_all(1:Ntot);
        % Row dependencies on the output variables: all, or (for Ntot=2) row 1
        % on s1 only and row 2 on s2 only
        outsels = {true(2,Ntot)};
        if Ntot==2,     outsels{2} = [true,false; false,true];      end
        for iout = 1:numel(outsels)
            outsel = outsels{iout};
            for isz = 1:size(sizes,1)
                ma = sizes(isz,1);      mb = sizes(isz,2);      my = sizes(isz,3);
                sz = [ma,mb,my];
                for kind = {'pw','kron','functional'}
                    if strcmp(kind{1},'functional') && iout>1,  continue,   end
                    for dval = 1:2
                        seed = seed+1;      rng(11000+seed);
                        lbl = sprintf('Ntot=%d,out=%d,[ma,mb,my]=[%d,%d,%d],%s,dval=%d',Ntot,iout,ma,mb,my,kind{1},dval);
                        if strcmp(kind{1},'functional'),    mono = 'a';     else,   mono = 'ab';    end
                        f = vertcat(build_f_multirow(Ntot,[ma,mb],kind{1},mono,outsel(1,:)),...
                                    build_f_multirow(Ntot,[ma,mb],kind{1},mono,outsel(2,:)));
                        if size(f.C.ops,1)~=2
                            error("case_multirow: test setup failed (%s): f should have 2 rows.",lbl)
                        end
                        val = build_val(Ntot,ma,my,dval);
                        y = rand_pol_multirow(svs,my);       b = rand_pol_multirow(svs,mb);
                        vy = eval_polyopvar(val,{y});
                        Yref = eval_polyopvar(f,{vy;b});
                        % (1) by name
                        fnew = subs(f,{'a'},val);
                        check(fnew,Yref,{'b',b,'y',y},sz,tol,lbl+" (1)");
                        % (2) by position (first factor)
                        fnew = subs(f,1,val);
                        check(fnew,Yref,{'b',b,'y',y},sz,tol,lbl+" (2)");
                        ncases = ncases+2;
                    end
                end

                %% (3) constants
                for cform = {'double','polynomial'}
                    rng(12000+100*Ntot+10*iout+isz+numel(cform{1}));
                    ca = rand_val_multirow(svs,ma,cform{1});     cb = rand_val_multirow(svs,mb,cform{1});
                    b = rand_pol_multirow(svs,mb);
                    lbl = sprintf('Ntot=%d,out=%d,[ma,mb]=[%d,%d],const %s',Ntot,iout,ma,mb,cform{1});
                    % two-factor pointwise function, all factors replaced
                    f = vertcat(build_f_multirow(Ntot,[ma,mb],'pw','ab',outsel(1,:)),build_f_multirow(Ntot,[ma,mb],'pw','ab',outsel(2,:)));
                    Yref = eval_polyopvar(f,{ca;cb});
                    fnew = subs(f,{'a';'b'},[ca;cb]);
                    check(fnew,Yref,{},sz,tol,lbl+" (3a) by name");
                    fnew = subs(f,[1,2],[ca;cb]);
                    check(fnew,Yref,{},sz,tol,lbl+" (3a) by position");
                    % single-factor function (pointwise, Kronecker) and functional
                    for kind = {'pw','kron','functional'}
                        if strcmp(kind{1},'functional') && iout>1,  continue,   end
                        f = vertcat(build_f_multirow(Ntot,[ma,mb],kind{1},'a',outsel(1,:)),build_f_multirow(Ntot,[ma,mb],kind{1},'a',outsel(2,:)));
                        Yref = eval_polyopvar(f,{ca;b});
                        if strcmp(kind{1},'kron')
                            Yref.varname = regexprep(Yref.varname,"_1$","");
                        end
                        fnew = subs(f,{'a'},ca);
                        check(fnew,Yref,{},sz,tol,lbl+" (3b) "+kind{1});
                        ncases = ncases+1;
                    end
                    % several monomials, of which only a is substituted: f = [a-part; b-part] terms
                    fa = vertcat(build_f_multirow(Ntot,[ma,mb],'pw','a',outsel(1,:)),build_f_multirow(Ntot,[ma,mb],'pw','a',outsel(2,:)));
                    fb = vertcat(build_f_multirow(Ntot,[ma,mb],'pw','b',outsel(1,:)),build_f_multirow(Ntot,[ma,mb],'pw','b',outsel(2,:)));
                    f = fa + fb;
                    Yref = eval_polyopvar(f,{ca;b});
                    fnew = subs(f,{'a'},ca);
                    check(fnew,Yref,{'b',b},sz,tol,lbl+" (3c) kept monomials");
                    ncases = ncases+3;
                end
            end
        end
    end

    fprintf('case_multirow passed (%d substitution cases).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_name()
% case_name -- substitution of states by name, and of
% several factors at once, in polynomial functionals
% f = sum_k Aop_k*Z_k(x) (intop coefficients, as V(x,u) of piesos_LF),
% see subs_var, subs_fctrs, subs_factor in @polyopvar/subs.m:
%   (i)-(iii) the states a and/or b of f(a,b) replaced by a function
%         val(y) of a new state y, by name, for monomials with repeated
%         factors, e.g. a o a, and several states, a o b o a (i), of
%         different degree, with Z padded with zeros (ii), and kernels
%         that are 'dpvar' with a different list of decision variables
%         in each region of integration, as produced by piesos_quadvar
%         (iii);
%   (iv)  several factors in a monomial, by number;
%   (v)   the new value containing the substituted state itself,
%         a --> val(a);
%   (vi)  renaming of states, also to the name of another state;
%   (vii) repeated factors in polynomial functions (tensopvar
%         coefficients), f(x) = Kop*(x o x);
%   errors for repeated or too large factor numbers, and for states on
%   different function spaces.
%
% Verification is semantic and independent of subs: fnew(y,b) is compared
% to the original f evaluated at the states a=val(y), b (evaluating the
% functional by integration, apply_op, of the kernels at the products of
% the polynomial states), for random polynomial states.
%

% NOTES:
% - All kernels of an 'intop' share one list of decision variables Zd, so
%   "a different list in each region" is realized by thin_dvars: a random
%   subset of Zd with nonzero coefficients.
% - States are evaluated with eval_polyopvar (Testfolder/polyopvar), not
%   apply_poly, since apply_poly evaluates by subs, the routine under
%   test.
% - eval_functional must name the dummy variables s_k_j exactly as
%   apply_op does (matched by name).
% - assert_close and (vii) reject a near-zero reference value, so no case
%   can pass vacuously.

    pvar s1 s2
    svar_all = {s1,s2};
    tol = 1e-8;
    ncases = 0;

    %% (i)-(iii) substitution of one or both states of f(a,b), for a function val(y), by name
    %   (Z rows padded with zeros mix the degrees, rows with repeated a, b)
    Zlist = {   [1,1];
                [1,2;  2,1];
                [1,1,2;  2,0,0;  1,2,0];
                [1,1,2;  2,2,1;  1,1,1]   };
    for Ntot = 1:2
        for iz = 1:numel(Zlist)
            Z = Zlist{iz};
            if Ntot==2 && size(Z,2)>2,  continue,   end         % skip monomials of 3 factors in 2D
            for q = [0,7]
                for dval = 1:2
                    seed = 1000*Ntot + 100*iz + 10*dval + q;
                    rng(seed);
                    [svars,dom] = spatial_space(Ntot);
                    f = build_f_name(svars,dom,{'a';'b'},Z,2,q);
                    val = [build_val_name(svars,dom,'y',dval); build_val_name(svars,dom,'y',dval)];
                    svs = svar_all(1:Ntot);
                    d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:max(q,1),'UniformOutput',false),num2cell(randn(max(q,1),1)));
                    label = sprintf('Ntot=%d,Z=%d,q=%d,dval=%d',Ntot,iz,q,dval);
                    y = rand_pol_name(svs);      b = rand_pol_name(svs);
                    va = eval_polyopvar(val(1,:),{y});     vb = eval_polyopvar(val(2,:),{y});

                    % (a) substitute a only: val is a single row
                    fnew = subs(f,{'a'},val(1,:));
                    no_state_left(fnew,'a',label);
                    Yref = eval_functional_name(f,struct('a',va,'b',b),d0map);
                    Ynew = eval_functional_name(fnew,struct('y',y,'b',b),d0map);
                    assert_close_name(Ynew,Yref,tol,['(a) ',label]);
                    ncases = ncases+1;

                    % (b) substitute both a and b, listed in a different order than the states
                    fnew = subs(f,{'b','a'},[val(2,:);val(1,:)]);
                    no_state_left(fnew,'a',label);       no_state_left(fnew,'b',label);
                    Yref = eval_functional_name(f,struct('a',va,'b',vb),d0map);
                    Ynew = eval_functional_name(fnew,struct('y',y),d0map);
                    assert_close_name(Ynew,Yref,tol,['(b) ',label]);
                    ncases = ncases+1;
                end
            end
        end
    end

    %% (iv) substitution of several factors in a single monomial by index
    for Ntot = 1:2
        for Zc = {[1,2,0], [1,1], [2,1,2]}
            Z = Zc{1};
            if Ntot==2 && size(Z,2)>2,  continue,   end
            for q = [0,5]
                rng(2000+10*Ntot+q+size(Z,2));
                [svars,dom] = spatial_space(Ntot);
                svs = svar_all(1:Ntot);
                f = build_f_name(svars,dom,{'a';'b'},Z,2,q);
                d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:max(q,1),'UniformOutput',false),num2cell(randn(max(q,1),1)));
                y = rand_pol_name(svs);      b = rand_pol_name(svs);
                nfac = nnz(Z);
                val = build_val_name(svars,dom,'y',1);
                for k=2:nfac,   val = [val; build_val_name(svars,dom,'y',1)];     end %#ok<AGROW>
                vy = cell(nfac,1);
                for k=1:nfac,   vy{k} = eval_polyopvar(val(k,:),{y});    end
                label = sprintf('Ntot=%d,Z=[%s],q=%d',Ntot,num2str(Z),q);
                % all factors
                fnew = subs(f,1:nfac,val);
                Yref = eval_functional_fctrs(f,Z,vy,1:nfac,struct('a',b,'b',b),d0map);
                Ynew = eval_functional_name(fnew,struct('y',y),d0map);
                assert_close_name(Ynew,Yref,tol,['(iv) all factors, ',label]);
                ncases = ncases+1;
                % the last factor only, the others remain states a, b (both given the value b)
                fnew = subs(f,nfac,val(nfac,:));
                Yref = eval_functional_fctrs(f,Z,vy(nfac),nfac,struct('a',b,'b',b),d0map);
                Ynew = eval_functional_name(fnew,struct('y',y,'a',b,'b',b),d0map);
                assert_close_name(Ynew,Yref,tol,['(iv) last factor, ',label]);
                ncases = ncases+1;
            end
        end
    end

    %% (v) the new value contains the substituted state itself, a --> val(a), so that all of the
    %   factors must be substituted at once (the substituted factors may not be substituted again)
    for q = [0,5]
        for Zc = {[1,1], [1,2,1]}
            Z = Zc{1};
            rng(2500+q+size(Z,2));
            [svars,dom] = spatial_space(1);
            svs = svar_all(1);
            f = build_f_name(svars,dom,{'a';'b'},Z,2,q);
            d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:max(q,1),'UniformOutput',false),num2cell(randn(max(q,1),1)));
            a = rand_pol_name(svs);      b = rand_pol_name(svs);
            val = build_val_name(svars,dom,'a',2);
            va = eval_polyopvar(val,{a});
            label = sprintf('Z=[%s],q=%d',num2str(Z),q);
            fnew = subs(f,{'a'},val);
            Yref = eval_functional_name(f,struct('a',va,'b',b),d0map);
            Ynew = eval_functional_name(fnew,struct('a',a,'b',b),d0map);
            assert_close_name(Ynew,Yref,tol,['(v) a --> val(a), ',label]);
            ncases = ncases+1;
        end
    end

    %% (vi) renaming of states. For Ntot=2, the state a depends on s1 and s2, but b only on s1,
    %   so that the function spaces must move with the names
    for Ntot = 1:2
        for q = [0,5]
            rng(3000+10*Ntot+q);
            [svars,dom] = spatial_space(Ntot);
            svs = svar_all(1:Ntot);
            depmat = true(2,Ntot);
            if Ntot==2,     depmat(2,2) = false;    end
            Zren = [1,2,0;  1,1,0;  1,1,2];         % each monomial contains a, so depends on all of the variables
            f = build_f_name(svars,dom,{'a';'b'},Zren(1:(2+(Ntot==1)),:),2,q,depmat);
            d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:max(q,1),'UniformOutput',false),num2cell(randn(max(q,1),1)));
            a = rand_pol_name(svs);      b = rand_pol_name(svs(depmat(2,:)));
            label = sprintf('Ntot=%d,q=%d',Ntot,q);
            % a --> z moves the state after b
            fnew = subs(f,{'a'},{'z'});
            assert(isequal(fnew.varname,{'b';'z'}),"case_name (%s): states renamed or not sorted.",label)
            Yref = eval_functional_name(f,struct('a',a,'b',b),d0map);
            Ynew = eval_functional_name(fnew,struct('z',a,'b',b),d0map);
            assert_close_name(Ynew,Yref,tol,['(vi) a --> z, ',label]);
            ncases = ncases+1;
            % a,b --> b,a swaps the states
            fnew = subs(f,{'a','b'},{'b','a'});
            Yref = eval_functional_name(f,struct('a',a,'b',b),d0map);
            Ynew = eval_functional_name(fnew,struct('b',a,'a',b),d0map);
            assert_close_name(Ynew,Yref,tol,['(vi) a,b --> b,a, ',label]);
            ncases = ncases+1;
            if Ntot==1
                % a --> b merges the states
                fnew = subs(f,{'a'},{'b'});
                assert(isequal(fnew.varname,{'b'}),"case_name (%s): states not merged.",label)
                Yref = eval_functional_name(f,struct('a',b,'b',b),d0map);
                Ynew = eval_functional_name(fnew,struct('b',b),d0map);
                assert_close_name(Ynew,Yref,tol,['(vi) a --> b, ',label]);
                ncases = ncases+1;
            else
                % states on different function spaces cannot be merged
                expect_error_name(@() subs(f,{'a'},{'b'}),'same function space');
            end
        end
    end

    %% (vii) polynomial FUNCTIONS (tensopvar coefficients) with repeated factors, f(x) = Kop*(x o ... o x),
    %   x --> val(y) by name, for pointwise and Kronecker outputs, compared to Kop applied to val(y)
    for Ntot = 1:2
        for d = 2:(4-Ntot)          % evaluating the products of 3 factors in 2D takes long
            for otc = {[1,2],[2,2]}
                ot = otc{1};
                for dval = 1:2
                    rng(100*Ntot+10*d+ot(1)+dval);
                    [f,val,Kop] = build_function_case(Ntot,d,ot,dval);
                    x = rand_pol_name(svar_all(1:Ntot));         y = rand_pol_name(svar_all(1:Ntot));
                    vy = eval_polyopvar(val,{y});
                    Yref = apply_op(Kop,repmat({vy},1,d));
                    fnew = subs(f,{'x'},val);
                    xv = cell(numel(fnew.varname),1);
                    for k=1:numel(xv)
                        if strcmp(fnew.varname{k},'y'),     xv{k} = y;
                        else,                               xv{k} = x;      end
                    end
                    Ynew = eval_polyopvar(fnew,xv);
                    label = sprintf('(vii) function, Ntot=%d,d=%d,otype=[%d,%d],dval=%d',Ntot,d,ot(1),ot(2),dval);
                    Yc = polynomial(Yref).coeff;
                    if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
                        error("case_name (%s): reference value is zero; test would be vacuous.",label)
                    end
                    err = polynomial(Ynew-Yref);
                    if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
                        error("case_name failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
                    end
                    ncases = ncases+1;
                end
            end
        end
    end

    %% Errors
    rng(4000);
    [svars,dom] = spatial_space(1);
    f = build_f_name(svars,dom,{'a';'b'},[1,2],2,0);
    val = build_val_name(svars,dom,'y',1);
    expect_error_name(@() subs(f,[1,1],[val;val]),'unique');
    expect_error_name(@() subs(f,[1,3],[val;val]),'degree');
    % states on different function spaces cannot be merged; depmat declared consistently with the
    % operators (assigning depmat.in alone sets chkval=false and initialize rejects the object)
    fc = build_f_name(svars,dom,{'a';'b'},[1,2],2,0,[true;false]);
    expect_error_name(@() subs(fc,{'a'},{'b'}),'same function space');

    fprintf('case_name passed (%d substitution cases).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_nested()
% case_nested -- substitution of some, or all, of the
% states in a monomial acted on by a nested operator (a factor of the
% 'tensopvar' is itself a 'tensopvar' acting on several of the factors of
% the monomial): the sum of monomials with the same Z and the same number
% of factors but multipliers m_k from R depending on different variables,
% x*m1*y + x*m2*y. subs_factor composes the coefficient operator with an
% operator that has a factor for each leaf factor of the monomial
% (@tensopvar/mtimes, nested A).
%
% Verification is semantic and independent of subs: the polynomial is
% evaluated (eval_polyopvar, test-only, one input for each leaf factor of
% the operator) after the substitution at the remaining states, and
% compared with the product of the polynomials x*(m1+m2)*y with the
% substituted values, as polynomials.
% Parts:
%   (1) substitution of x, y, both, and x then y, by constants
%       (double/polynomial);
%   (2) substitution by position, by a constant;
%   (3) substitution by polynomial functions of new states (a sum of
%       degree-one monomials in the new state, so that the operator of
%       the substituted factor has several terms), of x, y, and x then y;
%   (4) a constant for one state, a function for the other.
%

    tol = 1e-8;
    ncases = 0;
    pvar s1 s2
    m1 = 1-s1^2;        m2 = 2+s2+s1*s2;        m3 = 3-s2;      m4 = 1+s1+s2^2;
    polyopvar x(s1,s2) y(s1,s2)
    Id = struct('x',x.C.ops{1,1},'y',y.C.ops{1,1});
    % Each row: Z of the monomials, factors of each monomial, coefficients
    specs = {   [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'}},                  [1,1]
                [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'}},                  [2,-1]
                [0,1,2],    {{m1,'x','y'},{m3,'x','y'}},                  [1,1]
                [1,2,0],    {{'x','y',m1},{'x','y',m3}},                  [1,1]
                [1,0,0,2],  {{'x',m1,m3,'y'},{'x',m3,m1,'y'}},            [1,1]
                [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'},{'x',m4,'y'}},     [1,1,1]   };
    rng(1);

    %% (1) constants, by name
    for k=1:size(specs,1)
        [f,facs,cf,Zrow] = build_nested(specs(k,:),Id);
        cx = rand_const_nested([s1,s2]);       cy = rand_const_nested([s1,s2]);
        lbl = sprintf('(1.%d)',k);
        g = subs(f,{'x'},cx);
        check_nested(g,struct('y',rand_const_nested([s1,s2])),Zrow,facs,cf,struct('x',cx),tol,[lbl ' x']);
        g = subs(f,{'y'},cy);
        check_nested(g,struct('x',rand_const_nested([s1,s2])),Zrow,facs,cf,struct('y',cy),tol,[lbl ' y']);
        % both, at once, and one after the other
        g = subs(f,{'x';'y'},[cx;cy]);
        check_nested(g,struct(),Zrow,facs,cf,struct('x',cx,'y',cy),tol,[lbl ' x,y']);
        g = subs(subs(f,{'x'},cx),{'y'},cy);
        check_nested(g,struct(),Zrow,facs,cf,struct('x',cx,'y',cy),tol,[lbl ' x then y']);
        g = subs(subs(f,{'y'},cy),{'x'},cx);
        check_nested(g,struct(),Zrow,facs,cf,struct('x',cx,'y',cy),tol,[lbl ' y then x']);
        ncases = ncases+5;
    end

    %% (2) by position in Z
    for k=1:size(specs,1)
        [f,facs,cf,Zrow] = build_nested(specs(k,:),Id);
        cx = rand_const_nested([s1,s2]);       cy = rand_const_nested([s1,s2]);
        lbl = sprintf('(2.%d)',k);
        g = subs(f,find(Zrow==1),cx);
        check_nested(g,struct('y',rand_const_nested([s1,s2])),Zrow,facs,cf,struct('x',cx),tol,[lbl ' x']);
        g = subs(f,find(Zrow==2),cy);
        check_nested(g,struct('x',rand_const_nested([s1,s2])),Zrow,facs,cf,struct('y',cy),tol,[lbl ' y']);
        ncases = ncases+2;
    end

    %% (3) polynomial functions of new states: x = val1(w1) (a sum of dval operators acting on w1), y = val2(w2)
    for k=1:size(specs,1)
        for dval=1:2
            rng(100*k+dval);
            [f,facs,cf,Zrow] = build_nested(specs(k,:),Id);
            valx = build_val_nested('w1',dval);        valy = build_val_nested('w2',dval);
            w1 = rand_const_nested([s1,s2]);           w2 = rand_const_nested([s1,s2]);
            vx = eval_polyopvar(valx,{w1});     vy = eval_polyopvar(valy,{w2});
            lbl = sprintf('(3.%d,%d)',k,dval);
            g = subs(f,{'x'},valx);
            check_nested(g,struct('w1',w1,'y',rand_const_nested([s1,s2])),Zrow,facs,cf,struct('x',vx),tol,[lbl ' x']);
            g = subs(f,{'y'},valy);
            check_nested(g,struct('w2',w2,'x',rand_const_nested([s1,s2])),Zrow,facs,cf,struct('y',vy),tol,[lbl ' y']);
            g = subs(subs(f,{'x'},valx),{'y'},valy);
            check_nested(g,struct('w1',w1,'w2',w2),Zrow,facs,cf,struct('x',vx,'y',vy),tol,[lbl ' x then y']);
            % by position
            g = subs(f,find(Zrow==1),valx);
            check_nested(g,struct('w1',w1,'y',rand_const_nested([s1,s2])),Zrow,facs,cf,struct('x',vx),tol,[lbl ' x by position']);
            ncases = ncases+4;
        end
    end

    %% (4) a constant for one state, a function for the other
    for k=1:size(specs,1)
        rng(300+k);
        [f,facs,cf,Zrow] = build_nested(specs(k,:),Id);
        valy = build_val_nested('w2',2);
        cx = rand_const_nested([s1,s2]);           w2 = rand_const_nested([s1,s2]);
        vy = eval_polyopvar(valy,{w2});
        lbl = sprintf('(4.%d)',k);
        g = subs(subs(f,{'x'},cx),{'y'},valy);
        check_nested(g,struct('w2',w2),Zrow,facs,cf,struct('x',cx,'y',vy),tol,[lbl ' x const, y function']);
        valx = build_val_nested('w1',2);           cy = rand_const_nested([s1,s2]);
        w1 = rand_const_nested([s1,s2]);           vx = eval_polyopvar(valx,{w1});
        g = subs(subs(f,{'y'},cy),{'x'},valx);
        check_nested(g,struct('w1',w1),Zrow,facs,cf,struct('x',vx,'y',cy),tol,[lbl ' y const, x function']);
        ncases = ncases+2;
    end

    fprintf('case_nested passed (%d cases).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_ode_rows()
% case_ode_rows -- substitution of the PDE states u in the
% monomials Z(u) by the rows of Tx=T*x, for PIEs with ODE and PDE states
% (rows of different spatial dependence: the ODE rows depend on no variable
% and the PDE rows on s), also with vector-valued states:
%   - library example 9:  ODE state (size 1) and PDE state (size 1);
%   - library example 8:  ODE state (size 4) and PDE state (size 2);
% for monomials of degree 1 (all rows) and of degree 1 in a subset of the
% states, by name (@polyopvar/subs.m, subs_var with val_rows), and of the
% inner product <Tx,Tx> of the rows (sosprogramming/innerprod.m: ODE rows
% contribute u'*u, PDE rows the integral of u'*u).
%
% Verification is independent of subs: the PDE states u0=Tx(x0) and the
% monomials Z(u0) are evaluated by apply_op (eval_polyopvar), for random
% PIE states x0, and compared with the substituted polynomial Z(Tx)
% evaluated at x0.
%

    tol = 1e-8;
    sv = polynomial({'s'});     sv = sv(1);
    ncases = 0;
    for ex = [9,8]
        rng(500+ex);
        PDE = examples_PDE_library_PIETOOLS_nouinput(ex);
        evalc('PIE = convert_PIETOOLS_PDE_nonlinear(PDE,''silent'');');
        sz = PIE.x_tab(:,2);
        nst = numel(sz);
        % Random PIE states: constants for ODE states, polynomials in s for PDE states
        X0 = cell(nst,1);
        for i=1:nst
            if PIE.x_tab(i,3)==0
                X0{i} = polynomial(randn(sz(i),1));
            else
                X0{i} = polynomial(randn(sz(i),1)) + randn(sz(i),1)*sv + randn(sz(i),1)*sv^2;
            end
        end
        Tx = PIE.T*PIE.x;
        u0 = eval_polyopvar(Tx,X0);
        row0 = cumsum([0;sz(:)]);
        U0 = cell(nst,1);
        for i=1:nst
            U0{i} = u0(row0(i)+1:row0(i+1));
        end
        xu = polyopvar(PIE.varname_pde,sz,{'s'},PIE.dom,logical(PIE.x_tab(:,3)));
        Zu = dmonomials(xu,1);
        % (a) all states
        Zref = eval_polyopvar(Zu,U0);
        Znew = subs(Zu,Zu.varname,Tx);
        Y = eval_polyopvar(Znew,state_values(Znew,PIE,X0,U0));
        check_equal(Y,Zref,tol,sprintf('example %d, all states',ex));
        ncases = ncases+1;
        % (b) substitution of a subset of the states: the others remain
        for i=1:nst
            vals = U0;      vals(i) = [];
            names = PIE.varname_pde(i);
            ridx = row0(i)+1:row0(i+1);
            % substitute only state i by its row(s) of Tx: the remaining PDE states stay
            Tx_i = Tx(ridx,:);
            Zi = subs(Zu,names,Tx_i);
            % evaluate: remaining PDE states are given by their values u0, the PIE states by x0
            Yi = eval_polyopvar(Zi,state_values(Zi,PIE,X0,U0));
            check_equal(Yi,Zref,tol,sprintf('example %d, only state %d substituted',ex,i));
            ncases = ncases+1;
        end
        % (c) inner product of the rows of Tx: ODE rows u'*u, PDE rows int u'*u ds
        V = innerprod(Tx,Tx);
        Vref = polynomial(0);
        for i=1:nst
            q = U0{i}.'*U0{i};
            if PIE.x_tab(i,3)~=0
                q = int(q,sv,0,1);
            end
            Vref = Vref + q;
        end
        vals = struct();
        for i=1:nst,     vals.(PIE.varname_pie{i}) = X0{i};     end
        Y = eval_V_ode_rows(V,vals);
        check_equal(Y,Vref,tol,sprintf('example %d, innerprod',ex));
        ncases = ncases+1;
    end

    fprintf('case_ode_rows passed (%d cases).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_var()
% case_var -- exercise the subs_var (variable-name
% substitution) path of @polyopvar/subs.m, as distinct from
% Test_polyopvar_subs/case_fctr, which only exercises subs_factor
% (position-index substitution).
%
% Oracle: fnew(y) is compared with Aop applied to the same y (apply_op),
% and f's own output space must be unchanged by the substitution.
%
%
% NOTES:
% - subs_var must align state names with common_states, not common_vars:
%   common_vars also unifies the output spaces, which would overwrite the
%   (empty) vars.out of the functional f with that of val (check (1)).
% - Evaluation uses eval_polyopvar (test-only, does not use subs), not
%   apply_poly, which evaluates by subs and would make the check circular.

    pvar s1
    svars = {'s1'};
    dom = [0,1];

    % f(x) = int_0^1 x(s1) ds1 : scalar functional, linear in state x
    rng(8001);
    Aop = rand_intop([1,1],svars,dom,1,3);
    f = polyopvar();
    f.Z = 1;
    f.C = tensopmat(Aop);
    f.varname = {'x'};
    f.dims.in = 1;          f.dims.out = 1;
    f.vars.in = svars;      f.vars.out = cell(1,0);
    f.dom.in = dom;         f.dom.out = zeros(0,2);
    f.depmat.in = true;     f.depmat.out = zeros(1,0);
    f.chkval = true;

    y = polyopvar('y',svars,dom);

    % fnew(y) = int_0^1 y(s1) ds1
    fnew = subs(f,{'x'},y);

    % (1) f's own output space must be preserved exactly (f is still a
    % scalar functional, regardless of val's own output space)
    if ~isequal(fnew.vars.out,cell(1,0)) || ~isequal(size(fnew.dom.out),[0,2]) ...
            || ~isequal(size(fnew.depmat.out),[1,0])
        error("case_var failed: substitution corrupted f's own output space (vars.out/dom.out/depmat.out).")
    end

    % (2) the substituted-out variable 'x' must no longer appear in Z
    x_idx = find(strcmp(fnew.varname,'x'));
    if ~isempty(x_idx) && any(fnew.Z(:)==x_idx)
        error("case_var failed: fnew.Z still references the substituted-out variable x.")
    end

    % (3) the result must independently validate
    fnew.chkval = false;
    fnew = initialize(fnew);

    % (4) numerically: fnew(y) must equal Aop*y for random y
    yval = randn*s1 + randn;
    xvals = cell(numel(fnew.varname),1);
    for kk = 1:numel(fnew.varname)
        if strcmp(fnew.varname{kk},'y')
            xvals{kk} = yval;
        else
            xvals{kk} = polynomial(0);   % unused slot (x), value irrelevant
        end
    end
    Yfnew = eval_polyopvar(fnew,xvals);
    Yref = apply_op(Aop,yval);
    err = double(Yfnew-Yref);
    if norm(err) > 1e-8
        error("case_var failed: fnew(y) does not match Aop*y numerically (residual %.3e).",norm(err))
    end

    disp('case_var passed.')
end

%% --------------------------------------------------------------------------
function case_var_function()
% case_var_function -- substitution of a state variable
% by name (subs_var in @polyopvar/subs.m) in
%   (1) polynomial functions f(a,b)(s) = Kop*(a o b), tensopvar
%       coefficient,
%   (2) polynomial functionals f(a,b) = Aop*(a o b), intop coefficient,
% by (i) a polynomial function val(y) in a new state y, (ii) a constant.
% The variable is given as cellstr or as 'polyopvar' (polyopvar('a',...)).
%
% Verification is semantic and independent of subs: fnew(b,y) is compared
% to the original f evaluated (eval_polyopvar) at a=val(y) -- itself
% evaluated by eval_polyopvar -- and the unchanged b.
%
%
% NOTES:
% - Evaluation uses eval_polyopvar (test-only, does not use subs), not
%   apply_poly, which evaluates by subs and would make the check circular.
% - Kronecker outputs: the reference evaluation renames output variables
%   (suffix "_1"); stripped before comparing in (ii).

    pvar s1 s2
    svar_all = {s1,s2};
    tol = 1e-8;
    ncases = 0;

    %% (i) substitution by a polynomial function of a new state y
    for Ntot = 1:2
        for kind = {'function_pw','function_kron','functional'}
            for dval = 1:2
                for varform = {'cellstr','polyopvar'}
                    rng(1000*Ntot+100*dval+10*numel(kind{1})+numel(varform{1}));
                    [f,val,svars,dom] = build_f_var_function(Ntot,kind{1},false,dval);
                    if strcmp(varform{1},'cellstr')
                        var = {'a'};
                    else
                        var = polyopvar('a',svars,dom.in);
                    end
                    fnew = subs(f,var,val);
                    check_subs(f,fnew,val,svar_all(1:Ntot),tol,...
                        sprintf('(i) Ntot=%d,%s,dval=%d,var=%s',Ntot,kind{1},dval,varform{1}));
                    ncases = ncases+1;
                end
            end
        end
    end

    %% (ii) substitution by a constant, for a state a not depending on any spatial variable
    % The constant replaces the entire (degree-one) monomial, giving a degree-0 polynomial,
    % or some factors of a monomial of higher degree, which remain as factors that are multipliers
    % from R (see polyopvar), also for a 'double' constant (partial substitution in subs_const).
    % Verified against the original polynomial evaluated at the constant.
    for Ntot = 1:2
        for ot = {[1,2],[2,2]}
            for cform = {'polyopvar','double'}
                rng(5000+100*Ntot+10*ot{1}(1)+numel(cform{1}));
                c = randn;
                if strcmp(cform{1},'polyopvar')
                    val = constant2polyopvar(c);
                else
                    val = c;
                end
                % Degree one in a: f(a)(s) = Kop*a
                f = build_const_f(Ntot,1,ot{1});
                fnew = subs(f,{'a'},val);
                Yref = eval_polyopvar(f,{polynomial(c)});
                if strcmp(cform{1},"double")    % subs_const returns the output in the original variables
                    Yref.varname = regexprep(Yref.varname,"_1$","");   % apply_tensopvar renames output variables for Kronecker outputs
                end
                Ynew = eval_polyopvar(fnew,cell(numel(fnew.varname),1));
                assert_close_var_function(Ynew,Yref,tol,sprintf('(ii) Ntot=%d,ot=[%d,%d],val=%s',Ntot,ot{1}(1),ot{1}(2),cform{1}));
                ncases = ncases+1;
                % Degree two in a: a constant polyopvar replaces both factors, leaving a degree-0
                % term with two factors (a 'double' constant is handled by subs_const instead,
                % see Test_polyopvar_subs/case_const)
                if strcmp(cform{1},'polyopvar')
                    f2 = build_const_f(Ntot,2,ot{1});
                    fnew2 = subs(f2,{'a'},val);
                    Yref2 = eval_polyopvar(f2,{polynomial(c)});
                    Ynew2 = eval_polyopvar(fnew2,cell(numel(fnew2.varname),1));
                    assert_close_var_function(Ynew2,Yref2,tol,sprintf('(ii) a o a, Ntot=%d,ot=[%d,%d]',Ntot,ot{1}(1),ot{1}(2)));
                    if any(any(fnew2.Z))
                        error("case_var_function failed: a o a with a constant should be a degree-0 term.")
                    end
                    ncases = ncases+1;
                end
                % Constant in a monomial also containing b: the constant remains as a factor, Z=[0,b],
                % for a constant polyopvar as well as for a 'double' (partial substitution in subs_const)
                for kind3 = {'function_pw','function_kron'}
                    f3 = build_f_var_function(Ntot,kind3{1},true,1);
                    fnew3 = subs(f3,{'a'},val);
                    b3 = rand_pol_var_function(svar_all(1:Ntot));
                    Yref3 = eval_polyopvar(f3,{polynomial(c);b3});
                    xv3 = cell(numel(fnew3.varname),1);
                    for k3=1:numel(xv3)
                        if strcmp(fnew3.varname{k3},'b'),   xv3{k3} = b3;
                        else,                               xv3{k3} = polynomial(0);    end
                    end
                    assert_close_var_function(eval_polyopvar(fnew3,xv3),Yref3,tol,sprintf('(ii) a o b, Ntot=%d,%s,val=%s',Ntot,kind3{1},cform{1}));
                    if ~isequal(full(fnew3.Z),[0,2]) && ~isequal(full(fnew3.Z),[0,1])   % constant, then b
                        error("case_var_function failed: a o b with a constant should have Z=[0,b], not [%s].",num2str(full(fnew3.Z)))
                    end
                    ncases = ncases+1;
                end
            end
        end
    end

    fprintf('case_var_function passed (%d substitution cases).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_vector_states()
% case_vector_states -- substitution (all routes of
% @polyopvar/subs.m) in polynomials of vector-valued states:
%   (1) subs_var by name, by a polynomial function val(y) of a vector
%       state y;
%   (2) subs_factor by position, same;
%   (3) subs_const by name and by position, with 'double'/'polynomial'
%       vectors.
% States a,b,y have sizes ma,mb,my; f is a function (pointwise or
% Kronecker output type) or a functional (degree one), in 1 or 2 spatial
% variables.
%
% Verification is independent of subs: the substituted polynomial is
% evaluated at the remaining states and compared with the original f
% evaluated (eval_polyopvar, apply_op) at a = val(y) or a = constant.
%
%
% NOTES:
% - Evaluation uses eval_polyopvar (test-only, does not use subs), not
%   apply_poly, which evaluates by subs and would make the check circular.

    pvar s1 s2
    svar_all = {s1,s2};
    tol = 1e-8;
    ncases = 0;
    sizes = [2,1,1; 1,3,2; 2,3,2];          % [ma,mb,my]
    seed = 0;

    for Ntot = 1:2
        svs = svar_all(1:Ntot);
        for isz = 1:size(sizes,1)
            ma = sizes(isz,1);      mb = sizes(isz,2);      my = sizes(isz,3);
            for kind = {'pw','kron','functional'}
                for dval = 1:2
                    seed = seed+1;      rng(7000+seed);
                    lbl = sprintf('Ntot=%d,[ma,mb,my]=[%d,%d,%d],%s,dval=%d',Ntot,ma,mb,my,kind{1},dval);
                    if strcmp(kind{1},'functional')
                        f = build_f_vector_states(Ntot,[ma,mb],kind{1},'a');
                    else
                        f = build_f_vector_states(Ntot,[ma,mb],kind{1},'ab');
                    end
                    val = build_val_vector_states(Ntot,ma,my,dval);
                    y = rand_pol_vector_states(svs,my);       b = rand_pol_vector_states(svs,mb);
                    vy = eval_polyopvar(val,{y});
                    if strcmp(kind{1},'functional')
                        % degree-one functional in a: f(a,b) = Aop*a
                        f = build_f_vector_states(Ntot,[ma,mb],kind{1},'a');
                        Yref = eval_polyopvar(f,{vy;b});
                        % (1) by name
                        fnew = subs(f,{'a'},val);
                        check_vector_states(fnew,Yref,{'b',b,'y',y},[ma,mb,my],tol,lbl+" (1)");
                        % (2) by position
                        fnew = subs(f,1,val);
                        check_vector_states(fnew,Yref,{'b',b,'y',y},[ma,mb,my],tol,lbl+" (2)");
                        ncases = ncases+2;
                    else
                        Yref = eval_polyopvar(f,{vy;b});
                        % (1) by name
                        fnew = subs(f,{'a'},val);
                        check_vector_states(fnew,Yref,{'b',b,'y',y},[ma,mb,my],tol,lbl+" (1)");
                        % (2) by position (first factor of a o b)
                        fnew = subs(f,1,val);
                        check_vector_states(fnew,Yref,{'b',b,'y',y},[ma,mb,my],tol,lbl+" (2)");
                        ncases = ncases+2;
                    end
                end
            end

            %% (2b) by position in a o a: replace one factor, keep the other
            for kind = {'pw','kron'}
                rng(8000+10*Ntot+isz);
                lbl = sprintf('Ntot=%d,[ma,mb,my]=[%d,%d,%d],%s,a o a',Ntot,ma,mb,my,kind{1});
                [f,Kop] = build_f_vector_states(Ntot,[ma,mb],kind{1},'aa');
                val = build_val_vector_states(Ntot,ma,my,1);
                x = rand_pol_vector_states(svs,ma);       y = rand_pol_vector_states(svs,my);
                vy = eval_polyopvar(val,{y});
                for fctr = 1:2
                    X = {x,x};      X{fctr} = vy;
                    Yref = apply_op(Kop,X);
                    fnew = subs(f,fctr,val);
                    check_vector_states(fnew,Yref,{'a',x,'y',y},[ma,mb,my],tol,lbl+sprintf(" fctr=%d",fctr));
                    ncases = ncases+1;
                end
            end

            %% (3) constants, by name and by position
            for cform = {'double','polynomial'}
                rng(9000+10*Ntot+isz+numel(cform{1}));
                ca = rand_val_vector_states(svs,ma,cform{1});     cb = rand_val_vector_states(svs,mb,cform{1});
                b = rand_pol_vector_states(svs,mb);
                lbl = sprintf('Ntot=%d,[ma,mb]=[%d,%d],const %s',Ntot,ma,mb,cform{1});
                % two-factor pointwise function, all factors replaced
                f = build_f_vector_states(Ntot,[ma,mb],'pw','ab');
                Yref = eval_polyopvar(f,{ca;cb});
                fnew = subs(f,{'a';'b'},[ca;cb]);
                check_vector_states(fnew,Yref,{},[ma,mb,my],tol,lbl+" (3a) by name");
                fnew = subs(f,[1,2],[ca;cb]);
                check_vector_states(fnew,Yref,{},[ma,mb,my],tol,lbl+" (3a) by position");
                fnew = subs(f,{'b';'a'},[cb;ca]);
                check_vector_states(fnew,Yref,{},[ma,mb,my],tol,lbl+" (3a) reordered names");
                % single-factor function and functional
                for kind = {'pw','kron','functional'}
                    f = build_f_vector_states(Ntot,[ma,mb],kind{1},'a');
                    Yref = eval_polyopvar(f,{ca;b});
                    if strcmp(kind{1},'kron')
                        Yref.varname = regexprep(Yref.varname,"_1$","");
                    end
                    fnew = subs(f,{'a'},ca);
                    check_vector_states(fnew,Yref,{},[ma,mb,my],tol,lbl+" (3b) "+kind{1});
                    ncases = ncases+1;
                end
                % partial substitution: the factor replaced remains as a known factor
                f = build_f_vector_states(Ntot,[ma,mb],'pw','ab');
                Yref = eval_polyopvar(f,{ca;b});
                fnew = subs(f,{'a'},ca);
                check_vector_states(fnew,Yref,{'b',b},[ma,mb,my],tol,lbl+" (3c) partial, by name");
                fnew = subs(f,1,ca);
                check_vector_states(fnew,Yref,{'b',b},[ma,mb,my],tol,lbl+" (3c) partial, by position");
                ncases = ncases+2;
                % wrong number of rows (scalar for a vector state) errors
                if ma>1
                    expect_error_vector_states(@() subs(f,{'a';'b'},[ca(1);cb]),'Number of rows');
                end
                ncases = ncases+5;
            end
        end
    end

    fprintf('case_vector_states passed (%d substitution cases).\n',ncases)
end

%% ---- helpers --------------------------------------------------------------

function [svars,dom] = spatial_space(N)
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
end

function p = rand_pol(svs,m)
% Random m x 1 polynomial vector in the variables svs (cell of pvar)
p = polynomial(randn(m,1));
for i=1:numel(svs)
    p = p + randn(m,1)*svs{i} + randn(m,1)*svs{i}^2;
end
end

function v = rand_const(svs,m,form)
% Random m x 1 constant function: 'double', 'polynomial' in svs, or 'dpvar' (no decision variables)
if strcmp(form,'double')
    v = randn(m,1);
else
    v = rand_pol(svs,m);
    if strcmp(form,'dpvar'),    v = dpvar(v);   end
end
end

function f = make_f(N,Zrow,sz,kind,raw)
% Function f = Kop*(x_{Z1} o ... o x_{Zd}) in the states a,b,c of sizes sz, all depending on the
% variables s1,...,sN, with Kop random, with pointwise ('pw', output size 2) or Kronecker ('kron') output
[svars,dom] = spatial_space(N);
vars = struct('in',{svars},'out',{svars});
d = numel(Zrow);
din = reshape(sz(Zrow),[],1);
if strcmp(kind,'pw')
    dims = [2*ones(d,1),din];       ot = [1,2];     mout = 2;
else
    dims = [(1+(1:d)'),din];        ot = [2,2];     mout = prod(dims(:,1));
end
dp = struct('in',true(d,N),'out',true(d,N));
Kop = rand_tensopvar(dims,vars,dom,dp,1,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = Zrow;                 f.varname = {'a';'b';'c'};
f.dims.in = sz(:);          f.dims.out = mout;
f.vars.in = svars;          f.dom.in = dom.in;          f.depmat.in = true(3,N);
f.C = tensopmat(Kop);
f.vars.out = svars;         f.dom.out = dom.out;        f.depmat.out = Kop.deps.out;
if nargin<5 || ~raw,      f = initialize(f);      end
end

function V = make_V(svars,dom,Z,q)
% Functional V(a,b,c) = sum_k Aop_k*Z_k, Aop_k random intop, for the monomials in the rows of Z (padded with zeros)
N = numel(svars);
C = [];
for k=1:size(Z,1)
    idcs = Z(k,Z(k,:)~=0);
    Aop = rand_intop([1,1],svars,dom.in,numel(idcs)*ones(1,N),2,q);
    Ck = tensopmat(Aop);
    if k==1,    C = Ck;     else,   C = [C,Ck];     end %#ok<AGROW>
end
V = polyopvar();
V.Z = Z;                    V.C = C;                    V.varname = {'a';'b';'c'};
V.dims.in = ones(3,1);      V.dims.out = 1;
V.vars.in = svars;          V.vars.out = cell(1,0);
V.dom.in = dom.in;          V.dom.out = zeros(0,2);
V.depmat.in = true(3,N);    V.depmat.out = zeros(1,0);
V = initialize(V);
end

function Y = eval_V_tensor(V,cvals,svars)
% Value of the functional V at the constant states cvals: each monomial is applied to the Kronecker
% product of its factors (first factor slowest), the variable s of the j-th factor being named s_j.
% Independent of subs and of the handling of the cell of factors by apply_op of 'intop'
Y = 0;
Zf = full(V.Z);
for k=1:size(Zf,1)
    idcs = Zf(k,Zf(k,:)~=0);
    v = 1;
    for j=1:numel(idcs)
        xj = polynomial(cvals{idcs(j)});
        if numel(idcs)>1
            for i=1:numel(svars)
                xj.varname(strcmp(xj.varname,svars{i})) = {[svars{i},'_',num2str(j)]};
            end
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(V.C.ops{1,k},v);
end
Y = polynomial(Y);
end

function val = make_val(svars,dom,yname,sza,szy,dval)
% val(y)(s): vector of size sza, dval degree-one monomials in the new state y of size szy, random operators
N = numel(svars);
vars = struct('in',{svars},'out',{svars});
dep1 = struct('in',true(1,N),'out',true(1,N));
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([sza,szy],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok<AGROW>
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {yname};
val.dims.in = szy;          val.dims.out = sza;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = combine_terms(initialize(val));       % one monomial: vertcat (common_basis) drops duplicate monomials
end

function Y = eval_new(f,xvals,sz,allnames)
% Evaluate f at the values xvals.(name) of its states; states without a value get zeros
xv = cell(numel(f.varname),1);
for k=1:numel(xv)
    nm = f.varname{k};
    if isfield(xvals,nm)
        xv{k} = xvals.(nm);
    else
        xv{k} = polynomial(zeros(sz(strcmp(allnames,nm)),1));
    end
end
Y = eval_polyopvar(f,xv);
end

function n = check_const(f,sub,cvals,x,sz,tol,label)
% Substitute the states sub of f by the constants cvals, by name, and by position if f has a single monomial
% in which each is a single factor, and compare with f evaluated at the constants
allnames = {'a','b','c'};
p = randperm(numel(sub));       sub = sub(p);       cvals = cvals(p);       % the order of var and val need not be that of the states
cmat = vertcat(cvals{:});
xs = x;
for k=1:numel(sub)
    xs.(sub{k}) = cvals{k};
    if isa(xs.(sub{k}),'dpvar'),    xs.(sub{k}) = dpvar2poly(xs.(sub{k}));      end
end
Yref = eval_polyopvar(f,{xs.a;xs.b;xs.c});
fnew = subs(f,sub,cmat);
assert_close(eval_new(fnew,rmfield(x,sub),sz,allnames),Yref,tol,label+" (by name)");
n = 1;
Zr = full(f.Z);
idx = cellfun(@(s) find(strcmp(allnames,s)),sub);
if size(Zr,1)==1 && all(arrayfun(@(i) sum(Zr==i)==1,idx))
    pos = arrayfun(@(i) find(Zr==i),idx);
    fnew = subs(f,pos(:)',cmat);
    assert_close(eval_new(fnew,rmfield(x,sub),sz,allnames),Yref,tol,label+" (by position)");
    n = 2;
end
end

function svs = pv_cell(N)
% Cell with the pvar objects s1,...,sN
pvar s1 s2 s3
svs = {s1,s2,s3};
svs = svs(1:N);
end

function n = check_poly(f,sub,x,sz,N,dval,tol,label)
% Substitute the states sub of f by functions val_k(y) of a new state y, by name, and compare with f evaluated
% at the values val_k(y)
allnames = {'a','b','c'};
sub = sub(randperm(numel(sub)));        % the order of var and val need not be that of the states
[svars,dom] = spatial_space(N);
svs = pv_cell(N);
szy = 1+mod(N,2);
y = rand_pol(svs,szy);
val = [];       xs = x;
for k=1:numel(sub)
    sza = sz(strcmp(allnames,sub{k}));
    vk = make_val(svars,dom,'b2',sza,szy,dval);
    xs.(sub{k}) = eval_polyopvar(vk,{y});
    if k==1,    val = vk;   else,   val = [val;vk];     end %#ok<AGROW>
end
Yref = eval_polyopvar(f,{xs.a;xs.b;xs.c});
fnew = subs(f,sub,val);
xr = rmfield(x,sub);        xr.b2 = y;
assert_close(eval_new(fnew,xr,[sz;szy],[allnames,{'b2'}]),Yref,tol,label+" (val(y) by name)");
n = 1;
end

function n = check_pos(f,Z,pos,x,sz,N,tol,label)
% Substitute the factor pos of the monomial Z of f (single monomial) by a function of a new state y, and by a
% constant, and compare with Kop applied to the factors in which factor pos is replaced
allnames = {'a','b','c'};
[svars,dom] = spatial_space(N);
svs = pv_cell(N);
szy = 1+mod(N,2);
y = rand_pol(svs,szy);
nm = allnames{Z(pos)};      sza = sz(Z(pos));
fac = arrayfun(@(i) x.(allnames{i}),Z,'UniformOutput',false);
% function of y
val = make_val(svars,dom,'b2',sza,szy,1+mod(pos,2));
fac_y = fac;        fac_y{pos} = eval_polyopvar(val,{y});
Yref = apply_op(f.C,{fac_y});
fnew = subs(f,pos,val);
xr = x;     xr.b2 = y;
assert_close(eval_new(fnew,xr,[sz;szy],[allnames,{'b2'}]),Yref,tol,label+" (val(y))");
% constant
cval = rand_const(svs,sza,'polynomial');
fac_c = fac;        fac_c{pos} = cval;
Yref = apply_op(f.C,{fac_c});
fnew = subs(f,pos,cval);
assert_close(eval_new(fnew,x,sz,allnames),Yref,tol,label+" (constant)");
n = 2;
end

function n = check_poly_V(V,N,tol,label)
% Functional V: all states (by name) by functions of y, evaluated through the kernels (eval_polyopvar)
allnames = {'a','b','c'};
[svars,dom] = spatial_space(N);
svs = pv_cell(N);
y = rand_pol(svs,1);
val = [];       xs = struct();
for k=1:3
    vk = make_val(svars,dom,'b2',1,1,1+mod(k,2));
    xs.(allnames{k}) = eval_polyopvar(vk,{y});
    if k==1,    val = vk;   else,   val = [val;vk];     end %#ok<AGROW>
end
Yref = eval_functional(V,xs);
fnew = subs(V,allnames,val);
Ynew = eval_functional(fnew,struct('b2',y));
assert_close(Ynew,Yref,tol,label+" (all states by val(y))");
n = 1;
end

function Y = eval_functional(V,vals)
% Evaluate the functional V (polyopvar of intop coefficients) at the states vals.(name), polynomials in the spatial
% variables: the factors of each monomial are multiplied (kron), with the variable s_k of the factors that depend
% on it renamed s_k_j, j numbering these factors, as apply_op of 'intop' requires
vnames = V.vars.in;
dep_state = double(V.depmat.in)>0;
Zf = full(V.Z);
Y = 0;
for k=1:size(Zf,1)
    op = V.C.ops{1,k};
    if isempty(op),     continue,   end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    slot = zeros(1,numel(vnames));
    v = 1;
    for j=1:d
        xj = vals.(V.varname{idcs(j)});
        for mm=1:numel(vnames)
            if dep_state(idcs(j),mm)
                slot(mm) = slot(mm)+1;
                if d>1
                    xj.varname(strcmp(xj.varname,vnames{mm})) = {[vnames{mm},'_',num2str(slot(mm))]};
                end
            end
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end

function assert_close(Ynew,Yref,tol,label)
Yref = polynomial(Yref);        Ynew = polynomial(Ynew);
Yc = Yref.coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("case_class (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol*max(1,max(abs(full(Yc(:)))))
    error("case_class failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if isempty(str) || ~isempty(regexp(ME.message,str,'once'))
        return
    end
    error("case_class failed: unexpected error message '%s', expected '%s'.",ME.message,str)
end
error("case_class failed: expected an error (%s), but none was thrown.",str)
end

function f = build_f(N,ot,mono,kind)
% Polynomial in states a,b (both scalar, depending on all N variables):
%   'ab':   f = Kop*(a o b)          'aa':   f = Kop*(a o a)
%   'a':    f = Kop*a                'a+b':  f = Kop1*a + Kop2*b
% as function (pointwise/Kronecker type OT) or, for 'a', functional.
if nargin<4,    kind = 'function';  end
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
switch mono
    case 'ab'
        f = make_f_const([1,2],2,N,svars,dom,vars,ot,kind);
    case 'aa'
        f = make_f_const([1,1],2,N,svars,dom,vars,ot,kind);
    case 'a'
        f = make_f_const(1,1,N,svars,dom,vars,ot,kind);
    case 'a+b'
        f = make_f_const(1,1,N,svars,dom,vars,ot,kind) + make_f_const(2,1,N,svars,dom,vars,ot,kind);
end
end

function f = make_f_const(Z,d,N,svars,dom,vars,ot,kind)
f = polyopvar();
f.Z = Z;                    f.varname = {'a';'b'};
f.dims.in = [1;1];          f.dims.out = 1;
f.vars.in = svars;          f.dom.in = dom.in;
f.depmat.in = true(2,N);
if strcmp(kind,'function')
    dep = struct('in',true(d,N),'out',true(d,N));
    Kop = rand_tensopvar(ones(d,2),vars,dom,dep,1,ot,2);
    if ~Kop.chkval,     Kop = initialize(Kop);      end
    f.C = tensopmat(Kop);
    f.vars.out = svars;     f.dom.out = dom.out;    f.depmat.out = Kop.deps.out;
else
    Aop = rand_intop([1,1],svars,dom.in,ones(1,N),3);
    f.C = tensopmat(Aop);
    f.vars.out = cell(1,0); f.dom.out = zeros(0,2); f.depmat.out = zeros(1,0);
end
f = initialize(f);
end

function v = rand_val(svs,form)
% Random scalar constant (function) of the given type
if strcmp(form,'double')
    v = randn;
else
    v = polynomial(randn);
    for i=1:numel(svs)
        v = v + randn*svs{i} + randn*svs{i}^2;
    end
end
end

function Y = eval_new_const(fnew,vals)
% Evaluate fnew, giving values {'name',value,...} to the remaining states
xv = cell(numel(fnew.varname),1);
for k=1:numel(xv)
    xv{k} = polynomial(0);
    for i=1:2:numel(vals)
        if strcmp(fnew.varname{k},vals{i}),     xv{k} = vals{i+1};     end
    end
end
Y = eval_polyopvar(fnew,xv);
end

function assert_deg0(fnew,label)
if any(any(fnew.Z))
    error("case_const failed (%s): result should be a constant (degree-0) term.",label)
end
end

function assert_close_const(Ynew,Yref,tol,label)
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("case_const (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("case_const failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function expect_error_const(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("case_const failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("case_const failed: expected an error containing '%s', but none was thrown.",str)
end

function f = make_f_const_partial(N,Zrow,dep,sz,kind)
% Function f = Kop*(x_{Z1} o ... o x_{Zd}) in the states a,b,c of sizes sz, depending on the
% variables of dep (3 x N, row i the variables of state i), with Kop random, with pointwise
% ('pw', output size 2) or Kronecker ('kron') output
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
d = numel(Zrow);
din = reshape(sz(Zrow),[],1);
if strcmp(kind,'pw')
    dims = [2*ones(d,1),din];       ot = [1,2];     mout = 2;
else
    dims = [(1+(1:d)'),din];        ot = [2,2];     mout = prod(dims(:,1));
end
dp = struct('in',dep(Zrow,:),'out',true(d,N));
Kop = rand_tensopvar(dims,vars,dom,dp,1,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = Zrow;                 f.varname = {'a';'b';'c'};
f.dims.in = sz(:);          f.dims.out = mout;
f.vars.in = svars;          f.dom.in = dom.in;          f.depmat.in = dep;
f.C = tensopmat(Kop);
f.vars.out = svars;         f.dom.out = dom.out;        f.depmat.out = Kop.deps.out;
f = initialize(f);
end

function v = rand_const_const_partial(svs,m,form)
% Random m x 1 constant function: 'double', or 'polynomial' in the variables svs (cell of pvar)
if strcmp(form,'double')
    v = randn(m,1);
else
    v = polynomial(randn(m,1));
    for i=1:numel(svs)
        v = v + randn(m,1)*svs{i} + randn(m,1)*svs{i}^2;
    end
end
end

function z = pad(z,n)
z(end+1:n) = 0;
end

function z = compress_zeros(z)
% Replace each run of consecutive zeros by a single zero
keep = true(size(z));
keep(2:end) = ~(z(2:end)==0 & z(1:end-1)==0);
z = z(keep);
end

function Y = eval_new_const_partial(f,xvals,sz,allnames)
% Evaluate f at the values xvals.(name) of its states; states without a value get zeros
xv = cell(numel(f.varname),1);
for k=1:numel(xv)
    nm = f.varname{k};
    if isfield(xvals,nm)
        xv{k} = xvals.(nm);
    else
        xv{k} = polynomial(zeros(sz(strcmp(allnames,nm)),1));
    end
end
Y = eval_polyopvar(f,xv);
end

function n = check_partial(f,subnames,cvals,xvals,sz,Zexp,tol,label)
% Substitute the states subnames of f (by name, and by position if the monomial is a single one
% in which they each appear once) by the constants cvals, and compare with f evaluated at the
% constants. Zexp: expected Z of the result, if f has a single monomial (otherwise []).
allnames = {'a','b','c'};
cmat = vertcat(cvals{:});
xs = xvals;
for k=1:numel(subnames),    xs.(subnames{k}) = cvals{k};    end
Yref = eval_polyopvar(f,{xs.a;xs.b;xs.c});
xrem = rmfield(xvals,subnames);
fnew = subs(f,subnames,cmat);
assert_close_const_partial(eval_new_const_partial(fnew,xrem,sz,allnames),Yref,tol,label+" (by name)");
if ~isempty(Zexp)
    % Expected Z: consecutive multipliers from R combined, for an operator with pointwise
    % output (not supported otherwise, so left), in the states of fnew, from which the states
    % that no longer appear have been removed (by name: Zexp is in terms of the states a,b,c)
    op = f.C.ops{1,1};
    if isa(op,'tensopvar') && op.otype(1)==1,   Zexp = compress_zeros(Zexp);    end
    for q=find(Zexp)
        idx = find(strcmp(fnew.varname,allnames{Zexp(q)}));
        if isempty(idx),    error("case_const_partial failed (%s): state %s was removed.",label,allnames{Zexp(q)});    end
        Zexp(q) = idx;
    end
    n = max(size(fnew.Z,2),numel(Zexp));
    if size(fnew.Z,1)~=1 || ~isequal(pad(full(fnew.Z),n),pad(Zexp,n))
        error("case_const_partial failed (%s): Z = [%s], expected [%s].",label,num2str(full(fnew.Z)),num2str(Zexp))
    end
end
n = 1;
Zr = full(f.Z);
idx = cellfun(@(s) find(strcmp(allnames,s)),subnames);
if size(Zr,1)==1 && all(arrayfun(@(i) sum(Zr==i)==1,idx))
    pos = arrayfun(@(i) find(Zr==i),idx);
    fnew = subs(f,pos(:)',cmat);
    assert_close_const_partial(eval_new_const_partial(fnew,xrem,sz,allnames),Yref,tol,label+" (by position)");
    n = 2;
end
end

function g = known_fun(m,vars)
% Degree-0 polyopvar representing the known function m (polynomial) of the variables vars
% (cellstr, in [0,1]): a multiplier from R to L2 on the dummy input 1
M = mat2tensopvar(m,struct('in',{cell(1,0)},'out',{vars}),struct('in',zeros(0,2),'out',repmat([0,1],numel(vars),1)));
g = polyopvar();
g.varname = cell(0,1);      g.Z = sparse(1,0);      g.C = tensopmat(M);
g = initialize(g);
end

function n = check_known(f,subnames,gs,cvals,xvals,sz,tol,label)
% Substitute the states subnames of f by name by the known functions gs (polyopvar, stacked), and
% compare with f evaluated at the constant functions cvals that these represent
allnames = {'a','b','c'};
xs = xvals;
for k=1:numel(subnames),    xs.(subnames{k}) = cvals{k};    end
Yref = eval_polyopvar(f,{xs.a;xs.b;xs.c});
fnew = subs(f,subnames,vertcat(gs{:}));
assert_close_const_partial(eval_new_const_partial(fnew,rmfield(xvals,subnames),sz,allnames),Yref,tol,label);
n = 1;
end

function Y = eval_V(V,vals)
% Evaluate the functional V (intop coefficients) at the states vals.(name); the dummy
% variables of s1 are numbered over the factors, as in apply_op of 'intop'
Y = 0;
Zf = full(V.Z);
for k=1:size(Zf,1)
    op = V.C.ops{1,k};
    if isempty(op),     continue,   end
    idcs = Zf(k,Zf(k,:)~=0);
    v = 1;
    for j=1:numel(idcs)
        xj = vals.(V.varname{idcs(j)});
        if numel(idcs)>1,   xj.varname(strcmp(xj.varname,'s1')) = {['s1_',num2str(j)]};     end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end

function assert_close_const_partial(Ynew,Yref,tol,label)
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("case_const_partial (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("case_const_partial failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function expect_error_const_partial(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("case_const_partial failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("case_const_partial failed: expected an error containing '%s', but none was thrown.",str)
end

function [f,val,Kop] = build_case(N,d,ot,dval)
% Function f(x)(s) = Kop*(x o ... o x), d factors, and val(y) with dval
% degree-one monomials in the state y.
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
dep = struct('in',true(d,N),'out',true(d,N));
Kop = rand_tensopvar(ones(d,2),vars,dom,dep,1,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = ones(1,d);            f.C = tensopmat(Kop);       f.varname = {'x'};
f.dims.in = 1;              f.dims.out = 1;
f.vars.in = svars;          f.vars.out = svars;
f.dom.in = dom.in;          f.dom.out = dom.out;
f.depmat.in = true(1,N);    f.depmat.out = Kop.deps.out;
f = initialize(f);
dep1 = struct('in',true,'out',true);
dep1.in = true(1,N);        dep1.out = true(1,N);
Bop = cell(dval,1);
for k=1:dval
    Bop{k} = tensopmat(rand_tensopvar([1,1],vars,dom,dep1,1,[1,2],2));
end
C = Bop{1};
for k=2:dval,   C = [C,Bop{k}];     end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {'y'};
val.dims.in = 1;            val.dims.out = 1;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function val = build_val_deg2(N)
% val(y)(s) = (B1 y)(s).*(B2 y)(s), pointwise degree-2 monomial
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
dep = struct('in',true(2,N),'out',true(2,N));
Bop = rand_tensopvar(ones(2,2),vars,dom,dep,1,[1,2],2);
val = polyopvar();
val.Z = [1,1];              val.C = tensopmat(Bop);     val.varname = {'y'};
val.dims.in = 1;            val.dims.out = 1;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function p = rand_pol_fctr_function(svs)
p = polynomial(randn);
for i=1:numel(svs)
    p = p + randn*svs{i} + randn*svs{i}^2;
end
end

function expect_error_fctr_function(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("case_fctr_function failed: unexpected error message '%s'.",ME.message)
    end
    return
end
error("case_fctr_function failed: expected an error, but none was thrown.")
end

function f = build_f_multirow(N,sz,kind,mono,outsel)
% Single-row polynomial in vector states a (size sz(1)) and b (size sz(2)),
% both depending on all N spatial variables, with output depending on the
% variables selected by the 1 x N logical OUTSEL:
%   mono 'ab': f = Kop*(a o b);  'a': f = Kop*a;  'b': f = Kop*b
% kind 'pw'/'kron': function with pointwise/Kronecker output;
% kind 'functional' (only 'a'): f = Aop*a, scalar.
ma = sz(1);     mb = sz(2);
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
switch mono
    case 'ab',  Z = [1,2];      din = [ma,mb];
    case 'a',   Z = 1;          din = ma;
    case 'b',   Z = 2;          din = mb;
end
d = numel(Z);
f = polyopvar();
f.Z = Z;                    f.varname = {'a';'b'};
f.dims.in = [ma;mb];
f.vars.in = svars;          f.dom.in = dom.in;
f.depmat.in = true(2,N);
if strcmp(kind,'functional')
    Aop = rand_intop([1,ma],svars,dom.in,ones(1,N),3);
    f.C = tensopmat(Aop);
    f.dims.out = 1;
    f.vars.out = cell(1,0); f.dom.out = zeros(0,2); f.depmat.out = zeros(1,0);
else
    outsel = logical(outsel(:)');
    vars = struct('in',{svars},'out',{svars(outsel)});
    domk = struct('in',dom.in,'out',dom.out(outsel,:));
    dep = struct('in',true(d,N),'out',true(d,sum(outsel)));
    if strcmp(kind,'pw')
        dims = [2*ones(d,1),din(:)];        ot = [1,2];     mout = 2;
    else
        dims = [(1+(1:d)'),din(:)];         ot = [2,2];     mout = prod(dims(:,1));
    end
    Kop = rand_tensopvar(dims,vars,domk,dep,1,ot,2);
    if ~Kop.chkval,     Kop = initialize(Kop);      end
    f.C = tensopmat(Kop);
    f.dims.out = mout;
    f.vars.out = svars(outsel);     f.dom.out = dom.out(outsel,:);      f.depmat.out = Kop.deps.out;
end
f = initialize(f);
end

function val = build_val(N,ma,my,dval)
% val(y) in R^ma, dval degree-one monomials in the vector state y (size my)
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
dep1 = struct('in',true(1,N),'out',true(1,N));
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([ma,my],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {'y'};
val.dims.in = my;           val.dims.out = ma;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function p = rand_pol_multirow(svs,m)
p = polynomial(randn(m,1));
for i=1:numel(svs)
    p = p + randn(m,1)*svs{i} + randn(m,1)*svs{i}^2;
end
end

function v = rand_val_multirow(svs,m,form)
if strcmp(form,'double')
    v = randn(m,1);
else
    v = rand_pol_multirow(svs,m);
end
end

function check(fnew,Yref,vals,sz,tol,label)
% Evaluate fnew at states given as {'name',value,...}; unused states of
% fnew (e.g. the substituted-out 'a') get zeros of the right size
sizes = struct('a',sz(1),'b',sz(2),'y',sz(3));
xv = cell(numel(fnew.varname),1);
for k=1:numel(xv)
    xv{k} = polynomial(zeros(sizes.(fnew.varname{k}),1));
    for i=1:2:numel(vals)
        if strcmp(fnew.varname{k},vals{i}),     xv{k} = vals{i+1};     end
    end
end
Ynew = eval_polyopvar(fnew,xv);
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("case_multirow (%s): reference value is zero; test would be vacuous.",label)
end
if ~isequal(size(Ynew),size(Yref))
    error("case_multirow failed (%s): output size [%s] differs from reference [%s].",label,num2str(size(Ynew)),num2str(size(Yref)))
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("case_multirow failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function f = build_f_name(svars,dom,names,Z,pdegs,q,depmat)
% Functional f(a,b) = sum_k Aop_k*Z_k, for monomials in the rows of Z (padded with zeros),
% states depending on the spatial variables indicated by depmat (default: all). If q>0, the
% kernels are 'dpvar' in q decision variables, with a different, random subset of them in
% each region of integration
N = numel(svars);
nx = numel(names);
if nargin<7
    depmat = true(nx,N);
end
C = [];
for k=1:size(Z,1)
    idcs = Z(k,Z(k,:)~=0);
    Aop = rand_intop([1,1],svars,dom.in,sum(depmat(idcs,:),1),pdegs,q);
    if q>0
        Aop = thin_dvars(Aop);
    end
    Ck = tensopmat(Aop);
    if k==1,    C = Ck;     else,   C = [C,Ck];     end %#ok<AGROW>
end
f = polyopvar();
f.Z = Z;                    f.C = C;                    f.varname = names;
f.dims.in = ones(nx,1);     f.dims.out = 1;
f.vars.in = svars;          f.vars.out = cell(1,0);
f.dom.in = dom.in;          f.dom.out = zeros(0,2);
f.depmat.in = depmat;        f.depmat.out = zeros(1,0);
f = initialize(f);
end

function Aop = thin_dvars(Aop)
% Retain a different random subset of the decision variables in each parameter of Aop
% (the rows of the others, in the blocks of rows of each row of the matrix, are set to zero)
q = numel(Aop.Zd);
m = size(Aop,1);
for k=1:numel(Aop.params)
    P = Aop.params{k};
    if q==0,        continue,       end
    J = sort(randperm(q,ceil(q/2)));
    rows = reshape(((0:m-1)'*(q+1) + [1,1+J]).',[],1);
    P(setdiff(1:size(P,1),rows),:) = 0;
    Aop.params{k} = P;
end
end

function val = build_val_name(svars,dom,yname,dval)
% val(y)(s): dval degree-one monomials in the new state y, random operators
N = numel(svars);
vars = struct('in',{svars},'out',{svars});
dep1 = struct('in',true(1,N),'out',true(1,N));
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([1,1],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok<AGROW>
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {yname};
val.dims.in = 1;            val.dims.out = 1;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function [f,val,Kop] = build_function_case(N,d,ot,dval)
% Function f(x)(s) = Kop*(x o ... o x), d factors, and val(y) with dval degree-one monomials
% in the new state y (as Test_polyopvar_subs/case_fctr_function)
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
dep = struct('in',true(d,N),'out',true(d,N));
Kop = rand_tensopvar(ones(d,2),vars,dom,dep,1,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = ones(1,d);            f.C = tensopmat(Kop);       f.varname = {'x'};
f.dims.in = 1;              f.dims.out = 1;
f.vars.in = svars;          f.vars.out = svars;
f.dom.in = dom.in;          f.dom.out = dom.out;
f.depmat.in = true(1,N);    f.depmat.out = Kop.deps.out;
f = initialize(f);
val = build_val_name(svars,dom,'y',dval);
end

function p = rand_pol_name(svs)
p = polynomial(randn);
for i=1:numel(svs)
    p = p + randn*svs{i} + randn*svs{i}^2;
end
end

function Y = eval_functional_name(V,vals,d0map)
% Evaluate the functional V (polyopvar of intop kernels) at the states vals.(name),
% polynomials in the spatial variables. The dummy variables of s_k are numbered over the
% factors of the monomial that depend on s_k, and named s_k_j (as in apply_op)
Zf = full(V.Z);
fac = cell(size(Zf,1),1);
for k=1:size(Zf,1)
    idcs = Zf(k,Zf(k,:)~=0);
    fac{k} = cellfun(@(nm) vals.(nm),reshape(V.varname(idcs),1,[]),'UniformOutput',false);
end
Y = eval_terms(V,fac,d0map);
end

function Y = eval_functional_fctrs(V,Z,vals,fctrs,rest,d0map)
% Evaluate the functional V, defined by the single monomial Z, with the factors fctrs
% replaced by the polynomials vals{:}, and the others by the states in the struct rest
idcs = Z(Z~=0);
fac = cell(1,numel(idcs));
for j=1:numel(idcs)
    fac{j} = rest.(V.varname{idcs(j)});
end
for k=1:numel(fctrs)
    fac{fctrs(k)} = vals{k};
end
Y = eval_terms(V,{fac},d0map);
end

function Y = eval_terms(V,fac,d0map)
% Sum over the monomials k of Aop_k applied to the product of the factors fac{k} (cell of
% polynomials) of the monomial
vnames = V.vars.in;
dep_state = double(V.depmat.in)>0;
Zf = full(V.Z);
Y = 0;
for k=1:size(Zf,1)
    op = V.C.ops{1,k};
    if isempty(op),     continue,   end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    slot = zeros(1,numel(vnames));
    v = 1;
    for j=1:d
        xj = fac{k}{j};
        for mm=1:numel(vnames)
            if dep_state(idcs(j),mm)
                slot(mm) = slot(mm)+1;
                if d>1
                    isv = strcmp(xj.varname,vnames{mm});
                    xj.varname(isv) = {[vnames{mm},'_',num2str(slot(mm))]};
                end
            end
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
if isa(Y,'dpvar')
    Y = eval_dpvar_at(Y,d0map);
end
Y = double(polynomial(Y));
end

function no_state_left(f,name,label)
% The substituted state may remain declared, but must not appear in Z
[~,k] = ismember(name,f.varname);
if k>0 && any(f.Z(:)==k)
    error("case_name failed (%s): Z still references the substituted state '%s'.",label,name)
end
end

function assert_close_name(Ynew,Yref,tol,label)
if abs(Yref)<1e-6
    error("case_name (%s): reference value %g is zero; test would be vacuous.",label,Yref)
end
if abs(Ynew-Yref)>tol*max(1,abs(Yref))
    error("case_name failed (%s): value %.10g, expected %.10g.",label,Ynew,Yref)
end
end

function expect_error_name(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("case_name failed: unexpected error message '%s'.",ME.message)
    end
    return
end
error("case_name failed: expected an error, but none was thrown.")
end

function [f,facs,cf,Zrow] = build_nested(spec,Id)
% Sum of the monomials, with Z=spec{1} and factors spec{2}{l}, coefficients spec{3}, as a nested operator
Zrow = spec{1};     facs = spec{2};     cf = spec{3};
f = cf(1)*build_f_nested(Zrow,facs{1},Id);
for l=2:numel(facs)
    f = f + cf(l)*build_f_nested(Zrow,facs{l},Id);
end
f = combine_terms(f);
Top = f.C.ops{1,1};
if ~isa(Top,'tensopvar') || ~any(cellfun('isclass',Top.ops(:),'tensopvar'))
    error("case_nested failed: setup: the operator should be nested.")
end
end

function check_nested(g,vals_g,Zrow,facs,cf,subvals,tol,label) %#ok
% Evaluate g at the values of its remaining states (vals_g: struct by name; substituted-out states that remain
% declared get the value 0, and must not appear in Z) and compare with the product of the factors, with x and y
% replaced by the substituted values subvals.x, subvals.y, if given, else by their value in vals_g
Yref = 0;
for l=1:numel(facs)
    Yref = Yref + cf(l)*prod_facs(facs{l},subvals,vals_g);
end
xv = cell(numel(g.varname),1);
for kk=1:numel(xv)
    nm = g.varname{kk};
    if isfield(vals_g,nm)
        xv{kk} = vals_g.(nm);
    else
        if any(g.Z(:)==kk)
            error("case_nested failed (%s): the substituted state '%s' still appears in Z.",label,nm)
        end
        xv{kk} = polynomial(0);
    end
end
Yg = eval_polyopvar(g,xv);
assert_close_nested(Yg,Yref,tol,label);
end

function Y = prod_facs(facs,subvals,vals_g)
% Product of the values of the factors, as polynomials: x -> its value, y -> its value, polynomial -> itself
Y = 1;
for j=1:numel(facs)
    if ischar(facs{j})
        if isfield(subvals,facs{j}) && ~isempty(subvals.(facs{j}))
            Y = Y*subvals.(facs{j});
        else
            Y = Y*vals_g.(facs{j});
        end
    else
        Y = Y*facs{j};
    end
end
end

function assert_close_nested(Y,Yref,tol,label)
% Output variables of the nested operator are named s_1,...: renamed back
Y = polynomial(Y);      Yref = polynomial(Yref);
Y.varname = regexprep(Y.varname,'_1$','');
D = Y-Yref;
if ~isequal(size(D),size(Yref))
    error("case_nested failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
err = full(max([0;abs(nonzeros(D.coefficient))]));
if err>tol*max(1,full(max(abs(nonzeros(Yref.coefficient)))))
    error("case_nested failed %s: residual %.3e.",label,err)
end
if isempty(nonzeros(Yref.coefficient)) || max(abs(nonzeros(Yref.coefficient)))<1e-6
    error("case_nested %s: reference value is zero; test would be vacuous.",label)
end
end

function p = rand_const_nested(svs)
% Random polynomial in the spatial variables svs (cell or array of pvars), degree 2 in each
p = polynomial(randn);
for i=1:numel(svs)
    si = svs(i);
    if iscell(si),  si = si{1};     end
    p = p + randn*si + randn*si^2;
end
end

function val = build_val_nested(name,dval)
% polyopvar val(w)(s) = sum_k B_k w(s), dval random degree-one monomials in the new state w, in s1,s2
svars = {'s1','s2'};
dom = struct('in',[0,1;0,1],'out',[0,1;0,1]);
vars = struct('in',{svars},'out',{svars});
dep1 = struct('in',true(1,2),'out',true(1,2));
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([1,1],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {name};
val.dims.in = 1;            val.dims.out = 1;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,2);  val.depmat.out = ones(1,2);
val = initialize(val);
end

function f = build_f_nested(Zrow,facs,Id)
% polyopvar with monomial Z=Zrow in the states x,y, acted on by the product of the factors facs: the
% identity operators of the states for 'x','y', multipliers from R for polynomials
K = [];
for j=1:numel(facs)
    if ischar(facs{j})
        Fj = Id.(facs{j});
    else
        vars = cell(1,0);
        if ismember('s1',facs{j}.varname),  vars = [vars,{'s1'}];   end
        if ismember('s2',facs{j}.varname),  vars = [vars,{'s2'}];   end
        Fj = mat2tensopvar(facs{j},struct('in',{cell(1,0)},'out',{vars}),struct('in',zeros(0,2),'out',repmat([0,1],numel(vars),1)));
    end
    if isempty(K),  K = Fj;
    else,           K = otimes(K,Fj,[1,2]);     end
end
f = polyopvar();
f.varname = {'x';'y'};      f.Z = sparse(Zrow);     f.C = tensopmat(K);
f.dims.in = [1;1];
f.vars.in = {'s1','s2'};    f.dom.in = [0,1;0,1];   f.depmat.in = true(2,2);
f = initialize(f);
end

function xv = state_values(P,PIE,X0,U0)
% Values of the states of P, by name: PIE states x0, PDE states u0 (states of P that are not used in P
% are given a value of the right size)
xv = cell(numel(P.varname),1);
for k=1:numel(xv)
    [tf,jx] = ismember(P.varname{k},PIE.varname_pie);
    if tf
        xv{k} = X0{jx};
    else
        [~,jx] = ismember(P.varname{k},PIE.varname_pde);
        xv{k} = U0{jx};
    end
end
end

function Y = eval_V_ode_rows(V,vals)
% Evaluate a 1x1 polyopvar whose coefficients are integral operators acting on d-fold tensor products
% of the states (as in Test_piesos_LF): the dummy variables of s_k are numbered over the factors that
% depend on s_k, as s_k_j, if the monomial has more than 1 factor
vnames = V.vars.in;
dep_state = double(V.depmat.in)>0;
Y = 0;
Zf = full(V.Z);
for k=1:size(V.Z,1)
    op = V.C.ops{1,k};
    if isempty(op),     continue,   end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    slot = zeros(1,numel(vnames));
    v = 1;
    for j=1:d
        xj = vals.(V.varname{idcs(j)});
        for mm=1:numel(vnames)
            if dep_state(idcs(j),mm)
                slot(mm) = slot(mm)+1;
                if d>1
                    isv = strcmp(xj.varname,vnames{mm});
                    xj.varname(isv) = {[vnames{mm},'_',num2str(slot(mm))]};
                end
            end
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end

function check_equal(Y,Yref,tol,label)
if ~isequal(size(Y),size(Yref))
    error("case_ode_rows failed (%s): size [%s] differs from reference [%s].",label,num2str(size(Y)),num2str(size(Yref)))
end
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("case_ode_rows (%s): reference is zero; test would be vacuous.",label)
end
err = polynomial(Y-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol*max(1,max(abs(full(Yc(:)))))
    error("case_ode_rows failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function check_subs(f,fnew,val,svs,tol,label)
a = rand_pol_var_function(svs);      b = rand_pol_var_function(svs);      y = rand_pol_var_function(svs);
vy = eval_polyopvar(val,{y});
Yref = eval_polyopvar(f,{vy;b});
xv = cell(numel(fnew.varname),1);
for k=1:numel(xv)
    if strcmp(fnew.varname{k},'y'),     xv{k} = y;
    elseif strcmp(fnew.varname{k},'b'), xv{k} = b;
    elseif strcmp(fnew.varname{k},'a')
        % Substituted-out state may remain declared, but must not appear in Z
        if any(fnew.Z(:)==k)
            error("case_var_function failed (%s): fnew.Z still references the substituted state 'a'.",label)
        end
        xv{k} = polynomial(0);
    else
        error("case_var_function failed (%s): unexpected state '%s' in result.",label,fnew.varname{k})
    end
end
Ynew = eval_polyopvar(fnew,xv);
assert_close_var_function(Ynew,Yref,tol,label)
end

function assert_close_var_function(Ynew,Yref,tol,label)
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("case_var_function (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("case_var_function failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function [f,val,svars,dom] = build_f_var_function(N,kind,a_const,dval)
% f(a,b) with monomial a o b; if a_const, state a does not depend on any
% spatial variable. val(y): dval degree-one monomials in a new state y.
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
depin = [true(1,N)&~a_const; true(1,N)];
f = polyopvar();
f.Z = [1,2];                f.varname = {'a';'b'};
f.dims.in = [1;1];          f.dims.out = 1;
f.vars.in = svars;          f.dom.in = dom.in;
f.depmat.in = depin;
switch kind
    case {'function_pw','function_kron'}
        if strcmp(kind,'function_pw'),  ot = [1,2];     else,   ot = [2,2];     end
        dep = struct('in',depin,'out',true(2,N));
        Kop = rand_tensopvar(ones(2,2),vars,dom,dep,1,ot,2);
        if ~Kop.chkval,     Kop = initialize(Kop);      end
        f.C = tensopmat(Kop);
        f.vars.out = svars;     f.dom.out = dom.out;    f.depmat.out = Kop.deps.out;
    case 'functional'
        % degree-one functional only: f(a)=Aop*a
        f.Z = 1;
        Aop = rand_intop([1,1],svars,dom.in,double(depin(1,:)),3);
        f.C = tensopmat(Aop);
        f.vars.out = cell(1,0); f.dom.out = zeros(0,2); f.depmat.out = zeros(1,0);
end
f = initialize(f);
% val
dep1 = struct('in',true(1,N),'out',depin(1,:));
if a_const
    val = [];   return
end
dep1.out = true(1,N);
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([1,1],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {'y'};
val.dims.in = 1;            val.dims.out = 1;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function p = rand_pol_var_function(svs)
p = polynomial(randn);
for i=1:numel(svs)
    p = p + randn*svs{i} + randn*svs{i}^2;
end
end

function f = build_const_f(N,d,ot)
% f(a)(s) = Kop*(a o ... o a) (d factors), a not depending on spatial variables
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',zeros(0,2),'out',repmat([0,1],N,1));
vars = struct('in',{cell(1,0)},'out',{svars});
dep = struct('in',false(d,0),'out',true(d,N));
Kop = rand_tensopvar(ones(d,2),vars,dom,dep,1,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = ones(1,d);            f.varname = {'a'};
f.dims.in = 1;              f.dims.out = 1;
f.vars.in = cell(1,0);      f.dom.in = zeros(0,2);      f.depmat.in = false(1,0);
f.C = tensopmat(Kop);
f.vars.out = svars;         f.dom.out = dom.out;        f.depmat.out = Kop.deps.out;
f = initialize(f);
end

function [f,Kop] = build_f_vector_states(N,sz,kind,mono)
% Polynomial in vector states a (size sz(1)) and b (size sz(2)), both
% depending on all N spatial variables:
%   mono 'ab': f = Kop*(a o b);   'aa': f = Kop*(a o a);   'a': f = Kop*a
% kind 'pw'/'kron': function with pointwise/Kronecker output (output size 2
% resp. 2*2 for the pointwise/Kronecker product of two factors; 2 for
% one factor); kind 'functional' (only mono 'a'): f = Aop*a, scalar.
ma = sz(1);     mb = sz(2);
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
switch mono
    case 'ab',  Z = [1,2];      din = [ma,mb];
    case 'aa',  Z = [1,1];      din = [ma,ma];
    case 'a',   Z = 1;          din = ma;
end
d = numel(Z);
f = polyopvar();
f.Z = Z;                    f.varname = {'a';'b'};
f.dims.in = [ma;mb];
f.vars.in = svars;          f.dom.in = dom.in;
f.depmat.in = true(2,N);
Kop = [];
if strcmp(kind,'functional')
    Aop = rand_intop([1,ma],svars,dom.in,ones(1,N),3);
    f.C = tensopmat(Aop);
    f.dims.out = 1;
    f.vars.out = cell(1,0); f.dom.out = zeros(0,2); f.depmat.out = zeros(1,0);
else
    dep = struct('in',true(d,N),'out',true(d,N));
    if strcmp(kind,'pw')
        dims = [2*ones(d,1),din(:)];        ot = [1,2];     mout = 2;
    else
        dims = [(1+(1:d)'),din(:)];         ot = [2,2];     mout = prod(dims(:,1));
    end
    Kop = rand_tensopvar(dims,vars,dom,dep,1,ot,2);
    if ~Kop.chkval,     Kop = initialize(Kop);      end
    f.C = tensopmat(Kop);
    f.dims.out = mout;
    f.vars.out = svars;     f.dom.out = dom.out;    f.depmat.out = Kop.deps.out;
end
f = initialize(f);
end

function val = build_val_vector_states(N,ma,my,dval)
% val(y) in R^ma, dval degree-one monomials in the vector state y (size my)
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
dep1 = struct('in',true(1,N),'out',true(1,N));
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([ma,my],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {'y'};
val.dims.in = my;           val.dims.out = ma;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function p = rand_pol_vector_states(svs,m)
p = polynomial(randn(m,1));
for i=1:numel(svs)
    p = p + randn(m,1)*svs{i} + randn(m,1)*svs{i}^2;
end
end

function v = rand_val_vector_states(svs,m,form)
if strcmp(form,'double')
    v = randn(m,1);
else
    v = rand_pol_vector_states(svs,m);
end
end

function check_vector_states(fnew,Yref,vals,sz,tol,label)
% Evaluate fnew at states given as {'name',value,...}; unused states of
% fnew (e.g. the substituted-out 'a') get zeros of the right size
sizes = struct('a',sz(1),'b',sz(2),'y',sz(3));
xv = cell(numel(fnew.varname),1);
for k=1:numel(xv)
    xv{k} = polynomial(zeros(sizes.(fnew.varname{k}),1));
    for i=1:2:numel(vals)
        if strcmp(fnew.varname{k},vals{i}),     xv{k} = vals{i+1};     end
    end
end
Ynew = eval_polyopvar(fnew,xv);
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("case_vector_states (%s): reference value is zero; test would be vacuous.",label)
end
if ~isequal(size(Ynew),size(Yref))
    error("case_vector_states failed (%s): output size [%s] differs from reference [%s].",label,num2str(size(Ynew)),num2str(size(Yref)))
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("case_vector_states failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function expect_error_vector_states(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("case_vector_states failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("case_vector_states failed: expected an error containing '%s', but none was thrown.",str)
end
