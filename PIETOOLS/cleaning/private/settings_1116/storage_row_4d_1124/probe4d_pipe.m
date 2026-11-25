function probe4d_pipe(label)
% PROBE4D_PIPE(LABEL) the container pipeline on the heat equation
% x_t = sum_k x_{s_k s_k} + x with Dirichlet conditions, in 3 and 4
% variables: declaration (pde_var), conversion (convert_PIETOOLS_PDE),
% exec_pie_sop, the free operator of the Q form (lpivar_cdopvar), the
% storage of the light tier (pos_row of exec_tools_sop), the composition
% D = T'PA + A'PT, and lpi_ineq_sop on -D as PIETOOLS_stability_sop 'P'
% poses it at its first attempt. In 3 variables PIESIM is called on the
% PDE to record its refusal. Times, decision-variable counts and error
% messages go to the diary LABEL_diary.log.
here = fileparts(mfilename('fullpath'));
diary(fullfile(here,[label '_diary.log']));
root = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
if isempty(which('poscopvar')),  run(fullfile(root,'pietools_path_update.m'));  end
w = strrep(which('poscopvar'),'\','/');
assert(strcmpi(w,[root '/lpi_programming_sopvar/poscopvar.m']),'poscopvar resolves to %s',w);
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
fprintf('--- %s (%s)\n',label,datestr(now,'yyyy-mm-dd HH:MM:SS'));
for n = [3 4]
    nm = arrayfun(@(k) sprintf('s%d',k),1:n,'UniformOutput',false);
    s = polynomial(nm(:));
    clear stateNameGenerator
    x = pde_var(s,repmat([0 1],n,1));
    rhs = diff(x,s(1),2);
    for k = 2:n,    rhs = rhs + diff(x,s(k),2);    end
    rhs = rhs + x;
    PDE = diff(x,'t')==rhs;
    for k = 1:n,    PDE = [PDE; subs(x,s(k),0)==0; subs(x,s(k),1)==0];    end
    t0 = tic;
    try
        PIE = convert_PIETOOLS_PDE(PDE);
        fprintf('n=%d convert %.1f s; PIE.dim %d; T %s %s; A %s\n',n,toc(t0), ...
                PIE.dim,class(PIE.T),mat2str(size(PIE.T)),class(PIE.A));
    catch ME
        fprintf('n=%d convert ERR %.1f s: %s\n',n,toc(t0),ME.message);   continue
    end
    if n==3
        t0 = tic;
        try
            PIESIM(PDE,struct('plot','no'));
            fprintf('n=%d PIESIM ran %.1f s\n',n,toc(t0));
        catch ME
            fprintf('n=%d PIESIM ERR %.1f s: %s\n',n,toc(t0),ME.message);
        end
    end
    PIE = exec_pie_sop(PIE);
    st = exec_sop_settings(lpisettings_sop('light'));
    tl = exec_tools_sop(PIE,st);    st.sop = tl.sop;
    prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);
    Tm = PIE.T;     Am = PIE.A;
    % The free operator of the Q form: deg_Q of this dimension
    q = st.sop.deg_Q{min(n,end)};
    t0 = tic;
    try
        [p2,~] = lpivar_cdopvar(prog,Tm,q);
        fprintf('n=%d light lpivar_cdopvar(T,%s) %.1f s, %d dvars\n',n,mat2str(q),toc(t0),size(p2.decvartable,1));
    catch ME
        fprintf('n=%d light lpivar_cdopvar ERR %.1f s: %s\n',n,toc(t0),ME.message);
    end
    % The storage of the stability executive, form 'P', and its margin
    t0 = tic;
    [prog,Pm] = tl.storage(prog,Tm,'out');
    Pm = Pm + tl.eye(1,tl.eppos);
    fprintf('n=%d light storage %.1f s, %d dvars\n',n,toc(t0),size(prog.decvartable,1));
    t0 = tic;
    Dm = Tm'*Pm*Am + Am'*Pm*Tm;
    fprintf('n=%d light D = T''PA + A''PT %.1f s\n',n,toc(t0));
    % The slack as the first attempt of the executive poses it: dw =
    % dw_feas of this dimension, the cap Jmax(T) + cap_delta
    o = struct('like',Pm,'dw',st.sop.dw_feas,'margin',st.sop.margin,'cap',tl.capfor(Dm));
    fprintf('n=%d light lpi_ineq_sop: dw %d, cap %s\n',n,o.dw,mat2str(o.cap));
    t0 = tic;
    try
        [p3,~,li] = lpi_ineq_sop(prog,-Dm,o);
        fprintf('n=%d light slack built %.1f s, %d dvars in the program, degrees %s\n',n,toc(t0), ...
                size(p3.decvartable,1),mat2str(li.degrees));
    catch ME
        fprintf('n=%d light lpi_ineq_sop ERR %.1f s: %s\n',n,toc(t0),ME.message);
    end
end
fprintf('--- done (%s)\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
diary off
end
