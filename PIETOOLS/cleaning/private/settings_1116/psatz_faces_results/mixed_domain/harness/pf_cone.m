function R = pf_cone(tag,ops,variants)
% R = PF_CONE(TAG,OPS,VARIANTS) tests whether fixed operators lie in the
% cone of a storage parameterization: for each named 'copvar' X in OPS and
% each variant {name, sepv, pos} in VARIANTS, declares P by the storage
% builder (name 'def': POSCOPVAR with DEG omitted, the uniform patterns of
% the tiers, degrees prog.settings.lpi.pos = pos; 'sep': PF_POS_SEP with the
% per-direction split flags sepv, every basis operator, degrees pos),
% imposes X = P (LPI_EQ_SOP) and solves: 'optimal' = X is in the cone,
% 'infeasible' = a MOSEK certificate that it is not. Results in
% results/cone_<TAG>.mat. SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
out = sprintf('%s/results/cone_%s.mat',S,tag);
R = struct('name',{},'variant',{},'status',{},'ndv',{},'t',{},'numerr',{},'feasratio',{},'relres',{},'eta',{},'err',{});
for i = 1:numel(ops)
    X = ops(i).X;
    [sp,dm] = copvar_space_list(X,'out');
    dom = struct('vars',{X.vars},'dom',X.dom);
    for j = 1:numel(variants)
        v = variants{j};    vlab = sprintf('%s_%s_%s',v{1},char('0'+v{2}),v{4});
        r = struct('name',ops(i).name,'variant',vlab,'status','','ndv',NaN,'t',NaN,'numerr',NaN,'feasratio',NaN,'relres',NaN,'eta',NaN,'err','');
        t0 = tic;
        try
            st = lpisettings_sop('light');  st.sos_opts.solver = 'mosek';
            st.lpi.pos = v{3};
            prog = lpiprogram(X.vars(:),[],X.dom,st);
            if strcmp(v{1},'def')
                [prog,P] = poscopvar(prog,dm,sp,dom);
            else
                setappdata(0,'PF_POS',v{3});
                [prog,P] = pf_pos_sep(prog,X,'out',dom,v{2});
                rmappdata(0,'PF_POS');
            end
            prog = lpi_eq_sop(prog,X - P);
            r.ndv = size(prog.decvartable,1);
            prog = lpisolve(prog,st.sos_opts);
            cl = lpi_classify_sop(prog,[]);
            r.t = toc(t0);  r.status = cl.status;   r.numerr = cl.numerr;
            r.feasratio = cl.feasratio;     r.relres = cl.relres;  r.eta = cl.eta;
        catch ME
            r.err = strrep(ME.message,newline,' ');     r.t = toc(t0);
        end
        R(end+1) = r;                                                       %#ok<AGROW>
        fprintf('@@ CONE %s variant=%s %s ndv=%d t=%.0f numerr=%g feas=%.3g relres=%.1e eta=%.1e %s\n', ...
                r.name,vlab,r.status,r.ndv,r.t,r.numerr,r.feasratio,r.relres,r.eta,r.err);
        save(out,'R');
    end
end
fprintf('## CONE DONE %s\n',tag);
end
