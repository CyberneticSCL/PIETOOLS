function R = pf_fixed(tag,ops,optl,tier)
% R = PF_FIXED(TAG,OPS,OPTL,TIER) certifies the positivity of fixed operators:
% for each named 'copvar' X in the struct array OPS (fields name, X) and
% each option struct in the cell OPTL, the constraint X >= 0 by
% LPI_INEQ_SOP (the slack reader sizes the certificate from the support of
% X unless the options say otherwise), solved by MOSEK and classified:
% 'optimal' = X is in the parameterized cone at that sizing; 'infeasible'
% = a MOSEK certificate that it is not; 'unknown' = no decision. Prints
% the reader's lift and weight degrees and the Gram sizes. Results in
% results/fixed_<TAG>.mat. SCRATCH, MMP, 11/16/2026.
if nargin<4 || isempty(tier),   tier = 'light';     end
S = pf_setup(false);
out = sprintf('%s/results/fixed_%s.mat',S,tag);
R = struct('name',{},'opts',{},'status',{},'ndv',{},'t',{},'D',{},'w',{},'N',{},'numerr',{},'feasratio',{},'relres',{},'eta',{},'err',{});
for i = 1:numel(ops)
    X = ops(i).X;
    for j = 1:numel(optl)
        o = optl{j};    olab = jsonencode(o);
        r = struct('name',ops(i).name,'opts',olab,'status','','ndv',NaN,'t',NaN,'D','','w','','N','', ...
                   'numerr',NaN,'feasratio',NaN,'relres',NaN,'eta',NaN,'err','');
        t0 = tic;
        try
            st = lpisettings_sop(tier);     st.sos_opts.solver = 'mosek';
            prog = lpiprogram(X.vars(:),[],X.dom,st);
            [prog,~,info] = lpi_ineq_sop(prog,X,o);
            r.ndv = size(prog.decvartable,1);
            if isfield(info,'degrees') && isfield(info.degrees,'D')
                r.D = mat2str(info.degrees.D);   r.w = mat2str(info.degrees.w);
            end
            if isfield(info,'N'),   r.N = mat2str(info.N);    end
            prog = lpisolve(prog,st.sos_opts);
            cl = lpi_classify_sop(prog,[]);
            r.t = toc(t0);  r.status = cl.status;   r.numerr = cl.numerr;
            r.feasratio = cl.feasratio;     r.relres = cl.relres;  r.eta = cl.eta;
        catch ME
            r.err = strrep(ME.message,newline,' ');     r.t = toc(t0);
        end
        R(end+1) = r;                                                       %#ok<AGROW>
        fprintf('@@ FIX %s opts=%s %s ndv=%d t=%.0f D=%s w=%s N=%s numerr=%g feas=%.3g relres=%.1e eta=%.1e %s\n', ...
                r.name,olab,r.status,r.ndv,r.t,r.D,r.w,r.N,r.numerr,r.feasratio,r.relres,r.eta,r.err);
        save(out,'R');
    end
end
fprintf('## FIX DONE %s\n',tag);
end
