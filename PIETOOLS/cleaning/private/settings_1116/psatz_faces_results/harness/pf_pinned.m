function R = pf_pinned(N,bc,sets,dws,rungs,tag,route)
% R = PF_PINNED(N,BC,SETS,DWS,RUNGS,TAG) the pinned-P certificates of the
% plan (C0 in one variable, C1 in two): x_t = sum_i x_{s_i s_i} + lam x on
% [0,1]^N with the boundary conditions BC (HEATND_PIE), P = I, so the
% Lyapunov derivative -(T'A0 + A0'T) - 2 lam T'T must be certified >= 0 by
% LPI_INEQ_SOP with the term set alone as the unknown. For every set in
% SETS (PF_TERMS names), every dw in DWS and every rung f in RUNGS the
% program is built (lift from the reader, weight w = rule + dw) and solved
% with HEATND_SOLVE (tight): st +1 certified, -1 Farkas, 0 uncertain.
% Results are appended to results/pinned_<TAG>.mat after every solve and
% printed as '@@ PIN' lines. SCRATCH, MMP, 11/15/2026.
if nargin<7 || isempty(route),  route = false;  end
S = pf_setup(route);
pie = heatNd_pie(N,0,bc);   lam1 = pie.exact.lambda1;
T = pie.T;  A0 = pie.A0;    vars = pie.vars;    dom = pie.dom;
X0 = -(T'*A0 + A0'*T);      Y = -2*(T'*T);
out = sprintf('%s/results/pinned_%s.mat',S,tag);
R = struct('N',{},'bc',{},'set',{},'o',{},'pd',{},'dw',{},'D',{},'w',{},'f',{},'st',{},'why',{}, ...
           'rel_b',{},'psd',{},'viol',{},'Ngram',{},'ndec',{},'m',{},'nx',{},'tb',{},'ts',{},'iter',{},'thr',{},'err',{});
for is = 1:numel(sets)
    [rows,offs,pd] = pf_terms(sets{is},N,[]);
    o = 1;  if ~isempty(offs),  o = offs(1);    end
    for dw = dws
        for f = rungs
            r = struct('N',N,'bc',{bc},'set',sets{is},'o',o,'pd',pd,'dw',dw,'D',NaN,'w',NaN,'f',f,'st',NaN,'why','', ...
                       'rel_b',NaN,'psd',NaN,'viol',NaN,'Ngram',[],'ndec',NaN,'m',NaN,'nx',NaN,'tb',NaN,'ts',NaN,'iter',NaN,'thr',NaN,'err','');
            try
                opts = struct('psatz',0,'dw',dw);
                if ~isempty(rows)
                    opts.psatz_rows = rows;     opts.psatz_rows_offset = offs;  opts.psatz_rows_perdir = pd;
                end
                prog = lpiprogram(polynomial(vars(:)),[],dom);
                t0 = tic;
                [prog,~,info] = lpi_ineq_sop(prog,X0 + f*lam1*Y,opts);
                r.tb = toc(t0);
                r.Ngram = info.N;   r.ndec = info.ndec;
                if isfield(info.degrees,'D'),   r.D = mat2str(info.degrees.D);  r.w = mat2str(info.degrees.w);   end
                t0 = tic;
                v = heatNd_solve(prog,struct('tight',true));
                r.ts = toc(t0);
                r.st = v.st;    r.why = v.why;  r.rel_b = v.rel_b;  r.psd = v.psd_relmin;  r.viol = v.cert_viol;
                r.m = v.m;      r.nx = v.nx;    r.iter = v.iter;
                if isfield(v,'mosek_threads'),  r.thr = v.mosek_threads;    end
            catch ME
                r.err = strrep(ME.message,newline,' ');
            end
            R(end+1) = r;                                                   %#ok<AGROW>
            fprintf('@@ PIN N=%d set=%s o=%d pd=%d dw=%d D=%s w=%s f=%.4f st=%+d relb=%.1e psd=%.1e viol=%.1e Ngram=%s ndec=%d m=%d nx=%d tb=%.0f ts=%.0f it=%d thr=%d %s %s\n', ...
                    N,sets{is},o,pd,dw,num2str(r.D),num2str(r.w),f,r.st,r.rel_b,r.psd,r.viol,mat2str(r.Ngram),r.ndec,r.m,r.nx,r.tb,r.ts,r.iter,r.thr,r.why,r.err);
            save(out,'R');
        end
    end
end
fprintf('## PIN DONE %s\n',tag);
end
