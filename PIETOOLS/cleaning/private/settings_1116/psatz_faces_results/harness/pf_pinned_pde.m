function R = pf_pinned_pde(name,sets,dws,lams,tag)
% R = PF_PINNED_PDE(NAME,SETS,DWS,LAMS,TAG) pinned-P certificates on the
% converter-built plants of E2 (evidence note): 'L1' reaction lam (1 - 2 s1) x,
% 'L2' diffusion (0.05 + s1 (1 - s1)) x_{s1 s1} + x_{s2 s2} + lam x, 'ODE'
% x_t = Lap x + lam x + X, X_t = -2 X + int int x; Dirichlet on every edge
% of [0,1]^2. The PIE is affine in lam: A = A0 + lam A1 from the PIEs at
% lam = 0 and 1. For P = I the derivative -(T'A + A'T) is certified by
% LPI_INEQ_SOP with the term set alone as the unknown, at the absolute
% rungs LAMS, verdicts by HEATND_SOLVE (tight). SCRATCH, MMP, 11/15/2026.
S = pf_setup();
pvar s1 s2 t
PIE = cell(1,2);
for i = 1:2
    lam = i-1;
    switch name
        case 'L1'
            x = pde_var(1,[s1;s2],[0,1;0,1]);
            PDE = [diff(x,t)==diff(x,s1,2)+diff(x,s2,2)+lam*(1-2*s1)*x; bc4(x,s1,s2)];
        case 'L2'
            x = pde_var(1,[s1;s2],[0,1;0,1]);
            PDE = [diff(x,t)==(0.05+s1*(1-s1))*diff(x,s1,2)+diff(x,s2,2)+lam*x; bc4(x,s1,s2)];
        case 'ODE'
            X = pde_var(1);     x = pde_var(1,[s1;s2],[0,1;0,1]);
            PDE = [diff(X,t)==-2*X+int(x,[s1;s2],[0,1;0,1]);
                   diff(x,t)==diff(x,s1,2)+diff(x,s2,2)+lam*x+X; bc4(x,s1,s2)];
        otherwise,  error('pf_pinned_pde:name','unknown plant %s',name);
    end
    PIE{i} = initialize(convert(PDE,{'pie','silent'}));
end
T = PIE{1}.T;   A0 = PIE{1}.A;  A1 = PIE{2}.A - PIE{1}.A;
X0 = -(T'*A0 + A0'*T);      X1 = -(T'*A1 + A1'*T);
vars = T.vars;  dom = T.dom;    N = numel(vars);
out = sprintf('%s/results/pinpde_%s.mat',S,tag);
R = struct('plant',{},'set',{},'o',{},'pd',{},'dw',{},'D',{},'w',{},'lam',{},'st',{},'why',{}, ...
           'rel_b',{},'psd',{},'viol',{},'Ngram',{},'ndec',{},'m',{},'nx',{},'tb',{},'ts',{},'iter',{},'err',{});
for is = 1:numel(sets)
    [rows,offs,pd] = pf_terms(sets{is},N,[]);
    o = 1;  if ~isempty(offs),  o = offs(1);    end
    for dw = dws
        for lam = lams
            r = struct('plant',name,'set',sets{is},'o',o,'pd',pd,'dw',dw,'D','','w','','lam',lam,'st',NaN,'why','', ...
                       'rel_b',NaN,'psd',NaN,'viol',NaN,'Ngram',[],'ndec',NaN,'m',NaN,'nx',NaN,'tb',NaN,'ts',NaN,'iter',NaN,'err','');
            try
                opts = struct('psatz',0,'dw',dw);
                if ~isempty(rows)
                    opts.psatz_rows = rows;     opts.psatz_rows_offset = offs;  opts.psatz_rows_perdir = pd;
                end
                prog = lpiprogram(polynomial(vars(:)),[],dom);
                t0 = tic;   [prog,~,info] = lpi_ineq_sop(prog,X0 + lam*X1,opts);   r.tb = toc(t0);
                r.Ngram = info.N;   r.ndec = info.ndec;
                if isfield(info.degrees,'D'),   r.D = mat2str(info.degrees.D);  r.w = mat2str(info.degrees.w);   end
                t0 = tic;   v = heatNd_solve(prog,struct('tight',true));    r.ts = toc(t0);
                r.st = v.st;    r.why = v.why;  r.rel_b = v.rel_b;  r.psd = v.psd_relmin;  r.viol = v.cert_viol;
                r.m = v.m;      r.nx = v.nx;    r.iter = v.iter;
            catch ME
                r.err = strrep(ME.message,newline,' ');
            end
            R(end+1) = r;                                                   %#ok<AGROW>
            fprintf('@@ PDE %s set=%s o=%d pd=%d dw=%d D=%s w=%s lam=%g st=%+d relb=%.1e psd=%.1e viol=%.1e Ngram=%s ndec=%d m=%d nx=%d tb=%.0f ts=%.0f %s %s\n', ...
                    name,sets{is},o,pd,dw,r.D,r.w,lam,r.st,r.rel_b,r.psd,r.viol,mat2str(r.Ngram),r.ndec,r.m,r.nx,r.tb,r.ts,r.why,r.err);
            save(out,'R');
        end
    end
end
fprintf('## PDE DONE %s %s\n',name,tag);
end

function B = bc4(x,s1,s2)
% Dirichlet on the four edges of [0,1]^2
B = [subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
end
