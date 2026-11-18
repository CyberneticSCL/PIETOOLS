% PF_PROBE  Phase 0 of the psatz_faces run: the instrument checks and one
% unit-cost probe. (a) 1-D: the row 'box' at offset 1 reproduces the
% repo's Markov-Lukacs pair ('auto'); (b) the per-direction rule changes
% the Gram sizes in 2-D as intended; (c) one 2-D pinned-P solve with the
% faces at dw 0 for the unit cost. SCRATCH, MMP, 11/15/2026.
S = pf_setup();
%% (a) 1-D: box rows at o = 1 against the repo pair
pie = heatNd_pie(1,0,{'DD'});   lam1 = pie.exact.lambda1;
T = pie.T;  A0 = pie.A0;    X = -(T'*A0 + A0'*T) - 2*0.5*lam1*(T'*T);
vars = pie.vars;    dom = pie.dom;
p0 = lpiprogram(polynomial(vars(:)),[],dom);
[p1,~,i1] = lpi_ineq_sop(p0,X,struct('psatz','auto','dw',1));
[rows,offs,pd] = pf_terms('box',1,1);
[p2,~,i2] = lpi_ineq_sop(p0,X,struct('psatz',0,'psatz_rows',rows,'psatz_rows_offset',offs,'psatz_rows_perdir',pd,'dw',1));
D1 = heatNd_sdp(p1);    D2 = heatNd_sdp(p2);
same = isequal(i1.N,i2.N) && isequal(size(D1.At),size(D2.At)) && nnz(D1.At)==nnz(D2.At) ...
       && norm(sort(nonzeros(D1.At))-sort(nonzeros(D2.At)))<1e-12 && norm(D1.b-D2.b)<1e-12;
fprintf('## (a) 1-D pair vs box rows: N %s vs %s, m %d vs %d, nx %d vs %d, SAME %d\n', ...
        mat2str(i1.N),mat2str(i2.N),size(D1.At,2),size(D2.At,2),size(D1.At,1),size(D2.At,1),same);
%% (b) 2-D Gram sizes per set and offset rule (estimate only: build without solve)
pie = heatNd_pie(2,0,{'DD','DD'});  lam1 = pie.exact.lambda1;
T = pie.T;  A0 = pie.A0;    X = -(T'*A0 + A0'*T) - 2*0.5*lam1*(T'*T);
vars = pie.vars;    dom = pie.dom;
p0 = lpiprogram(polynomial(vars(:)),[],dom);
for s = {'none','box','box0','product','faces','faces0','facesS','corners'}
    for dw = [0 1]
        [rows,offs,pd] = pf_terms(s{1},2,[]);
        opts = struct('psatz',0,'dw',dw);
        if ~isempty(rows),  opts.psatz_rows = rows;  opts.psatz_rows_offset = offs;  opts.psatz_rows_perdir = pd;  end
        t0 = tic;   [~,~,info] = lpi_ineq_sop(p0,X,opts);  tb = toc(t0);
        ints = '';
        for c = 1:numel(info.direct.terms)
            tc = info.direct.terms(c);
            if isfield(tc,'row'),   ints = [ints sprintf(' [%s]',num2str(tc.row))];    end %#ok<AGROW>
        end
        fprintf('## (b) set=%-8s dw=%d D=%s w=%s Ngram=%s ndec=%d build %.1fs rows%s\n', ...
                s{1},dw,mat2str(info.degrees.D),mat2str(info.degrees.w),mat2str(info.N),info.ndec,tb,ints);
    end
end
%% (c) unit cost: faces o = 1 per direction, dw 0, f = 0.5
R = pf_pinned(2,{'DD','DD'},{'faces'},0,0.5,'probe');
fprintf('## (c) unit cost: build %.0f s, solve %.0f s, st %+d, m %d, nx %d, Ngram %s\n',R.tb,R.ts,R.st,R.m,R.nx,mat2str(R.Ngram));
disp('## PROBE DONE');
