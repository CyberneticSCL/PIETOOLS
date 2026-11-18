function H = pf_heat(N,bc,d,ep,cfgs,RD,wQs,tag,bopts_in)
% H = PF_HEAT(N,BC,D,EP,CFGS,RD,WQS,TAG,BOPTS) the heat benchmark (Cor. 35,
% HEATND_LPI with RQ 'custom') bisected on kappa (HEATND_BISECT, certified
% verdicts) for every term-set configuration in CFGS (PF_CFG) and every
% weight degree wQ in WQS: R at (D, w) = RD (1 x 2, e.g. [1 1]), Q at
% (D, w) = (RD(1) + 1, wQ) per direction, the configuration's rows on both
% R and Q (psatzon 'RQ'), offsets applied per direction inside the scratch
% poscopvar_direct. Results are appended to results/heat_<TAG>.mat after
% every bisection and printed as '@@ HEAT' lines. SCRATCH, MMP, 11/15/2026.
S = pf_setup();
pie = heatNd_pie(N,0,bc);   lam1 = pie.exact.lambda1;
bopts = struct('rtol',1e-3,'maxsolve',16,'tmax',900,'solve',struct('tight',true),'verbose',false);
if nargin>8 && ~isempty(bopts_in)
    fn = fieldnames(bopts_in);  for i = 1:numel(fn),    bopts.(fn{i}) = bopts_in.(fn{i});   end
end
out = sprintf('%s/results/heat_%s.mat',S,tag);
H = struct('N',{},'bc',{},'d',{},'cfg',{},'o',{},'pd',{},'RD',{},'QD',{},'lo',{},'hi',{},'gap',{},'reach',{}, ...
           'resolved',{},'stop',{},'nsolve',{},'t',{},'m',{},'nx',{},'Ks',{},'nnz',{},'err',{});
for ic = 1:numel(cfgs)
    c = cfgs(ic);   o = NaN;    if ~isempty(c.offs),  o = c.offs(1);    end
    for wQ = wQs
        QD = [RD(1)+1, wQ];
        h = struct('N',N,'bc',{bc},'d',d,'cfg',c.name,'o',o,'pd',c.pd,'RD',RD,'QD',QD,'lo',NaN,'hi',NaN,'gap',[NaN NaN],'reach',NaN, ...
                   'resolved',false,'stop','','nsolve',NaN,'t',NaN,'m',NaN,'nx',NaN,'Ks',[],'nnz',NaN,'err','');
        try
            opts = struct('preset','bench','RQ','custom', ...
                          'Rdeg',struct('int',RD(2)*ones(1,N),'mult',RD(1)*ones(1,N)), ...
                          'Qdeg',struct('int',QD(2)*ones(1,N),'mult',QD(1)*ones(1,N)));
            if isempty(c.rows),     opts.psatz = 'none';
            else
                opts.psatz = 'rows';    opts.psatz_rows = c.rows;
                opts.psatz_rows_offset = c.offs;    opts.psatz_rows_perdir = c.pd;
            end
            B = heatNd_bisect(pie,d,ep,opts,bopts);
            h.lo = B.lo;    h.hi = B.hi;    h.gap = B.gap;  h.reach = B.lo/lam1;
            h.resolved = B.resolved;    h.stop = B.stop;    h.nsolve = B.nsolve;    h.t = B.t_total;
            if isfield(B,'sdp') && ~isempty(B.sdp)
                h.m = B.sdp.m;  h.nx = B.sdp.nx;    h.Ks = B.sdp.Ks;
                if isfield(B.sdp,'nnz'),    h.nnz = B.sdp.nnz;  end
            end
        catch ME
            h.err = strrep(ME.message,newline,' ');
        end
        H(end+1) = h;                                                       %#ok<AGROW>
        fprintf('@@ HEAT N=%d d=%d cfg=%s o=%g pd=%d RD=%s QD=%s lo=%.7f hi=%.7f reach=%.5f gap=[%.2e %.2e] res=%d nsolve=%d t=%.0f m=%d nx=%d Ks=%s nnz=%d stop=%s %s\n', ...
                N,d,c.name,o,c.pd,mat2str(RD),mat2str(QD),h.lo,h.hi,h.reach,h.gap(1),h.gap(2),h.resolved,h.nsolve,h.t,h.m,h.nx,mat2str(h.Ks),h.nnz,h.stop,h.err);
        save(out,'H');
    end
end
fprintf('## HEAT DONE %s\n',tag);
end
