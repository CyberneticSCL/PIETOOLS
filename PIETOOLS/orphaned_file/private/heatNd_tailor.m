function T = heatNd_tailor(N,dlist,bopts,bc,grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = HEATND_TAILOR(N,DLIST,BOPTS) does the proof program's degree
% tailoring of R and Q reach the stock bounds of the heat benchmark with a
% smaller SDP? For every P degree d in DLIST:
%   - the stock presets 'bench' (Tspan dp 0, two faces at the same degree)
%     and 'heavy' (dp 1) are bisected on kappa (HEATND_BISECT);
%   - the targets of R and Q, E1 = P*T - ep^2 T*T and X1 + X2, are read
%     for their lower-kernel support {s^i s'^j} and multiplier degree:
%       Dmin = max min(i,j),  Dmax = max max(i,j),  Mdeg = max multiplier i;
%   - the tailored grid: lift D = Dmin + dD, weight w = max(ceil((Dmax-D-1)
%     /2), ceil(Mdeg/2), 0) + dw, (dD,dw) in {0,1}^2 (the necessary
%     conditions of sopvar_lift_notes Sec. 7: lower kernel in {min(i,j) <= D,
%     max(i,j) <= D+2w+1}, multiplier degree <= 2w), each R and Q from its
%     own target, no joint cap, with the Psatz choices
%       none | product at w-1 (Markov-Lukacs pair) | product at w |
%       faces at w (stock) | faces at w-1,
%     through HEATND_LPI opts.RQ 'custom' and opts.psatz_offset.
% Every configuration is bisected with the same BOPTS; the SDP shape (m
% rows, Gram sizes Ks, nx, decision variables) comes from HEATND_BISECT.
% At the end, for each d, the smallest SDP (by nx, then m) reaching the
% bench bound and the heavy bound within 1e-3 in kappa is printed.
%
% INPUT: N (default 1), DLIST (default 0:2), BOPTS for HEATND_BISECT
%   (defaults here: rtol 1e-4, atol 1e-6, maxsolve 24, retry tight/rows,
%   verbose false), BC (cell of per-direction boundary conditions for
%   HEATND_PIE, default its own: DD in s1, DN elsewhere), GRID (struct with
%   fields dD, dw (vectors, default [0 1]), psz (cell rows {name, offset},
%   default the five above) and stock (cellstr of presets, default bench
%   and heavy) to select part of the grid, e.g. for 2-D under a time
%   budget). ep = 0.1 (paper), r = 0. In N-D every rule is per direction and
%   dD, dw add to every direction.
% OUTPUT: struct array T, one row per configuration: d, label, RD, Rw, QD,
%   Qw, psatz, offset, lo, hi, gap, stop, resolved, m, Ksmax, sumKs2, nx,
%   ndec, nsolve, t.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(N),      N = 1;      end
if nargin<2 || isempty(dlist),  dlist = 0:2;    end
if nargin<3 || isempty(bopts),  bopts = struct();   end
bdef = struct('rtol',1e-4,'atol',1e-6,'maxsolve',24,'verbose',false,'retry',{{'tight','rows'}});
fn = fieldnames(bdef);
for i = 1:numel(fn),    if ~isfield(bopts,fn{i}),   bopts.(fn{i}) = bdef.(fn{i});   end,    end
ep = 0.1;
if nargin<4 || isempty(bc),     pie = heatNd_pie(N,0);  else,   pie = heatNd_pie(N,0,bc);   end
lam1 = pie.exact.lambda1;
fprintf('heatNd_tailor: N = %d, %s, lambda1 = %.6f, ep = %g\n',N,strjoin(pie.bc,'/'),lam1,ep);
psz = {'none',0; 'product',1; 'product',0; 'linear',0; 'linear',1};
dDl = 0:1;  dwl = 0:1;  stock = {'bench','heavy'};
if nargin>=5 && ~isempty(grid)
    if isfield(grid,'dD'),      dDl = grid.dD;      end
    if isfield(grid,'dw'),      dwl = grid.dw;      end
    if isfield(grid,'psz'),     psz = grid.psz;     end
    if isfield(grid,'stock'),   stock = grid.stock; end
end
T = struct('d',{},'label',{},'RD',{},'Rw',{},'QD',{},'Qw',{},'psatz',{},'offset',{}, ...
           'lo',{},'hi',{},'gap',{},'stop',{},'resolved',{},'m',{},'Ksmax',{},'sumKs2',{}, ...
           'nx',{},'ndec',{},'nsolve',{},'t',{});
for d = dlist
    % % % Stock presets.
    for pr = reshape(stock,1,[])
        try
            T(end+1) = run_one(pie,d,ep,struct('preset',pr{1}),bopts,pr{1},NaN(1,2*N),NaN(1,2*N)); %#ok<AGROW>
            show(T(end));
        catch ME
            fprintf('  %-34s FAILED: %s\n',pr{1},ME.message);
        end
    end
    % % % The targets of R and Q, from one base build.
    [~,meta] = heatNd_lpi(pie,d,[],ep,struct('preset','bench'));
    sE = support_degrees(meta.base.E1.C{1,1});
    XQ = meta.base.X1 + meta.base.X2;
    sX = support_degrees(XQ.C{1,1});
    fprintf('  d = %d targets: E1 (R): Dmin %s Dmax %s Mdeg %s;  X1+X2 (Q): Dmin %s Dmax %s Mdeg %s\n', ...
            d,mat2str(sE.Dmin),mat2str(sE.Dmax),mat2str(sE.Mdeg),mat2str(sX.Dmin),mat2str(sX.Dmax),mat2str(sX.Mdeg));
    % % % The tailored grid.
    for dD = reshape(dDl,1,[])
        for dw = reshape(dwl,1,[])
            [DR,wR] = rule(sE,dD,dw);   [DQ,wQ] = rule(sX,dD,dw);
            for ip = 1:size(psz,1)
                o = struct('preset','bench','RQ','custom', ...
                           'Rdeg',struct('int',wR,'mult',DR),'Qdeg',struct('int',wQ,'mult',DQ), ...
                           'psatz',psz{ip,1},'psatz_offset',psz{ip,2});
                lbl = sprintf('tailored dD%d dw%d %s-%d',dD,dw,psz{ip,1},psz{ip,2});
                try
                    T(end+1) = run_one(pie,d,ep,o,bopts,lbl,[DR wR],[DQ wQ]);     %#ok<AGROW>
                catch ME
                    fprintf('  %-34s FAILED: %s\n',lbl,ME.message);
                    continue
                end
                show(T(end));
            end
        end
    end
end
% % % Summary: smallest SDP reaching each stock bound.
fprintf('\nSummary (within 1e-3 in kappa of the stock bound; size by nx, then m):\n');
for d = dlist
    rows = T([T.d]==d);
    for pr = {'bench','heavy'}
        ref = rows(strcmp({rows.label},pr{1}));
        if isempty(ref) || isnan(ref.lo),   continue,   end
        ok = rows([rows.lo]>=ref.lo-1e-3);
        [~,o] = sortrows([[ok.nx]',[ok.m]']);
        best = ok(o(1));
        fprintf('  d = %d, %-5s bound %.6f: stock nx %d m %d ndec %d;  smallest reaching it: %-30s nx %d m %d ndec %d (kappa %.6f)\n', ...
                d,pr{1},ref.lo,ref.nx,ref.m,ref.ndec,best.label,best.nx,best.m,best.ndec,best.lo);
    end
    [~,o] = max([rows.lo]);
    fprintf('  d = %d, highest kappa: %-30s %.6f (nx %d, m %d), gap to lambda1 %.2e\n',d,rows(o).label,rows(o).lo,rows(o).nx,rows(o).m,lam1-rows(o).lo);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = run_one(pie,d,ep,o,bopts,lbl,RDw,QDw)
% RDw = [D (1 x N), w (1 x N)] of R, QDw of Q; NaN for a preset.
B = heatNd_bisect(pie,d,ep,o,bopts);
N = numel(RDw)/2;
r = struct('d',d,'label',lbl,'RD',RDw(1:N),'Rw',RDw(N+1:end),'QD',QDw(1:N),'Qw',QDw(N+1:end), ...
           'psatz',o_field(o,'psatz','preset'),'offset',o_field(o,'psatz_offset',0), ...
           'lo',B.lo,'hi',B.hi,'gap',B.gap(2),'stop',B.stop,'resolved',B.resolved, ...
           'm',B.sdp.m,'Ksmax',max([B.sdp.Ks,0]),'sumKs2',sum(B.sdp.Ks.^2),'nx',B.sdp.nx, ...
           'ndec',B.meta.ndec,'nsolve',B.nsolve,'t',B.t_total);
end


function v = o_field(o,f,dflt)
if isfield(o,f),    v = o.(f);  else,   v = dflt;   end
if ischar(dflt) && ~ischar(v),  v = dflt;   end
end


function show(r)
if ischar(r.psatz),     ps = r.psatz;   else,   ps = 'preset';  end
fprintf('  %-34s R(D,w)=(%s,%s) Q(D,w)=(%s,%s) %-7s off %d | kappa [%.6f, %.6f) gap %.2e %-9s | m %5d Ks<=%3d sumKs2 %6d nx %6d ndec %6d | %d solves %.1f s\n', ...
        r.label,nstr(r.RD),nstr(r.Rw),nstr(r.QD),nstr(r.Qw),ps,r.offset,r.lo,r.hi,r.gap,r.stop(1:min(9,end)), ...
        r.m,r.Ksmax,r.sumKs2,r.nx,r.ndec,r.nsolve,r.t);
end


function s = nstr(x)
if any(isnan(x)),   s = '-';
elseif isscalar(x), s = sprintf('%d',x);
else,               s = mat2str(x);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = support_degrees(B)
% Per direction d: the support {(i_d,j_d)} of the cells that are lower in
% d, Dmin(d) = max min, Dmax(d) = max max, and the degree Mdeg(d) of the
% cells that are multipliers in d, over the constant part and every
% decision-variable row of the sdopvar block B (one space, all variables
% shared). Lower and upper cells are adjoint images, so the lower ones
% suffice.
nv = numel(B.ZL);
EL = multiindex_grid(B.ZL,'first_slowest');     NL = size(EL,1);
ER = multiindex_grid(B.ZR,'first_slowest');     NR = size(ER,1);
m = B.dims(1);  n = B.dims(2);
s = struct('Dmin',zeros(1,nv),'Dmax',zeros(1,nv),'Mdeg',-ones(1,nv));
for q = 1:numel(B.params.A)
    A = B.params.A{q};  Bq = B.params.B{q};
    cols = find(A);
    if ~isempty(Bq),    [~,cB] = find(Bq);   cols = unique([cols(:); cB(:)]);  end
    if isempty(cols),   continue,   end
    [ridx,cidx] = ind2sub([m*NL,n*NR],cols);
    l = mod(ridx-1,NL)+1;   e = mod(cidx-1,NR)+1;
    gam = gamma_of_cell(q,nv);
    for d = 1:nv
        i = EL(l,d);    j = ER(e,d);
        if gam(d)==1
            s.Mdeg(d) = max([s.Mdeg(d); i(:)]);
        elseif gam(d)==2
            s.Dmin(d) = max([s.Dmin(d); min(i(:),j(:))]);
            s.Dmax(d) = max([s.Dmax(d); max(i(:),j(:))]);
        end
    end
end
end


function [D,w] = rule(s,dD,dw)
% Per direction: lift D = Dmin + dD; weight w from the reach D + 2w + 1 >=
% Dmax and the multiplier degree <= 2w, plus dw.
D = s.Dmin + dD;
w = max([ceil((s.Dmax-D-1)/2); ceil(s.Mdeg/2); zeros(size(D))],[],1) + dw;
end
