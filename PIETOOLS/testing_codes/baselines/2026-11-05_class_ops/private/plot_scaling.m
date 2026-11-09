function plot_scaling(datadir,outdir)
% PLOT_SCALING(DATADIR,OUTDIR) draws scaling plots in the style of Figs. 2-4
% of Jagt et al., IEEE Control Systems Letters 6 (2022) 3493-3498, from the
% CSV files of EXPORT_SCALING_CSV, and saves them as PNG and PDF in OUTDIR:
%   fig2_composition_vs_q   time of the composition D*A of a decision
%                           operator D with a constant operator A against
%                           the number q of decision variables: (a) nnz
%                           fixed, (b) nnz about 4q; dopvar and sdopvar;
%   fig3_ops_vs_q           (a) addition D + E, disjoint and shared lists of
%                           decision variables, (b) adjoint D', against q;
%   fig4_stability_parse    top: parsing time of the stability LPI against
%                           the degree d of the Lyapunov variable, 1-D and
%                           2-D (filled markers, dotted: runs with a full
%                           solve; open markers, dashed: parse-only runs);
%                           bottom: percentage of the total time spent
%                           parsing, full solves only;
%   fig5_of_parse           the same for the output feedback LPI, when
%                           of_parse_vs_degree.csv exists.
% INPUT
% - datadir: char, folder with ops_vs_q.csv, stability_parse_vs_degree.csv
%            and optionally of_parse_vs_degree.csv;
% - outdir:  char, folder for the figures.
% OUTPUT
% none; the figure files.

% Initial coding MMP, 11/06/2026: scaling plots (README.md).
% MMP, 11/06/2026: the LPI figures by one local function LPI_FIGURE, which
%                replaces the inline Fig. 4 block; parse-only runs (stub
%                solver) on their own dashed line with open markers and no
%                percentage bar, since they exclude what a full solve counts
%                as parsing after the solve; y-limits padded so no marker sits
%                on the frame, 1-2-5 ticks under two decades; one legend
%                above the panels; Fig. 5 for the output feedback LPI.

Tq = readtable(fullfile(datadir,'ops_vs_q.csv'),'TextType','string');
Td = readtable(fullfile(datadir,'stability_parse_vs_degree.csv'),'TextType','string');
sty = struct('dopvar',{{'o',[0 0.45 0.74]}},'sdopvar',{{'v',[0.20 0.63 0.17]}});

% Fig. 2: composition against q
f = figure('Visible','off','Position',[100 100 900 340]);
tl = tiledlayout(f,1,2,'TileSpacing','compact');
ttl = {'nnz fixed','nnz \approx 4q'};    sw = ["nnz_fixed","nnz_4q"];
for p = 1:2
    ax = nexttile(tl);  hold(ax,'on');
    for c = ["dopvar","sdopvar"]
        sel = Tq.sweep==sw(p) & Tq.op=="mtimes_DA" & Tq.class==c;
        s = sty.(char(c));
        loglog(ax,Tq.q(sel),Tq.t(sel),[':' s{1}],'Color',s{2},'MarkerFaceColor',s{2},'DisplayName',char(c));
    end
    set(ax,'XScale','log','YScale','log');  grid(ax,'on');
    xlabel(ax,'Number of Decision Variables (q)');   title(ax,['Composition D\cdotA, ' ttl{p}]);
    if p==1,    ylabel(ax,'Operation Execution Time (s)');  legend(ax,'Location','northwest');  end
end
save_fig(f,outdir,'fig2_composition_vs_q');

% Fig. 3: addition (disjoint and shared lists) and adjoint against q
f = figure('Visible','off','Position',[100 100 900 340]);
tl = tiledlayout(f,1,2,'TileSpacing','compact');
ax = nexttile(tl);  hold(ax,'on');
for c = ["dopvar","sdopvar"]
    s = sty.(char(c));
    sel = Tq.sweep=="nnz_fixed" & Tq.op=="plus" & Tq.class==c;
    loglog(ax,Tq.q(sel),Tq.t(sel),[':' s{1}],'Color',s{2},'MarkerFaceColor',s{2},'DisplayName',char(c)+", disjoint lists");
    sel = Tq.sweep=="nnz_fixed_shared" & Tq.op=="plus" & Tq.class==c;
    loglog(ax,Tq.q(sel),Tq.t(sel),['-' s{1}],'Color',s{2},'DisplayName',char(c)+", shared list");
end
set(ax,'XScale','log','YScale','log');  grid(ax,'on');  title(ax,'Addition D + E');
xlabel(ax,'Number of Decision Variables (q)');   ylabel(ax,'Operation Execution Time (s)');
legend(ax,'Location','northwest');
ax = nexttile(tl);  hold(ax,'on');
for c = ["dopvar","sdopvar"]
    s = sty.(char(c));
    sel = Tq.sweep=="nnz_fixed" & Tq.op=="ctranspose" & Tq.class==c;
    loglog(ax,Tq.q(sel),Tq.t(sel),[':' s{1}],'Color',s{2},'MarkerFaceColor',s{2},'DisplayName',char(c));
end
set(ax,'XScale','log','YScale','log');  grid(ax,'on');  title(ax,'Adjoint D''');
ta = Tq.t(Tq.sweep=="nnz_fixed" & Tq.op=="ctranspose");     % whole decades, as in the other panels
ylim(ax,[10^floor(log10(min(ta))), 10^(ceil(log10(max(ta)))+1)]);
xlabel(ax,'Number of Decision Variables (q)');
save_fig(f,outdir,'fig3_ops_vs_q');

% MMP, 11/06/2026: start of the LPI figures. Replaces the inline Fig. 4 block
% of 7688dcf2 (one dotted line per series through every row, two panels for
% nv = 1:2).
% Fig. 4: the stability LPI; Fig. 5: the output feedback LPI (if measured)
lab = struct('legacy','opvar/dopvar (legacy)','sop_stock','sopvar/sdopvar, stock SDP','sop','sopvar/sdopvar');
lpi_figure(Td,'Stability LPI, %d-D',lab,outdir,'fig4_stability_parse');
fo = fullfile(datadir,'of_parse_vs_degree.csv');
if isfile(fo)
    To = readtable(fo,'TextType','string');
    lab.sop_stock = 'sopvar/sdopvar, stock-type slack';
    lpi_figure(To,'Output feedback LPI, %d-D',lab,outdir,'fig5_of_parse');
end
end

function lpi_figure(Td,ttl,lab,outdir,name)
% Top: parsing time against the Lyapunov degree, one panel per dimension.
% Runs with a full solve: filled markers, dotted. Parse-only runs (the stub
% solver): open markers, dashed; a separate line, since parse-only time
% excludes the post-solve work a full solve counts as parsing. Bottom:
% percentage of the total time spent parsing, full solves only.
po = Td.parse_only;
if ~islogical(po)
    if isnumeric(po),   po = po~=0;     else,   po = strcmpi(string(po),"true") | string(po)=="1";  end
end
ser = unique(Td.series,'stable');   nvs = unique(Td.nv)';
cols = struct('legacy',[0 0.45 0.74],'sop_stock',[0 0 0],'sop',[0.20 0.63 0.17]);
mk = struct('legacy','o','sop_stock','s','sop','v');
f = figure('Visible','off','Position',[100 100 450*numel(nvs) 600]);
tl = tiledlayout(f,2,numel(nvs),'TileSpacing','compact');
for k = 1:numel(nvs)
    nv = nvs(k);
    ax = nexttile(tl,k);    hold(ax,'on');
    for s = ser'
        sc = char(s);   shown = false;      % one legend entry per series
        for p = [false true]    % full solves, then parse-only runs
            sel = find(Td.nv==nv & Td.series==s & po==p & isfinite(Td.t_parse));
            if isempty(sel),    continue,   end
            [ds,ix] = sort(Td.d(sel));
            if p,   lsty = '--';    fc = 'w';   else,   lsty = ':';     fc = cols.(sc);     end
            h = semilogy(ax,ds,Td.t_parse(sel(ix)),[lsty mk.(sc)],'Color',cols.(sc),'MarkerFaceColor',fc,'DisplayName',lab.(sc));
            if shown,   h.HandleVisibility = 'off';     end
            shown = true;
        end
    end
    set(ax,'YScale','log');     grid(ax,'on');
    tp = Td.t_parse(Td.nv==nv);     tp = tp(isfinite(tp) & tp>0);
    yl = [min(tp)/1.5, 1.5*max(tp)];    ylim(ax,yl);    % no marker on the frame
    if yl(2)/yl(1) < 100    % under two decades MATLAB labels one decade at most
        yt = 10.^(floor(log10(yl(1))):ceil(log10(yl(2))))'*[1 2 5];
        yt = sort(yt(:));   yt = yt(yt>=yl(1) & yt<=yl(2));
        yticks(ax,yt);  yticklabels(ax,compose('%g',yt));
    end
    dsv = unique(Td.d(Td.nv==nv));  xticks(ax,dsv);     xlim(ax,[min(dsv)-0.5, max(dsv)+0.5]);
    title(ax,sprintf(ttl,nv));
    xlabel(ax,'Degree of the Lyapunov Variable (d)');
    if k==1,    ylabel(ax,'Parsing Time (s)');  end
    ax = nexttile(tl,k+numel(nvs));
    ds = unique(Td.d(Td.nv==nv));   P = nan(numel(ds),numel(ser));
    for j = 1:numel(ser)
        for i = 1:numel(ds)
            sel = Td.nv==nv & Td.series==ser(j) & Td.d==ds(i) & ~po;
            if any(sel),    P(i,j) = Td.pct_parse(find(sel,1));  end
        end
    end
    b = bar(ax,ds,P,'grouped');
    for j = 1:numel(ser),   b(j).FaceColor = cols.(char(ser(j)));   end
    ylim(ax,[0 100]);   grid(ax,'on');  xlabel(ax,'Degree of the Lyapunov Variable (d)');
    if k==1,    ylabel(ax,'Time Spent Parsing (%)');  end
end
% One legend for the figure, above the panels: a legend inside a 2-D panel
% covers data, whose lines cross both the upper left and the lower right
ax = nexttile(tl,1);    hs = gobjects(0);
for s = ser'
    sc = char(s);
    hs(end+1) = plot(ax,NaN,NaN,[':' mk.(sc)],'Color',cols.(sc),'MarkerFaceColor',cols.(sc),'DisplayName',lab.(sc)); %#ok<AGROW>
end
if any(po)
    hs(end+1) = plot(ax,NaN,NaN,'--o','Color',[0.4 0.4 0.4],'MarkerFaceColor','w','DisplayName','open, dashed: parse only, no solve');
end
lg = legend(ax,hs,'NumColumns',numel(hs));     lg.Layout.Tile = 'north';
save_fig(f,outdir,name);
end
% MMP, 11/06/2026: end of the LPI figures.

function save_fig(f,outdir,name)
% PNG preview and vector PDF
exportgraphics(f,fullfile(outdir,[name '.png']),'Resolution',150);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector');
close(f);
end
