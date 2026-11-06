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
%                           2-D; bottom: percentage of the total time spent
%                           parsing.
% INPUT
% - datadir: char, folder with ops_vs_q.csv and
%            stability_parse_vs_degree.csv;
% - outdir:  char, folder for the figures.
% OUTPUT
% none; the figure files.

% Initial coding MMP, 11/06/2026: scaling plots (README.md).

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

% Fig. 4: parsing time and percentage against the Lyapunov degree
ser = unique(Td.series,'stable');
cols = struct('legacy',[0 0.45 0.74],'sop_stock',[0 0 0],'sop',[0.20 0.63 0.17]);
mk = struct('legacy','o','sop_stock','s','sop','v');
lab = struct('legacy','opvar/dopvar (legacy)','sop_stock','sopvar/sdopvar, stock SDP','sop','sopvar/sdopvar');
f = figure('Visible','off','Position',[100 100 900 600]);
tl = tiledlayout(f,2,2,'TileSpacing','compact');
for nv = 1:2
    ax = nexttile(tl,nv);   hold(ax,'on');
    for s = ser'
        sel = Td.nv==nv & Td.series==s & isfinite(Td.t_parse);
        if ~any(sel),   continue,   end
        sc = char(s);
        semilogy(ax,Td.d(sel),Td.t_parse(sel),[':' mk.(sc)],'Color',cols.(sc),'MarkerFaceColor',cols.(sc),'DisplayName',lab.(sc));
    end
    set(ax,'YScale','log');     grid(ax,'on');
    dsv = unique(Td.d(Td.nv==nv));  xticks(ax,dsv);     xlim(ax,[min(dsv)-0.5, max(dsv)+0.5]);
    title(ax,sprintf('Stability LPI, %d-D',nv));
    xlabel(ax,'Degree of the Lyapunov Variable (d)');
    if nv==1,   ylabel(ax,'Parsing Time (s)');  legend(ax,'Location','northwest');  end
    ax = nexttile(tl,nv+2);
    ds = unique(Td.d(Td.nv==nv));   P = nan(numel(ds),numel(ser));
    for j = 1:numel(ser)
        for i = 1:numel(ds)
            sel = Td.nv==nv & Td.series==ser(j) & Td.d==ds(i);
            if any(sel),    P(i,j) = Td.pct_parse(find(sel,1));  end
        end
    end
    b = bar(ax,ds,P,'grouped');
    for j = 1:numel(ser),   b(j).FaceColor = cols.(char(ser(j)));   end
    ylim(ax,[0 100]);   grid(ax,'on');  xlabel(ax,'Degree of the Lyapunov Variable (d)');
    if nv==1,   ylabel(ax,'Time Spent Parsing (%)');  end
end
save_fig(f,outdir,'fig4_stability_parse');
end

function save_fig(f,outdir,name)
% PNG preview and vector PDF
exportgraphics(f,fullfile(outdir,[name '.png']),'Resolution',150);
exportgraphics(f,fullfile(outdir,[name '.pdf']),'ContentType','vector');
close(f);
end
