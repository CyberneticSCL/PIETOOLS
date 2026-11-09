function bl_bisdump(ids)
% bl_bisdump -- write the FEASIBILITY programs a bisection would visit.
%
% For each objective-form case, gamma is pinned at fractions of the value the
% objective form returned, straddling it. Two of those fractions are below
% gamma* and therefore INFEASIBLE by construction, which is the point: a
% bisection is only as good as its ability to recognise the infeasible side, and
% Mosek recognises it with a certificate (PRIMAL_INFEASIBLE_CER, verified on
% hinf_rd1: rel_b exactly 1.0 and ||x|| ~ 1e-18 below gamma*, OPTIMAL above).
% cuADMM has no infeasibility certificate at all, so whether it can drive a
% bisection is an open question these dumps exist to settle.
%
% Fractions are relative to each case's own gamma*, so every case is probed at
% the same relative distance from its own boundary rather than at absolute
% gammas that would mean different things per case.
%
% CC, 09/26/2026: H2 dropped from the testing regime (maintainer decision):
%   h2c_rd1, h2oco_rd1 and h2est_rd1 removed from the gamma* map; old map kept
%   in a comment above it. The gam^2 note on the coercive H2 executives stays,
%   for anyone restoring them.

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
DMP  = fullfile(cuadmm_outdir(),'baseline','dumps');
BIS  = fullfile(cuadmm_outdir(),'baseline','bisect');  if ~exist(BIS,'dir'), mkdir(BIS); end
TSV  = fullfile(cuadmm_outdir(),'baseline','bl_bisect.tsv');
if ~exist(TSV,'file')
    fid=fopen(TSV,'w');
    fprintf(fid,'id\tfrac\tgam\tm\tstatus\tprosta\tsolsta\tt_mosek\trel_b\tpsd_relmin\tnormx\tfeas\n');
    fclose(fid);
end

% gamma* as returned by the objective form (bl_mosek.tsv), per case
% These are values of the SDP OBJECTIVE VARIABLE -- the thing bl_fix pins -- which
% is NOT always the gamma the executive returns. PIETOOLS_H2_norm_o_coercive.m:172
% and _c_coercive.m:171 return gam = sqrt(objective), so for those the variable
% is gam^2. h2oco_rd1 was first pinned at its returned gam (9.03e-05), i.e. ~9000x
% above its optimum (8.16e-09): all four rungs sat deep in the feasible region and
% the "feasible even below gamma*" reading was an artefact of that. Every other
% executive listed here returns the objective directly (checked, 2026-09-25).
% CC, 09/26/2026: H2 entries removed (see header). Was:
%  {'hinf_rd1','hinfdu_rd1','hinf_rd1_hv','h2c_rd1','h2oco_rd1','est_rd1','h2est_rd1','hinf2_rd'}, ...
%  { 0.182626,  0.182551,    0.182482,     0.288786, 9.03171e-05^2, 0.000211481, 0.000171429, 1.27221});
G = containers.Map( ...
  {'hinf_rd1','hinfdu_rd1','hinf_rd1_hv','est_rd1','hinf2_rd'}, ...
  { 0.182626,  0.182551,    0.182482,     0.000211481, 1.27221});           % CC, 09/26/2026
FR = [0.80 0.99 1.01 1.20];     % two infeasible, two feasible, by construction

if nargin<1 || isempty(ids), ids = G.keys; end
for q = 1:numel(ids)
    id = ids{q};
    if ~isKey(G,id), fprintf('BS skip %s (no gamma*)\n',id); continue; end
    gs = G(id);
    for f = FR
        g = f*gs;
        lab = sprintf('%s_g%03.0f',id,f*100);
        try
            D = bl_fix(fullfile(DMP,[id '.mat']),g,fullfile(BIS,lab));
            prob = Sedumi2Mosek(D.At',full(D.b),D.c,D.K);
            t0=tic; [~,res]=mosekopt('minimize info echo(0)',prob); tm=toc(t0);
            x = MosekSol2SedumiSol(D.K,res); x = x(:);
            rb = norm(full(D.At'*x-D.b))/max(norm(full(D.b)),eps);
            off=D.K.f; pmin=inf; pmax=-inf;
            for k=1:numel(D.Ks)
                N=double(D.Ks(k)); Xk=reshape(x(off+(1:N^2)),N,N);
                ev=eig((Xk+Xk')/2); pmin=min(pmin,min(ev)); pmax=max(pmax,max(ev)); off=off+N^2;
            end
            pr = res.sol.itr.prosta; so = res.sol.itr.solsta;
            feas = double(isempty(strfind(pr,'INFEASIBLE')));
            st_s = 'ok';
        catch ME
            pr=''; so=''; tm=NaN; rb=NaN; pmin=NaN; pmax=1; x=NaN; feas=NaN;
            st_s = ['ERR ' regexprep(ME.message,'[\t\r\n]+',' ')];
        end
        fid=fopen(TSV,'a');
        fprintf(fid,'%s\t%.2f\t%.8g\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
            lab,f,g,num2str(D.m),st_s,pr,so,num2str(tm),num2str(rb), ...
            num2str(pmin/max(pmax,eps)),num2str(norm(x)),num2str(feas));
        fclose(fid);
        fprintf('BS %-22s f=%.2f gam=%-12.6g %-26s feas=%s t=%.3f rel_b=%.3e\n', ...
            lab,f,g,pr,num2str(feas),tm,rb);
    end
end
fprintf('BSDONE\n');
end
