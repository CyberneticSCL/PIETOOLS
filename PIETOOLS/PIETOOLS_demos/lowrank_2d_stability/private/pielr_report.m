function pielr_report(cert)                                                 % CC, 09/23/2026
% PIELR_REPORT  Print a pielr_solve result.
%
% Prints the maxima the residual ratio is built from, not only the ratio.
% opcheck_2d computed those on every call and returned the ratio alone, which
% is why no log this package produced can say whether two residuals differ
% because the numerator moved or because the denominator did -- and the
% denominator is free to move along any direction the equality system does
% not see.

fprintf('\n================ pielr_solve report ================\n');
fprintf('lpi    : %s   (%d-D)\n',cert.lpi,cert.dim);
fprintf('blocks : N = [%s]   free = %d   unknowns = %d\n', ...
    num2str(cert.Ns),cert.Kf,cert.unknowns_full);
if cert.ok
    fprintf('ranks  : r = [%s] of [%s]   (unknowns %d -> %d, %.0fx fewer)\n', ...
        num2str(cert.r),num2str(cert.Ns),cert.unknowns_full,cert.unknowns_face, ...
        cert.unknowns_full/max(cert.unknowns_face,1));
    fprintf('verdict: CERTIFIES   rel = %.4e   PSD = %d   mineig = %s\n', ...
        cert.rel,cert.psd,mat2str(cert.mineig,3));
    if isfield(cert,'score') && isfinite(cert.score)
        fprintf(['score  : rel/ref = %.4g   (ref %.4e, threshold %.4e = ' ...
                 'max(%.1e, %g x ref))\n'], ...
            cert.score,cert.ref_rel,cert.thresh,cert.gate_abs,cert.gate_k);
    else
        fprintf('score  : no reference supplied; threshold is the absolute %.1e\n', ...
            cert.thresh);
    end
    if isfield(cert,'rel_i') && numel(cert.rel_i) > 1
        fprintf(['         per-equality rel = %s   (the gate is the WORST of ' ...
                 'these; equality %d decided)\n'],mat2str(cert.rel_i,4),cert.worst_eq);
    end
    fprintf('         rel_d = %.4e  (against the operator that must vanish)\n',cert.rel_d);
    fprintf('         maxRes = %.4e   maxNrm = %.4e   maxDop = %.4e\n', ...
        cert.maxRes,cert.maxNrm,cert.maxDop);
    if ~isnan(cert.gam)
        fprintf('gamma  : %.6g   (UPPER bound: a face can only lose feasible points)\n',cert.gam);
    end
else
    fprintf('verdict: NOT CERTIFIED at rank <= the ladder ceiling -- see notes\n');
end
fprintf('times  : setup %.1fs | discovery %.1fs\n',cert.t_setup,cert.t_discover);
if isfield(cert,'why') && ~isempty(cert.why)
    n = numel(cert.why);
    ex = {cert.why.lm_exit};
    u = unique(ex);
    cnt = cellfun(@(e)sum(strcmp(ex,e)),u);
    fprintf('attempts: %d   LM exits: %s\n',n, ...
        strjoin(arrayfun(@(k)sprintf('%s x%d',u{k},cnt(k)),1:numel(u),'uni',0),', '));
    rM = [cert.why.rM];  np = [cert.why.np];
    if any(~isnan(rM))
        fprintf('face    : np in [%s]  rM in [%s]  (rM <= np always: rank of an m x np matrix)\n', ...
            rng_str(np),rng_str(rM));
    end
    re = [cert.why.rel_eq_full];
    if any(~isnan(re))
        fprintf('          full-row equality residual on the face: min %.3e  max %.3e\n', ...
            min(re(~isnan(re))),max(re(~isnan(re))));
    end
end
if isfield(cert,'nrm_why') && ~isempty(cert.nrm_why)
    fprintf('denom  : %s\n',wrapline(cert.nrm_why,9));
end
if isfield(cert,'notes')
    for i = 1:numel(cert.notes)
        fprintf('note   : %s\n',wrapline(cert.notes{i},9));
    end
end
fprintf('====================================================\n\n');
end

function s = rng_str(v)
v = v(~isnan(v));
if isempty(v), s = '-'; elseif numel(unique(v))==1, s = sprintf('%g',v(1));
else, s = sprintf('%g..%g',min(v),max(v)); end
end

function s = wrapline(s,ind)
% keep long notes readable in a terminal without truncating them
w = 72;  pad = repmat(' ',1,ind);
out = {};  while numel(s) > w
    k = find(s(1:w)==' ',1,'last');   if isempty(k), k = w; end
    out{end+1} = s(1:k-1); s = s(k+1:end);   %#ok<AGROW>
end
out{end+1} = s;
s = strjoin(out,[newline pad]);
end
