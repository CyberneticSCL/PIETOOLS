function R = mx_exec(cases,tier)
% R = MX_EXEC(CASES,TIER) runs PIETOOLS_stability_sop or
% PIETOOLS_Hinf_gain_sop on systems of MX_SYS, on the tree set by
% MX_FROZEN_PATH, and prints one line per case (parsed by MX_ROWS2CSV.AWK).
% INPUT
% - cases:  n x 5 cell, one row per case: system name (MX_SYS), 'struct' of
%           parameters, 'stab' or 'hinf', the form, and the expected value
%           (true/false, stable or not, for 'stab'; the exact gain for
%           'hinf');
% - tier:   (optional) char, the lpisettings_sop tier (default 'light').
% OUTPUT
% - R:      n x 1 struct array: status, certified, gam, gap, verdict
%           ('pass', 'conservative', 'not certified', 'UNSOUND', 'ERROR'),
%           time, rows m, decision variables, PSD block sizes, dw.

% Initial coding MMP, 11/12/2026: runner of N3-N10 (README.md). Converts
% inside evalc; MX_EXEC2 converts with 'silent' instead (R2025b crashed in
% evalc in N2).
if nargin<2, tier = 'light'; end
R = struct('sys',cases(:,1),'par','','exec',cases(:,3),'form',cases(:,4),'expect',NaN, ...
           'status','','certified',false,'gam',NaN,'gap',NaN,'verdict','','t',NaN, ...
           'm',NaN,'ndv',NaN,'Ks','','dw',NaN,'msg','');
for i=1:size(cases,1)
    R(i).par = jsonencode(cases{i,2});
    R(i).expect = double(cases{i,5});
    try
        S = mx_sys(cases{i,1},cases{i,2});
        evalc('PIE = initialize(convert(S.PDE));');
        st = lpisettings_sop(tier);
        if any(strcmp(getappdata(0,'mx_tree'),{'t2','t3'}))
            st.lpi.ineq.psatz = 'corners';  % opt-in corner rows (live lpi_ineq_sop)
        end
        t0 = tic;
        switch cases{i,3}
            case 'stab'
                [~,~,info] = PIETOOLS_stability_sop(PIE,st,cases{i,4});
                gam = NaN;
            case 'hinf'
                [~,~,gam,info] = PIETOOLS_Hinf_gain_sop(PIE,st,cases{i,4});
        end
        R(i).t = toc(t0);
        R(i).status = info.status;  R(i).certified = info.certified;  R(i).gam = gam;
        if isfield(info,'shape') && ~isempty(info.shape)
            R(i).m = info.shape.m;  R(i).ndv = info.shape.ndv;  R(i).Ks = mat2str(info.shape.Ks(1:min(end,6)));
        end
        if isfield(info,'dw'), R(i).dw = info.dw; end
        switch cases{i,3}
            case 'stab'
                if R(i).certified && ~cases{i,5}
                    R(i).verdict = 'UNSOUND';           % certified an unstable PIE
                elseif R(i).certified
                    R(i).verdict = 'pass';
                elseif ~cases{i,5}
                    R(i).verdict = 'pass';              % no certificate for an unstable PIE
                else
                    R(i).verdict = 'conservative';      % stable, not certified
                end
            case 'hinf'
                R(i).gap = gam/cases{i,5}-1;
                if R(i).certified && R(i).gap < -1e-6
                    R(i).verdict = 'UNSOUND';           % bound below the exact gain
                elseif R(i).certified
                    R(i).verdict = 'pass';
                else
                    R(i).verdict = 'not certified';
                end
        end
    catch ME
        R(i).msg = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
        R(i).verdict = 'ERROR';
    end
    r = R(i);
    fprintf('%-26s %-28s %-4s %-2s | %-11s cert %d | gam %10.6g exact %10.6g gap %+9.2e | %-12s | %6.1fs m %6g ndv %8g Ks %s dw %g %s\n', ...
        r.sys,r.par,r.exec,r.form,r.status,r.certified,r.gam,r.expect,r.gap,r.verdict,r.t,r.m,r.ndv,r.Ks,r.dw,r.msg);
end
end
