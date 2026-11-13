function R = mx_exec3(cases,tier)
% R = MX_EXEC3(CASES,TIER) is MX_EXEC2 with the systems of MX_SYS3.
% INPUT
% - cases:  n x 5 cell as in MX_EXEC2;
% - tier:   (optional) char, the lpisettings_sop tier (default 'light').
% OUTPUT
% - R:      n x 1 struct array as in MX_EXEC.

% Initial coding MMP, 11/12/2026: runner of N12, N13b, N15, N17 (README.md).
if nargin<2, tier = 'light'; end
R = struct('sys',cases(:,1),'par','','exec',cases(:,3),'form',cases(:,4),'expect',NaN, ...
           'status','','certified',false,'gam',NaN,'gap',NaN,'verdict','','t',NaN, ...
           'm',NaN,'ndv',NaN,'Ks','','dw',NaN,'msg','');
for i=1:size(cases,1)
    R(i).par = jsonencode(cases{i,2});
    R(i).expect = double(cases{i,5});
    try
        S = mx_sys3(cases{i,1},cases{i,2});
        PIE = initialize(convert(S.PDE,{'pie','silent'}));
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
            case 'h2'
                [~,~,gam,~,~,info] = PIETOOLS_H2_norm_sop(PIE,st,cases{i,4});
        end
        R(i).t = toc(t0);
        R(i).status = info.status;  R(i).certified = info.certified;  R(i).gam = gam;
        if isfield(info,'shape') && ~isempty(info.shape)
            R(i).m = info.shape.m;  R(i).ndv = info.shape.ndv;  R(i).Ks = mat2str(info.shape.Ks(1:min(end,6)));
        end
        if isfield(info,'dw'), R(i).dw = info.dw; end
        if strcmp(cases{i,3},'stab')
            if R(i).certified && ~cases{i,5}
                R(i).verdict = 'UNSOUND';           % certified an unstable PIE
            elseif R(i).certified || ~cases{i,5}
                R(i).verdict = 'pass';
            else
                R(i).verdict = 'conservative';      % stable, not certified
            end
        else
            R(i).gap = gam/cases{i,5}-1;
            if R(i).certified && R(i).gap < -1e-6
                R(i).verdict = 'UNSOUND';           % bound below the exact value
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
    fprintf('%-26s %-28s %-4s %-3s | %-11s cert %d | gam %10.6g exact %10.6g gap %+9.2e | %-13s | %6.1fs m %6g ndv %8g Ks %s dw %g %s\n', ...
        r.sys,r.par,r.exec,r.form,r.status,r.certified,r.gam,r.expect,r.gap,r.verdict,r.t,r.m,r.ndv,r.Ks,r.dw,r.msg);
end
end
