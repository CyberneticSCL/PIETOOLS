function [q,inf_] = sed_ref_pos(Atf,bf,Ns,Kf,Pst)
% Full-rank SeDuMi reference solve for the bare-positivity program, returning
% a Gram vector in NORMALISED-b units (opcheck_pos2d convention).  Same shape
% as pm_rank's sed_ref; failure returns [] and the CALLER decides what that
% means (never treated as an infeasibility proof by itself).
q = [];  inf_ = struct('numerr',-1); %#ok<NASGU> % kept for the catch path
try
    if Kf ~= 0, error('Kf=%d unexpected',Kf); end
    [RRm,Ksed] = cone_map(Ns,'psd');  pars.fid = 0;
    [xs,~,inf_] = sedumi(RRm'*Atf,bf,RRm'*sparse(size(Atf,1),1),Ksed,pars);
    xs = full(xs)/Pst.nb0;
    q = zeros(Pst.Ntot,1);  off = Kf;
    for i = 1:numel(Ns)
        Q = reshape(xs(off+1:off+Ns(i)^2),Ns(i),Ns(i));  off = off+Ns(i)^2;
        q(Pst.rows{i}) = reshape((Q+Q')/2,[],1);
    end
catch
    q = [];
end
end
