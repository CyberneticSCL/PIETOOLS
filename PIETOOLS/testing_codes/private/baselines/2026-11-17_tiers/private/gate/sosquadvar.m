function varargout = sosquadvar(sos,Z1c,Z2c,varargin)
% SOSQUADVAR size gate of the library pass at every Gram declaration
% (baselines/2026-11-17_tiers; MMP, 11/23/2026): the repository SOSQUADVAR,
% then the cap LB_NDV_CAP of gate/lpisolve.m on the rows of decvartable.
% The gate of lpisolve.m acts once the whole program is built: PDE 55 on
% 'light' and PDE 54 on 'extreme' (28 M and 24.9 M decision variables)
% took about 40 GB before it (stress run of 11/22/2026). Here a program
% stops at the first Gram declaration that takes it above the cap: that
% Gram is declared, no positive operator is built on it. decvartable only
% grows, so a program stopped here is one the solve gate refuses. Free
% decision variables (lpidecvar) and LB_M_CAP are gated at the solve only.
% Without the appdata entry nothing is capped. lb_run puts this folder
% first on the path after lb_setup.
persistent repo
if isempty(repo)
    % A handle made while the current folder holds the repository file
    % stays bound to that file after the folder and the path change
    w = which('sosquadvar','-all');
    w = w(endsWith(w,'.m') & ~strcmpi(w,[mfilename('fullpath'),'.m']));
    o = cd(fileparts(w{1}));    repo = @sosquadvar;     cd(o);
end
[varargout{1:max(nargout,1)}] = repo(sos,Z1c,Z2c,varargin{:});
ndv = size(varargout{1}.decvartable,1);
if isappdata(0,'LB_NDV_CAP') && ndv > getappdata(0,'LB_NDV_CAP')
    error('lb:size','LB_SIZE ndv=%d after a Gram declaration, above the gate (ndv %g): not built', ...
          ndv,getappdata(0,'LB_NDV_CAP'));
end
end
