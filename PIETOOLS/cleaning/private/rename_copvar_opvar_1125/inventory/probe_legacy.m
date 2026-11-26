% Probe for the Phase 0 review: (1) what the legacy 'state'/'terms' format does at HEAD;
% (2) what the legacy lines "opvar T; T.<field> = ..." would do after a plain rename,
% simulated with today's name 'copvar' (the rename changes names only).
root = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
addpath(genpath(root));
fprintf('which copvar : %s\n', which('copvar'));
fprintf('which opvar  : %s\n', which('opvar'));
fprintf('which state  : %s\n', which('state'));

% (1) the legacy format at HEAD
try
    x = state('pde'); fprintf('L1 state(''pde''): ok, class %s\n', class(x));
catch e, fprintf('L1 state(''pde''): error [%s] %s\n', e.identifier, e.message); end
try
    y = x + x; fprintf('L2 x+x: ok, class %s\n', class(y));
catch e, fprintf('L2 x+x: error [%s] %s\n', e.identifier, e.message); end
try
    y = diff(x,'s',2); fprintf('L3 diff(x,s,2): ok, class %s\n', class(y));
catch e, fprintf('L3 diff(x,s,2): error [%s] %s\n', e.identifier, e.message); end
try
    y = -x; fprintf('L4 -x: ok, class %s\n', class(y));
catch e, fprintf('L4 -x: error [%s] %s\n', e.identifier, e.message); end

% (2) the post-rename behaviour of the legacy lines, through today's name
try
    copvar T1; T1.R.R0 = eye(2); fprintf('R1 T.R.R0 = eye(2): ok, size %s\n', mat2str(size(T1)));
catch e, fprintf('R1 T.R.R0 = eye(2): error [%s] %s\n', e.identifier, e.message); end
try
    copvar T2; T2.I = [0 1]; fprintf('R2 T.I = [0 1]: ok\n');
catch e, fprintf('R2 T.I = [0 1]: error [%s] %s\n', e.identifier, e.message); end
try
    copvar T3; s = pvar('s'); T3.var2 = s; fprintf('R3 T.var2 = s: ok\n');
catch e, fprintf('R3 T.var2 = s: error [%s] %s\n', e.identifier, e.message); end
try
    copvar T4; T4.dim = [1 1; 0 1]; fprintf('R4 T.dim = [1 1;0 1]: ok\n');
catch e, fprintf('R4 T.dim = [1 1;0 1]: error [%s] %s\n', e.identifier, e.message); end
try
    copvar T5; T5.P = 1; fprintf('R5 T.P = 1: ok\n');
catch e, fprintf('R5 T.P = 1: error [%s] %s\n', e.identifier, e.message); end
