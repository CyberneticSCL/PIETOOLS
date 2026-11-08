function [x,y,info] = cuadmm_sedumi(At,b,c,K,pars)
% CUADMM_SEDUMI --- Solve an SDP given in SeDuMi format with the cuADMM
% executable.
%
% [x,y,info] = cuadmm_sedumi(At,b,c,K,pars)
%
% Description: solves
%       min c'*x   s.t.   At'*x = b,   x in K,
% with K.f free variables and K.s positive semidefinite blocks, by
% cuADMM, a first-order (ADMM) SDP solver for CUDA GPUs that reads its
% problem from text files in a folder and writes its solution there. The
% PSD blocks are passed in cuADMM's svec coordinates (upper triangle by
% columns, off-diagonal entries scaled by sqrt(2)); the same map
% symmetrizes the constraint matrices, of which only the symmetric part
% acts on a symmetric X. cuADMM stops when its relative primal residual,
% dual residual and gap reach pars.tol, or after pars.maxit iterations.
% It does not detect infeasibility: info.pinf = info.dinf = 0 always, and
% the returned point is to be judged by its residual.
%
% Inputs:
% At:   n x m sparse matrix, the constraint matrix in SeDuMi's format (one
%       column per equality constraint), n = K.f + sum(K.s.^2);
% b:    m x 1 vector, the right-hand side;
% c:    n x 1 vector, the objective;
% K:    cone structure with the fields f (number of free variables) and s
%       (sizes of the PSD blocks); K.l, K.q and K.r must be absent, empty
%       or zero;
% pars: (optional) struct with the optional fields
%         exe       char, the path of cuadmm_exe as the launcher sees it,
%                   where a leading ~/ is the home folder (default: the
%                   environment variable CUADMM_EXE, else
%                   ~/cuADMM/build/cuadmm_exe, where cuadmm_install.sh
%                   builds it);
%         launcher  'auto' (default), 'wsl' or 'linux': 'auto' runs the
%                   executable under WSL on Windows and directly elsewhere;
%         tol       double, the stopping tolerance (default 1e-6);
%         maxit     double, the iteration limit (default 1e5);
%         timeout   double, seconds after which the run is stopped
%                   (default Inf);
%         workdir   char, the folder for the problem and solution files
%                   (default: a new folder in tempdir);
%         keep      logical, keep that folder (default false: deleted).
%
% Outputs:
% x:    n x 1 vector, the primal point in SeDuMi's coordinates;
% y:    m x 1 vector, the dual point (c - At*y is the dual slack);
% info: struct with fields
%         iter       double, the iterations run;
%         cpusec     double, cuADMM's solve time in seconds;
%         timing     1 x 2 double, cuADMM's setup and solve times;
%         pinf, dinf 0, since cuADMM decides neither infeasibility;
%         numerr     0 when cuADMM met pars.tol, 1 at the iteration limit;
%         feasratio  NaN, since cuADMM computes no feasibility ratio;
%         cuadmm     struct with fields pinf, dinf, relgap (cuADMM's last
%                    relative primal residual, dual residual and gap),
%                    converged (logical) and log (char, its log file).
%
% NOTES:
% The executable must take the arguments
%   <folder>/ <tol> <maxit> <eig_streams> <dump_sol>,
% write X_opt.txt, y_opt.txt, S_opt.txt into the folder when dump_sol is
% 1, and print a line "CUADMM_TIMING init <s> solve <s>". The upstream
% main.cu fixes the tolerance and the iteration limit and writes no
% solution, so it needs these additions: cuadmm_exe.patch, in this folder,
% makes them to cuADMM (github.com/ComputationalRobotics/cuADMM, commit
% fcadc30), and its first lines give the build steps. cuadmm_install.sh,
% in this folder, carries out those steps on Ubuntu or Debian, including
% Ubuntu under WSL. Without the CUADMM_TIMING line in its log the run is
% an error: the executable exits with status 0 even when CUDA fails.
% Under WSL, in the default distribution, the executable runs with
% LD_LIBRARY_PATH=/usr/lib/wsl/lib,
% where the WSL GPU driver is. On Linux, system() passes MATLAB's own
% LD_LIBRARY_PATH, whose libstdc++ can be older than the one the
% executable was built with; when the environment variable
% CUADMM_LD_LIBRARY_PATH is set, the executable runs with that path.
% Cost: the svec map and the files are O(nnz(At)) in time and memory; the
% text files take about 34 bytes per nonzero of At.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  M. Peet
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% 10/30/26 - MMP - initial coding
% 11/07/26 - MMP - Help text: the additions the executable needs ship in
%                  this folder as cuadmm_exe.patch, against cuADMM commit
%                  fcadc30; before, they existed only in local builds.
%                  Code unchanged.
% 11/08/26 - MMP - Without pars.exe and CUADMM_EXE, the executable is
%                  ~/cuADMM/build/cuadmm_exe, where the new
%                  cuadmm_install.sh builds it (was: an error). A leading
%                  ~/ in exe stays unquoted, so the shell expands it;
%                  quoted, it required an absolute path. The error for a
%                  failed run names cuadmm_install.sh.

% Options, with the defaults above
if nargin<5 || isempty(pars),   pars = struct();    end
dflt = struct('exe',getenv('CUADMM_EXE'),'launcher','auto','tol',1e-6,'maxit',1e5, ...
              'timeout',Inf,'workdir','','keep',false);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(pars,fn{i}) || isempty(pars.(fn{i})),   pars.(fn{i}) = dflt.(fn{i});    end
end
% if isempty(pars.exe)                                                      % MMP, 11/08/2026 (was)
%     error('cuadmm_sedumi:exe','No cuADMM executable: set pars.exe or the environment variable CUADMM_EXE.'); % MMP, 11/08/2026 (was)
% end                                                                       % MMP, 11/08/2026 (was)
% Neither set: the executable where cuadmm_install.sh builds it
if isempty(pars.exe),   pars.exe = '~/cuADMM/build/cuadmm_exe';     end     % MMP, 11/08/2026
L = pars.launcher;
if strcmp(L,'auto'),    if ispc, L = 'wsl'; else, L = 'linux'; end,     end

% Cones: free variables and PSD blocks only
Kf = 0;     if isfield(K,'f') && ~isempty(K.f),     Kf = double(K.f);   end
for f = {'l','q','r'}
    if isfield(K,f{1}) && any(K.(f{1})(:)~=0)
        error('cuadmm_sedumi:cone','cuADMM interface: K.%s cones are not supported.',f{1});
    end
end
Ks = [];    if isfield(K,'s') && ~isempty(K.s),     Ks = double(K.s(:)');   end
[n,m] = size(At);
if n ~= Kf+sum(Ks.^2)
    error('cuadmm_sedumi:size','size(At,1) = %d differs from K.f + sum(K.s.^2) = %d.',n,Kf+sum(Ks.^2));
end

% svec map T (vec_len x n): identity on the free variables; per PSD block,
% the upper entry (i,j) of each svec row, plus sqrt(2)/2 on both mirrored
% entries of an off-diagonal row, which symmetrizes a constraint matrix and
% scales by sqrt(2) in one step. Built once from triplets.
r2 = sqrt(2)/2;
I = cell(1,numel(Ks)+1);    J = I;  V = I;
I{1} = (1:Kf)';     J{1} = (1:Kf)';     V{1} = ones(Kf,1);
bs = Kf;    bv = Kf;                    % offsets in svec and in SeDuMi's vec
for k = 1:numel(Ks)
    N = Ks(k);
    idx = find(triu(true(N)));          % column-major upper triangle = svec order
    [ii,jj] = ind2sub([N N],idx);
    p = (1:numel(idx))';    od = (ii~=jj);
    I{k+1} = [bs+p;   bs+p(od)];
    J{k+1} = [bv+idx; bv+sub2ind([N N],jj(od),ii(od))];
    V{k+1} = [~od + od*r2;  repmat(r2,nnz(od),1)];
    bs = bs + N*(N+1)/2;    bv = bv + N^2;
end
vec_len = bs;
T = sparse(cat(1,I{:}),cat(1,J{:}),cat(1,V{:}),vec_len,n);
clear I J V

% Problem folder
wd = pars.workdir;
if isempty(wd),     wd = tempname;     end
if ~exist(wd,'dir'),    mkdir(wd);     end
for f = {'X_opt.txt','y_opt.txt','S_opt.txt'}  % no stale solution may survive
    if exist(fullfile(wd,f{1}),'file'),    delete(fullfile(wd,f{1}));     end
end
% blk.txt: 'u' for the free variables, 's N' per nonempty PSD block
fid = fopen(fullfile(wd,'blk.txt'),'w');
if Kf>0,    fprintf(fid,'u %d\n',Kf);   end
fprintf(fid,'s %d\n',Ks(Ks>0));
fclose(fid);
fid = fopen(fullfile(wd,'con_num.txt'),'w');   fprintf(fid,'%d\n',m);  fclose(fid);
% At.txt: 0-based COO of T*At (vec_len x m), rows ascending, as cuADMM's
% COO-to-CSC conversion assumes; the transpose's column order gives that
% order without a sort, and chunks bound the memory of the text buffer
[cj,ri,v] = find((T*At).');
fid = fopen(fullfile(wd,'At.txt'),'w');
ch = 5e6;
for s0 = 1:ch:numel(v)
    e0 = min(s0+ch-1,numel(v));
    fprintf(fid,'%d %d %.17g\n',[ri(s0:e0)-1, cj(s0:e0)-1, v(s0:e0)]');
end
fclose(fid);
clear cj ri v
write_vec(fullfile(wd,'b.txt'),b);
write_vec(fullfile(wd,'C.txt'),T*c(:));

% Run, wrapped per platform; dump_sol = 1, 15 eigensolver streams
logf = fullfile(wd,'cuadmm.log');
pre = '';   if isfinite(pars.timeout),  pre = sprintf('timeout %d ',max(1,ceil(pars.timeout)));  end
% The executable in single quotes, except a leading ~/, which the shell
% expands to the home folder only unquoted
exq = ['''' pars.exe ''''];                                                 % MMP, 11/08/2026
if strncmp(pars.exe,'~/',2),    exq = ['~/''' pars.exe(3:end) ''''];    end % MMP, 11/08/2026
% core = sprintf('%s''%s'' ''%s/'' %.6g %d 15 1 > ''%s'' 2>&1',pre,pars.exe,upath(L,wd), ... % MMP, 11/08/2026 (was)
core = sprintf('%s%s ''%s/'' %.6g %d 15 1 > ''%s'' 2>&1',pre,exq,upath(L,wd), ...
               pars.tol,round(pars.maxit),upath(L,logf));                   % MMP, 11/08/2026
switch L
    case 'wsl'
        cmd = sprintf('wsl -e bash -c "LD_LIBRARY_PATH=/usr/lib/wsl/lib %s"',core);
    case 'linux'
        lp = getenv('CUADMM_LD_LIBRARY_PATH');
        cmd = core;     if ~isempty(lp),    cmd = sprintf('LD_LIBRARY_PATH=''%s'' %s',lp,core);  end
    otherwise
        error('cuadmm_sedumi:launcher','Unknown launcher ''%s''.',L);
end
system(cmd);

% Success only with a CUADMM_TIMING line; last row of the iteration trace
% ' iter | pinf dinf | pobj dobj relgap | ...'
txt = '';   if exist(logf,'file'),  txt = fileread(logf);   end
tm = regexp(txt,'CUADMM_TIMING init (\S+) solve (\S+)','tokens','once');
if isempty(tm)
    tl = txt(max(1,end-1500):end);
%   error('cuadmm_sedumi:failed','cuADMM returned no solution (failed, or stopped by the timeout); end of %s:\n%s',logf,tl); % MMP, 11/08/2026 (was)
    % The hint names the script that builds a missing executable
    error('cuadmm_sedumi:failed',['cuADMM returned no solution (failed, or stopped by the timeout); end of %s:\n%s\n' ...
          'If %s is missing, cuadmm_install.sh in %s builds it.'],logf,tl,pars.exe,fileparts(mfilename('fullpath'))); % MMP, 11/08/2026
end
tok = regexp(txt,'^\s*(\d+)\s*\|\s*(\S+)\s+(\S+)\s*\|\s*(\S+)\s+(\S+)\s+(\S+)\s*\|','tokens','lineanchors');
tr = nan(1,6);  if ~isempty(tok),   tr = str2double(tok{end});  end
conv = ~isempty(strfind(txt,'Solver ended: converged'));

% Solution back in SeDuMi's coordinates: T'*T is the identity on symmetric
% blocks, so T' maps svec(X) to vec(X)
xs = read_vec(fullfile(wd,'X_opt.txt'));
ys = read_vec(fullfile(wd,'y_opt.txt'));
if numel(xs)~=vec_len || numel(ys)~=m
    error('cuadmm_sedumi:solution','cuADMM wrote %d primal and %d dual entries; expected %d and %d.', ...
          numel(xs),numel(ys),vec_len,m);
end
x = T'*xs;
y = ys;
info = struct('iter',tr(1),'cpusec',str2double(tm{2}),'timing',str2double(tm), ...
              'pinf',0,'dinf',0,'numerr',double(~conv),'feasratio',NaN, ...
              'cuadmm',struct('pinf',tr(2),'dinf',tr(3),'relgap',tr(6),'converged',conv,'log',logf));
if ~pars.keep && isempty(pars.workdir)
    rmdir(wd,'s');     info.cuadmm.log = '';
end
end


function write_vec(fn,v)
% WRITE_VEC(FN,V) writes the nonzeros of the vector V to the file FN as the
% 0-based COO lines "index 0 value" that cuADMM reads for b and C.
[i,~,w] = find(sparse(v(:)));
fid = fopen(fn,'w');    fprintf(fid,'%d 0 %.17g\n',[i(:)-1, w(:)]');   fclose(fid);
end

function v = read_vec(fn)
% V = READ_VEC(FN) reads the column of numbers in the file FN.
fid = fopen(fn,'r');
if fid<0,   error('cuadmm_sedumi:open','Cannot open %s.',fn);   end
v = fscanf(fid,'%f');   fclose(fid);
end

function p = upath(L,p)
% P = UPATH(L,P) converts a Windows path to the WSL mount (C:\a -> /mnt/c/a)
% when the launcher L is 'wsl', and returns P unchanged otherwise.
if strcmp(L,'wsl')
    p = strrep(p,'\','/');
    if numel(p)>=2 && p(2)==':',    p = ['/mnt/' lower(p(1)) p(3:end)];  end
end
end
