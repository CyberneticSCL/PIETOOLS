% search2d.m -- find a 2-D settings default tuned to cuADMM rather than to an
% interior-point method, WITHOUT giving up the stability threshold.
%
% Quality axis   : the largest lam that the configuration can certify, as a
%                  fraction of lam* = 2*pi^2 = 19.7392 (2-D Dirichlet Laplacian
%                  on the unit square).  A cheap configuration that cannot
%                  prove anything is worthless -- this is the sep=1 lesson.
% Cost axis      : sum(N^3) (cuADMM's dense eig per block per iteration) and
%                  nnz(At) (its memory), NOT m.  m is reported for contrast.
%
% Two facts from the shape map drive the candidate list:
%   - Dupx=Dupy=Dup2=3 is HARDCODED in both light and heavy, and dominates
%     sum(N^3): Dup 3->1 was 98x cheaper at light.
%   - light and heavy differ in the d2 cross-degrees AND in LF_opts.sep
%     (light all-ones, heavy all-zeros).  Since heavy certifies lam=2 and light
%     does not, LF separability is a prime suspect -- and turning it off costs
%     only +24% m with sum(N^3) UNCHANGED, against heavy's ~11x.
%
% Walks lam upward and stops at the first failure, so configurations that fail
% early cost little.

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
LOG  = fullfile(cuadmm_outdir(),'search2d_rows.txt');
fid  = fopen(LOG,'a');

LAMSTAR = 2*pi^2;
FRACS   = [0.10 0.25 0.50 0.75 0.90];

C = { 'light(base)'          , @() base('light')
      'light+LFsep0'         , @() sep0(base('light'))
      'light+LFsep0+Dup2'    , @() dup(sep0(base('light')),2)
      'light+LFsep0+Dup1'    , @() dup(sep0(base('light')),1)
      'light+Dup2'           , @() dup(base('light'),2)
      'heavy(ref)'           , @() base('heavy')
      'heavy+Dup2'           , @() dup(base('heavy'),2)
      'heavy+Dup1'           , @() dup(base('heavy'),1) };

fprintf('S2 label|m|nblk|Nmax|eigcost|nnzAt|frac|lam|cert|rel_b|feasratio|numerr|t\n');
for c = 1:size(C,1)
    lab = C{c,1};
    shape_done = false;  S = [];
    best = 0;
    for f = FRACS
        lam = f*LAMSTAR;
        try
            PIE = mk2d(lam);
            st  = C{c,2}();
            t0  = tic;
            evalc('prog = PIETOOLS_stability_2D(PIE,st);');
            tw  = toc(t0);
            if ~shape_done, S = cuadmm_private('sdpshape',prog); shape_done = true; end
            I  = prog.solinfo.info;
            xv = prog.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6) || norm(xv)<=1e-12;
            cert = (I.numerr==0) && (I.feasratio>0.5) && (rb<1e-5) && ~triv;
            row = sprintf('S2 %s|%d|%d|%d|%g|%d|%.2f|%.4f|%d|%.3e|%+.4f|%d|%.1f\n', ...
                lab,S.m,S.nblk,S.Nmax,S.eigcost,S.nnzAt,f,lam,cert,rb, ...
                I.feasratio,I.numerr,tw);
            fprintf('%s',row);  fprintf(fid,'%s',row);
            if cert, best = f; else, break; end       % stop at first failure
        catch ME
            row = sprintf('S2 %s|ERR|%.2f|%s\n',lab,f,strrep(ME.message,newline,' '));
            fprintf('%s',row);  fprintf(fid,'%s',row);
            break
        end
    end
    row = sprintf('S2REACH %s|reach=%.2f|eigcost=%g|nnz=%d|Nmax=%d|m=%d\n', ...
                  lab,best,S.eigcost,S.nnzAt,S.Nmax,S.m);
    fprintf('%s',row);  fprintf(fid,'%s',row);
end
fclose(fid);
fprintf('S2DONE\n');


function st = base(pre)
st = lpisettings(pre);
% eppos raised off the noise floor, per the maintainer's correction: at 1e-6
% the quantity deciding acceptance sits 2.6e-11 above roundoff.
st.settings_2d.eppos = 1e-2*[1;1;1;1];
st.eppos = 1e-2;  st.eppos2 = 1e-2;
% Solver PINNED.  sossolve defaults to the FIRST solver on the path and
% mosekopt is first in its list, so an unset solver silently selects Mosek --
% which is what produced the earlier 2-D numbers I had attributed to SeDuMi.
st.sos_opts.solver = 'mosek';
end

function st = sep0(st)
st.settings_2d.LF_opts.sep = zeros(6,1);
end

function st = dup(st,Dup)
s2 = st.settings_2d;
dx = s2.LF_deg.dx;  dy = s2.LF_deg.dy;  d2 = s2.LF_deg.d2;
s2.eq_deg.dx = {Dup+dx{1};  Dup+dx{2};  Dup+dx{3}};
s2.eq_deg.dy = {Dup+dy{1},  Dup+dy{2},  Dup+dy{3}};
s2.eq_deg.d2 = cellfun(@(c) Dup+c, d2, 'UniformOutput', false);
if isfield(s2,'eq_deg_psatz')
    for j = 1:numel(s2.eq_deg_psatz)
        s2.eq_deg_psatz{j}.dx = {Dup-1+dx{1};  Dup-1+dx{2};  Dup-1+dx{3}};
        s2.eq_deg_psatz{j}.dy = {Dup-1+dy{1},  Dup-1+dy{2},  Dup-1+dy{3}};
        s2.eq_deg_psatz{j}.d2 = cellfun(@(c) Dup-1+c, d2, 'UniformOutput', false);
    end
end
st.settings_2d = s2;
end

function PIE = mk2d(lam)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);
end
