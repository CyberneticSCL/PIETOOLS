% shapemap1d.m -- how do the PIETOOLS settings move the SDP's SHAPE?
% CPU only, no solves.  The defaults were chosen against an interior-point
% cost model (minimise m, the dense Schur is O(m^2) memory / O(m^3) time).
% cuADMM's cost model is different: memory ~ nnz(At), per-iteration time
% ~ sum(N^3) for a DENSE eig per block, and m is nearly free.  This maps every
% settings axis onto both cost models so the re-tuning is chosen, not guessed.
%
% Columns: m and nnz drive the IPM and cuADMM memory respectively; eigcost
% = sum(N^3) is cuADMM's per-iteration work; nblk vs Nmax says whether that
% work can use the 15 concurrent eig streams or serialises on one big block.

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');
global CENSUS_PROG

% --- the three executives worth separating: feasibility, cheap objective,
%     and the objective form that carries a large indefinite free block.
EX = {'PIETOOLS_PIE2PDEstability', @sysA
      'PIETOOLS_Hinf_gain_coercive', @sysB
      'PIETOOLS_Hinf_gain'         , @sysB};

% --- settings variants.  Each entry: label, and a mutator on a base struct.
V = { 'preset:stripped' , @() lpisettings('stripped')
      'preset:light'    , @() lpisettings('light')
      'preset:heavy'    , @() lpisettings('heavy')
      'preset:veryheavy', @() lpisettings('veryheavy')
      'preset:extreme'  , @() lpisettings('extreme')
      'Dup=2'           , @() setDup(lpisettings('light'),2)
      'Dup=3'           , @() setDup(lpisettings('light'),3)
      'sep=1'           , @() setsep(lpisettings('light'),1)
      'excl1=[0 0 1 1]' , @() setexcl(lpisettings('light'),[0 0 1 1])
      'excl1=[1 1 0 0]' , @() setexcl(lpisettings('light'),[1 1 0 0])
      'psatzLF(ovr1=0)' , @() setfld(lpisettings('light'),'override1',0)
      'noPsatzD(ovr2=1)', @() setfld(lpisettings('light'),'override2',1)
      'sosineq_on=1'    , @() setineq(lpisettings('light')) };

fprintf('SM header exec,variant,m,nvar,Kf,nblk,Nmax,Ns,eigcost,nnzAt,svec,normb_raw\n');
for e = 1:size(EX,1)
    for v = 1:size(V,1)
        CENSUS_PROG = [];
        try
            PIE = EX{e,2}();
            st  = V{v,2}();
            f   = str2func(EX{e,1});
            evalc('f(PIE,st);');
            S = cuadmm_private('sdpshape',CENSUS_PROG);
            fprintf('SM %s,%s,%d,%d,%d,%d,%d,[%s],%g,%d,%d,%.4e\n', ...
                EX{e,1}, V{v,1}, S.m, S.nvar, S.Kf, S.nblk, S.Nmax, ...
                strtrim(num2str(S.Ks)), S.eigcost, S.nnzAt, S.svec_len, S.normb);
        catch ME
            fprintf('SM %s,%s,FAILED,%s\n', EX{e,1}, V{v,1}, ...
                    strrep(ME.message,newline,' '));
        end
    end
end
fprintf('SMDONE\n');


% ---- settings mutators ---------------------------------------------------
cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function st = setDup(st,Dup)
% Rebuild the NEGATIVITY degrees dd2/dd3 for a different degree bump.  Dup is
% a local inside settings_PIETOOLS_light and is baked into dd2/dd3, so it can
% only be varied by regenerating them.  light has n_order = (1,1,1), n4 = 2.
n1=1; n2=1; n3=1; n4=n2+n3;
st.dd2 = {n1+Dup,   [n2+Dup-1, n3+Dup,   n4+Dup  ], [n2+Dup-1, n3+Dup,   n4+Dup  ]};
st.dd3 = {n1+Dup-1, [n2+Dup-2, n3+Dup-1, n4+Dup-1], [n2+Dup-2, n3+Dup-1, n4+Dup-1]};
end

function st = setsep(st,v)
st.options1.sep = v;  st.options12.sep = v;
if isfield(st,'options2'), st.options2.sep = v; end
if isfield(st,'options3'), st.options3.sep = v; end
end

function st = setexcl(st,e)
st.options1.exclude = e;
end

function st = setfld(st,f,v)
st.(f) = v;
end

function st = setineq(st)
% The lpi_ineq path is only wired up when sosineq_on was true at CONSTRUCTION,
% so opts must be supplied by hand here.
st.sosineq_on = 1;
st.opts.psatz = 0;
st.opts.pure  = 1;
end

% ---- test systems --------------------------------------------------------
function PIE = sysA()
pvar s t
x = pde_var('state',1,s,[0,1]);
sys = [diff(x,t,1) == diff(x,s,2) + 2*x;
       subs(x,s,0)==0; subs(x,s,1)==0];
PIE = convert(sys);
end

function PIE = sysB()
pvar s t
x = pde_var('state',1,s,[0,1]);
w = pde_var('input',1);
z = pde_var('output',1);
sys = [diff(x,t,1) == diff(x,s,2) + 2*x + s*w;
       z == int(x,s,[0,1]);
       subs(x,s,0)==0; subs(x,s,1)==0];
PIE = convert(sys);
end
