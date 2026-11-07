% psz_blocks.m -- what eq_use_psatz=[3;4;5;6] costs with a SEARCHED Lyapunov
% operator, which is the configuration the low-rank session actually runs.
%
% My own psatz reach numbers (0.00 stock -> 0.9999 with the four linear
% generators) were all measured with P PINNED TO THE IDENTITY -- the patch header
% says so explicitly and so does the note -- because the point was to vary the
% negativity certificate alone. The other session runs a SEARCHED P and reports
% a certificate at 0.10*lam* with psatz OFF, which is not a contradiction: with a
% full P search the shipped heavy settings were measured to reach 0.25*lam*, so
% 0.10 is inside that and 0.50 is not.
%
% What is NOT measured anywhere is the block count and reach for searched P WITH
% the linear generators, which is what they need to resize a per-block rank
% ladder. Measured here rather than inferred, since one poslpivar_2d call need
% not contribute exactly one Gram block.
cuadmm_path;
LAMSTAR = 2*pi^2;
CFG = { 'linear4' , [3;4;5;6] };
FR  = [0.10 0.50 0.90];
fprintf('PZ cfg|frac|m|nblk|Ks|nnzAt|rel_b|numerr|prosta|t_s\n');
for c = 1:size(CFG,1)
    for f = FR
        try
            clear stateNameGenerator
            pvar s1 s2 t
            x = pde_var(1,[s1;s2],[0,1;0,1]);
            PIE = initialize(convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+f*LAMSTAR*x;
                   subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0]));
            st = lpisettings('light');
            st.sos_opts.solver='mosek'; st.sos_opts.simplify=false;
            v = CFG{c,2};
            if isempty(v)
                st.settings_2d.eq_use_psatz = [0;0];
            else
                st.settings_2d.eq_use_psatz = v;
                % Shipped settings ship only TWO eq_opts_psatz entries (eq_use_psatz
                % is [0;0]), so entries 3 and 4 must be CLONED from an existing one:
                % assigning only .psatz creates a struct with no `exclude` field and
                % poslpivar_2d dies on it.
                otmpl = st.settings_2d.eq_opts_psatz{1};
                for j = 1:numel(v)
                    o = otmpl;  o.psatz = v(j);
                    st.settings_2d.eq_opts_psatz{j} = o;
                    % FULL eq_deg for the psatz blocks: an offset of -1 was measured
                    % to destroy the certificate outright (0.9999 -> 0).
                    st.settings_2d.eq_deg_psatz{j} = st.settings_2d.eq_deg;
                end
            end
            t0 = tic;
            evalc('[sol,Pop,Qop] = PIETOOLS_stability_2D(PIE,st);');
            tw = toc(t0);
            S = cuadmm_private('sdpshape',sol);
            Atf=[];bf=[];
            for q=1:sol.expr.num, Atf=[Atf,sol.expr.At{q}]; bf=[bf;sol.expr.b{q}]; end
            xr = sol.solinfo.RRx(:);
            rb = norm(full(Atf'*xr-bf))/max(norm(full(bf)),eps);
            I = sol.solinfo.info;
            fprintf('PZ %s|%.2f|%d|%d|%s|%d|%.4e|%d|%s|%.1f\n', ...
                CFG{c,1},f,S.m,S.nblk,mat2str(S.Ks),S.nnzAt,rb, ...
                gf(I,'numerr'),'-',tw);
        catch ME
            fprintf('PZ %s|%.2f|ERR|%s\n',CFG{c,1},f,regexprep(ME.message,'[\t\r\n]+',' '));
        end
    end
end
fprintf('PZDONE\n');

function v = gf(I,f), if isstruct(I)&&isfield(I,f), v=I.(f); else, v=-1; end, end
