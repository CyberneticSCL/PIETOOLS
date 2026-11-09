function regen_face()
% REGEN_FACE  Regenerate the shipped heat face (face_n1.mat) from scratch.
%
% The shipped .mat is not magic: this script rebuilds it with the package's own
% discovery on the n=1 heat equation and OVERWRITES face_n1.mat when the result
% certifies.  Budget ~15-25 min, dominated by the Burer-Monteiro search.
%
% Route note: 'bm' is the route that MEASURED to regenerate this face.  The
% deterministic alternative, 'mintrace' (one full SeDuMi solve + truncation),
% was tried on this exact program and did NOT certify: SeDuMi returned numerr=1
% after ~35 min, and the trace objective crushes the Lyapunov block to noise
% (top eigenvalues ~1e-9 against ~8e-6 on the slack block), so its truncated
% eigenvectors do not span a certifying face at any rank <= 6.  That is a
% statement about this route's solution quality, not a rank floor.

here = fileparts(mfilename('fullpath'));
cd(here);              % 2-D conversion measured to fail from a cluttered cwd
pielr_path();

PIE = nb_rd2d(1);
o = struct('settings',set2d_deg(4,[]),'route','bm');
cert = pielr_certify(PIE,o);

if ~cert.ok
    error(['regen_face: discovery did not certify -- face_n1.mat NOT ' ...
           'overwritten.  See the report above.']);
end

face = struct();
face.V    = cert.face;
face.S    = cert.S;         % ORIGINAL-b units: X_i = V_i*S_i*V_i' solves the
                            % UNNORMALISED program; pielr_certify owns the
                            % conversion to solver units (divides by its nb0)
face.Ns   = cert.Ns;
face.part = cert.part;
face.n    = 1;
face.meta = struct( ...
    'system','nb_rd2d(1): 2-D reaction-diffusion, lam=2, Dirichlet', ...
    'settings','set2d_deg(4,[])  (light 2D degrees, eq = LF + 4, psatz 0)', ...
    'route',cert.route, ...
    'op_rel',cert.op_rel, ...
    'rank',cert.rank, ...
    'gate',cert.gate, ...
    'units',['ORIGINAL-b units: S_i are coefficients of the unnormalised ' ...
             'program; do NOT rescale by any nb0 yourself'], ...
    'date',char(datetime('now','Format','yyyy-MM-dd HH:mm')));

save(fullfile(here,'face_n1.mat'),'face');
fprintf('regen_face: face_n1.mat regenerated (route %s, rank %s, op rel %.4g)\n', ...
    cert.route,mat2str(cert.r),cert.op_rel);
fprintf('REGEN_FACE_DONE\n');
end
