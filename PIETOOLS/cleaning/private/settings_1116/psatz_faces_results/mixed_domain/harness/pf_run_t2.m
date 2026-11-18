% PF_RUN_T2  Phase 1, T2: the direct PIE2PDE test by term set (PF_QDIRECT,
% Q at length 2, dw 0, eppos 1e-4) at FRAC of the margin on S3, S8, S17,
% S18, S26; FRAC set by the caller (default 0.5). SCRATCH, MMP, 11/16/2026.
if ~exist('FRAC','var'),    FRAC = 0.5;     end
S = pf_setup(false);
tau43 = 1.034726;   a26 = pi^2/4 - 0.5;
[optl,names] = pf_termsets({'none','facesS','faces0','faces','box','corners','faces_s1','faces_s2','box_s1','box_s1+faces'});
qd1 = {struct('mult',1,'int',[1 1],'out',1,'in',1)};
sys = {'S3_edge1d_heat2d',struct('r',FRAC*pi^2); 'S8_two_edges',struct('r3',FRAC*pi^2); ...
       'S17_transport_edge_diffusion',struct('r2',FRAC*pi^2); 'S18_lib43_interior_delay',struct('tau',FRAC*tau43,'sgn',-1); ...
       'S26_lib42_ode_heat_delay',struct('a',FRAC*a26)};
for i = 1:size(sys,1)
    for j = 1:numel(optl)
        try
            pf_qdirect(sys{i,1},{sys{i,2}},qd1,0,sprintf('t2_%s_%s_%g',sys{i,1}(1:3),names{j},FRAC),1e-4,0,optl{j});
        catch ME
            fprintf('@@ T2 ERROR %s %s: %s\n',sys{i,1},names{j},strrep(ME.message,newline,' '));
        end
    end
end
fprintf('## T2 DONE frac=%g\n',FRAC);
