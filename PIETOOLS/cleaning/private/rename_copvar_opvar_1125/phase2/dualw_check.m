function dualw_check(outfile)
% Each H-infinity gain form, each H2 norm form and the two synthesis
% executives on a 1-D reaction-diffusion PIE (lam = pi^2/2, Dirichlet),
% once with a distributed disturbance w(s,t) and once with a
% finite-dimensional w(t) entering through s(1-s); u(t) enters uniformly,
% z = [int x; u]. Tier 'light'. One line per run to OUTFILE.
w1 = which('-all','opvar');  w2 = which('-all','convert');
% convert: two class methods (@pde_struct, legacy @sys), both in this tree
if numel(w1)~=1 || numel(w2)~=2 || ~all(startsWith([w1;w2],pwd))
    error('dualw:path','opvar or convert resolves %d / %d times',numel(w1),numel(w2));
end
fprintf('opvar: %s\n',w1{1});
fid = fopen(outfile,'w');
fprintf(fid,'case\trun\tstatus\tcertified\tgam\tseconds\terror\n');
pvar s t
lam = pi^2/2;
for kase = 1:2
    x = pde_var(s,[0,1]);
    if kase==1
        w = pde_var('in',s,[0,1]);  Bw = 1;           cname = 'distributed';
    else
        w = pde_var('in');          Bw = s*(1-s);     cname = 'finite';
    end
    u = pde_var('control',1);
    z = pde_var('output',2);
    eqs = [diff(x,t) == diff(x,s,2) + lam*x + Bw*w + u;
           z == [int(x,s,[0,1]); u];
           subs(x,s,0) == 0;
           subs(x,s,1) == 0];
    PIE = convert(initialize(eqs));
    st = lpisettings_sop('light');
    runs = {'Hinf_gain Q',   @PIETOOLS_Hinf_gain_sop, {'Q'},   4;
            'Hinf_gain P',   @PIETOOLS_Hinf_gain_sop, {'P'},   4;
            'Hinf_gain Qd',  @PIETOOLS_Hinf_gain_sop, {'Qd'},  4;
            'Hinf_gain Pd',  @PIETOOLS_Hinf_gain_sop, {'Pd'},  4;
            'H2_norm c',     @PIETOOLS_H2_norm_sop,   {'c'},   6;
            'H2_norm o',     @PIETOOLS_H2_norm_sop,   {'o'},   6;
            'H2_norm cco',   @PIETOOLS_H2_norm_sop,   {'cco'}, 6;
            'H2_norm oco',   @PIETOOLS_H2_norm_sop,   {'oco'}, 6;
            'Hinf_control',  @PIETOOLS_Hinf_control_sop, {},   6;
            'H2_control',    @PIETOOLS_H2_control_sop,   {},   7};
    for k = 1:size(runs,1)
        out = cell(1,runs{k,4});  msg = '';  tic;
        try
            evalc('[out{:}] = runs{k,2}(PIE,st,runs{k,3}{:});');
            info = out{end};  stat = info.status;  cert = info.certified;
            gam = out{3};  if isempty(gam), gam = NaN; end
        catch err
            stat = 'ERROR';  cert = false;  gam = NaN;  msg = err.message;
        end
        sec = toc;
        fprintf(fid,'%s\t%s\t%s\t%d\t%.8g\t%.1f\t%s\n',cname,runs{k,1},stat,cert,double(gam),sec,strrep(msg,newline,' '));
        fprintf('%s\t%s\t%s\t%d\t%.8g\t%.1f\t%s\n',cname,runs{k,1},stat,cert,double(gam),sec,strrep(msg,newline,' '));
    end
end
fclose(fid);
end
