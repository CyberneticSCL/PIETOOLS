function cl = lpi_classify_sop(prog,obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CL = LPI_CLASSIFY_SOP(PROG[,OBJ]) reads the solver exit of a solved LPI
% program and names it. PROG.solinfo.info is what 'sossolve' keeps of
% SeDuMi or MOSEK (numerr, pinf, dinf, feasratio); OBJ is the objective
% decision variable (a 'dpvar' or name), read back when given.
%
% CL.status is one of
%   'optimal'     numerr 0, no infeasibility flag;
%   'inaccurate'  numerr 1 (the solver reached reduced accuracy);
%   'infeasible'  pinf 1 (a primal infeasibility certificate);
%   'unbounded'   dinf 1;
%   'unknown'     numerr 2: MOSEK status UNKNOWN. With an objective that has
%                 run away (|value| above CL.divtol, default 1e3 times the
%                 unit) CL.divergent is true: the trajectory of a program
%                 that is infeasible at every finite value of the objective
%                 without an exact Farkas witness, the fixed-degree
%                 alternative of the duality note (sopvar_lift_notes Sec. 10).
%                 The degree loop raises the slack weight on this status.
% CL also carries numerr, pinf, dinf, feasratio, obj (NaN without one) and
% cpusec.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cl = struct('status','','numerr',NaN,'pinf',NaN,'dinf',NaN,'feasratio',NaN,'obj',NaN, ...
            'divergent',false,'divtol',1e3,'cpusec',NaN);
if ~isfield(prog,'solinfo') || ~isfield(prog.solinfo,'info')
    cl.status = 'unsolved';     return
end
info = prog.solinfo.info;
g = @(f,dflt) local_get(info,f,dflt);
cl.numerr = g('numerr',0);  cl.pinf = g('pinf',0);  cl.dinf = g('dinf',0);
cl.feasratio = g('feasratio',NaN);  cl.cpusec = g('cpusec',NaN);
if nargin>=2 && ~isempty(obj)
    try
        if isa(obj,'dpvar') || isa(obj,'polynomial'),   cl.obj = double(lpigetsol_sop(prog,obj));
        else,                                           cl.obj = double(lpigetsol_sop(prog,obj));
        end
    catch
        cl.obj = NaN;
    end
end
if cl.pinf==1,          cl.status = 'infeasible';
elseif cl.dinf==1,      cl.status = 'unbounded';
elseif cl.numerr==2,    cl.status = 'unknown';
elseif cl.numerr==1,    cl.status = 'inaccurate';
else,                   cl.status = 'optimal';
end
if strcmp(cl.status,'unknown')
    cl.divergent = isnan(cl.obj) || abs(cl.obj)>cl.divtol;
end
end


function v = local_get(s,f,dflt)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)),  v = double(s.(f));    else,   v = dflt;   end
end
