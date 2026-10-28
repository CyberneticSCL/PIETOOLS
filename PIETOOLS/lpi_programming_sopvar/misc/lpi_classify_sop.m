function cl = lpi_classify_sop(prog,obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CL = LPI_CLASSIFY_SOP(PROG,OBJ) names the solver exit of the solved LPI
% program PROG from the flags of prog.solinfo.info: status 'infeasible' if
% pinf = 1, else 'unbounded' if dinf = 1, else 'unknown' if numerr = 2,
% else 'inaccurate' if numerr = 1, else 'optimal'; 'unsolved' when PROG
% carries no solver report. It reads back the value of the objective OBJ
% and marks an 'unknown' exit as divergent unless OBJ is given and its
% value is at most 1e3 in modulus.
%
% INPUT
% - prog:   'struct', an LPI program solved by lpisolve or sossolve; read
%           is the 'struct' prog.solinfo.info, with the optional fields
%           numerr, pinf, dinf (double, default 0), feasratio and cpusec
%           (double, default NaN);
% - obj:    (optional) 'dpvar' or 'polynomial' object, the objective of
%           PROG (a scalar function of its decision variables), or 'char',
%           the name of a decision variable (default []: no objective).
%
% OUTPUT
% - cl:     'struct' with fields
%             status     char, 'optimal', 'inaccurate', 'infeasible',
%                        'unbounded', 'unknown' or 'unsolved';
%             numerr, pinf, dinf, feasratio, cpusec
%                        double, the fields of prog.solinfo.info (the
%                        defaults above when absent; NaN for 'unsolved');
%             obj        double, the value of OBJ (NaN without OBJ, for
%                        'unsolved', or when it cannot be read back);
%             divergent  logical, true for status 'unknown' with obj NaN or
%                        |obj| > divtol;
%             divtol     double, 1e3.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_classify_sop
%
% Copyright (C) 2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 10/08/2026
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and a header for the local function; code
%                unchanged. The deleted header text was the interpretation
%                of a divergent exit, a notes reference and the use the
%                degree loop makes of the status.

% Every field unset; 'unsolved' without a solver report
cl = struct('status','','numerr',NaN,'pinf',NaN,'dinf',NaN,'feasratio',NaN,'obj',NaN, ...
            'divergent',false,'divtol',1e3,'cpusec',NaN);
if ~isfield(prog,'solinfo') || ~isfield(prog.solinfo,'info')
    cl.status = 'unsolved';     return
end
% Solver flags (0 when absent), feasibility ratio and time (NaN when absent)
info = prog.solinfo.info;
g = @(f,dflt) local_get(info,f,dflt);
cl.numerr = g('numerr',0);  cl.pinf = g('pinf',0);  cl.dinf = g('dinf',0);
cl.feasratio = g('feasratio',NaN);  cl.cpusec = g('cpusec',NaN);
% Value of the objective; NaN when it cannot be read back
if nargin>=2 && ~isempty(obj)
    try
        if isa(obj,'dpvar') || isa(obj,'polynomial'),   cl.obj = double(lpigetsol_sop(prog,obj));
        else,                                           cl.obj = double(lpigetsol_sop(prog,obj));
        end
    catch
        cl.obj = NaN;
    end
end
% Status, by precedence pinf, dinf, numerr 2, numerr 1
if cl.pinf==1,          cl.status = 'infeasible';
elseif cl.dinf==1,      cl.status = 'unbounded';
elseif cl.numerr==2,    cl.status = 'unknown';
elseif cl.numerr==1,    cl.status = 'inaccurate';
else,                   cl.status = 'optimal';
end
% An 'unknown' exit is divergent unless |obj| <= divtol
if strcmp(cl.status,'unknown')
    cl.divergent = isnan(cl.obj) || abs(cl.obj)>cl.divtol;
end
end


function v = local_get(s,f,dflt)
% V = LOCAL_GET(S,F,DFLT) returns the field F of the struct S as a double,
% or DFLT when S has no nonempty field F.
% INPUT  s: 'struct'; f: char, the field name; dflt: double, the default.
% OUTPUT v: double.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)),  v = double(s.(f));    else,   v = dflt;   end
end
