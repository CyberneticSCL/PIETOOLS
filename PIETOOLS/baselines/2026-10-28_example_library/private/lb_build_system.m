function [PIE, tg, sinfo] = lb_build_system(m, v, doconvert)
% [PIE, TG, SINFO] = LB_BUILD_SYSTEM(M, V, DOCONVERT) declares one library
% example, with its parameter set to V, and converts it to a PIE.
%
% INPUT
%   m:          one entry of lb_meta or of the catalog
%   v:          parameter value for m.setter; [] keeps the library default
%   doconvert:  true - convert; false - do not; 'auto' - convert unless the
%               PDE has 2 or more spatial variables
%
% OUTPUT
%   PIE:    the PIE (opvar or opvar2d based struct), [] if not converted
%   tg:     struct of executive toggles set by the example (lb_toggles)
%   sinfo:  struct with fields 'dim' (spatial dimension of the PDE, 1 for
%           DDE/NDS/DDF), 'kind' ('PDE','DDE','NDS','DDF') and 'size'
%           (char, first column of PIE.T.dim, '' if not converted)
%
% NOTES
%   - PDE: the example function is called as the library switch calls it,
%     feval(name,[index,]GUI=0,params), then initialize_PIETOOLS_PDE as in
%     the library, then convert. The toggles are read from the example
%     source, the evalin('base',...) lines that set them (lb_pde_toggles).
%   - DDE: initialize_PIETOOLS_DDE, convert_PIETOOLS_DDE(DDE,'pie').
%     NDS: initialize_PIETOOLS_NDS, convert_PIETOOLS_NDS(NDS,'pie').
%     DDF: initialize_PIETOOLS_DDF, minimize_PIETOOLS_DDF,
%     convert_PIETOOLS_DDF(DDF,'pie'), the DDF_min path of PIETOOLS_DDF.m.

if nargin<3, doconvert = 'auto'; end
names = lb_toggles();
PIE = [];
sinfo = struct('dim',1,'kind','','size','');
if strcmp(m.lib,'PDE')
    params = {};
    if ~isempty(v),  params = {sprintf(m.setter,v)};  end
    if m.idxarg
        PDE_t = feval(m.fn,str2double(m.id),0,params);
    else
        PDE_t = feval(m.fn,0,params);
    end
    PDE = initialize_PIETOOLS_PDE(PDE_t,true);
    tg = lb_pde_toggles(m.fn, names);
    sinfo.kind = 'PDE';
    sinfo.dim = PDE.dim;
    if isequal(doconvert,true) || (ischar(doconvert) && PDE.dim<=1)
        PIE = convert(PDE,'pie');
    end
else
    here = fileparts(mfilename('fullpath'));
    code = fileread(fullfile(here,'blocks',[m.id,'.m']));
    setter = '';
    if ~isempty(v),  setter = sprintf(m.setter,v);  end
    [sys, kind, tg] = lb_eval_block(code, setter);
    sinfo.kind = kind;
    if ~isequal(doconvert,false)
        switch kind
            case 'DDE'
                sys = initialize_PIETOOLS_DDE(sys);
                PIE = convert_PIETOOLS_DDE(sys,'pie');
            case 'NDS'
                sys = initialize_PIETOOLS_NDS(sys);
                PIE = convert_PIETOOLS_NDS(sys,'pie');
            case 'DDF'
                sys = initialize_PIETOOLS_DDF(sys);
                sys = minimize_PIETOOLS_DDF(sys);
                PIE = convert_PIETOOLS_DDF(sys,'pie');
        end
    end
end
if ~isempty(PIE)
    sinfo.size = mat2str(PIE.T.dim(:,1)');
end
end

function tg = lb_pde_toggles(fn, names)
% Toggles a PDE example sets, read from its source: every uncommented
% line evalin('base','<toggle> = 0|1;'). In all library example files
% these lines are unconditional and unindented (checked 10/28/2026), so
% the source gives the values the call leaves in the base workspace.
% Reading the base workspace back with evalin crashed R2025b.
tg = struct();
for k = 1:numel(names), tg.(names{k}) = 0; end
L = splitlines(fileread(which(fn)));
for k = 1:numel(L)
    t = regexp(L{k},'^\s*evalin\(\s*''base''\s*,\s*''(\w+)\s*=\s*([01])\s*;?\s*''\s*\)','tokens','once');
    if ~isempty(t) && isfield(tg,t{1})
        tg.(t{1}) = str2double(t{2});
    end
end
end
