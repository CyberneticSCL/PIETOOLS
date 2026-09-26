function [prog,P] = cx_h2_pos(prog,dm,sp,dom,dd,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = CX_H2_POS(PROG,DM,SP,DOM,DD,OPTIONS) is the container call
% standing in for a 1-D executive's  [prog,P] = poslpivar(prog,dim,dd,options):
% a poscopvar over the spaces SP (component counts DM), with poslpivar's
% degrees DD and options translated.
%
% Translation (poslpivar.m:143-262 defaults and gap filling):
%   dd       filled as poslpivar fills it, then cx_pl2pm; [] means
%            poslpivar's nargin==2 default {1,[1,1,1],[1,1,1]} (L143)
%   psatz    passed through; poscopvar accepts 0/1, as poslpivar does
%   exclude  (2:4) -> per-L2-space include mask over [Z1,Z2,Z3] with the
%            matching degree specs; exclude(1) (the R^n term) has no
%            poscopvar form, so it is an error when an R^n space exists
%   sep      NOT transcribed (error): poscopvar's alpha = 4 pairing is not
%            cross-checked against poslpivar sep=1, and no H2 setting used
%            here sets it
% Why: the H2 executives declare every positive operator this way; one
% translator keeps the mapping identical across the six transcriptions.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6 || isempty(options),    options = struct();     end
if ~isfield(options,'psatz'),       options.psatz = 0;      end
if ~isfield(options,'exclude'),     options.exclude = [0 0 0 0];    end
if ~isfield(options,'sep'),         options.sep = 0;        end
if options.sep
    error('cx_h2_pos:sep','poslpivar sep=1 is not transcribed (poscopvar alpha=4 unverified against it).')
end
dd = fill_degrees(dd);
deg = cx_pl2pm(dd,sp);
exc = logical(options.exclude(:).');
isR = cellfun(@isempty,sp);
if exc(1) && any(isR)
    error('cx_h2_pos:excludeR','exclude(1) (the R^n term) has no poscopvar form.')
end
opt = struct('psatz',options.psatz);
if any(exc(2:4))
    keep = ~exc(2:4);                       % over poslpivar's [Z1,Z2,Z3]
    if ~any(keep),  error('cx_h2_pos:empty','exclude removes every L2 term.'),  end
    inc = cell(1,numel(sp));
    for k = 1:numel(sp)
        if isR(k),  inc{k} = true;          % R^n: the identity basis only
        else,       inc{k} = keep;  deg{k} = deg{k}(keep);
        end
    end
    opt.include = inc;
end
[prog,P] = poscopvar(prog,dm,sp,dom,deg,opt);
end


function d = fill_degrees(d)
% poslpivar.m:143 (default) and L218-250 (gap filling), cell forms only.
if isempty(d),  d = {1,[1,1,1],[1,1,1]};    return,     end
if isnumeric(d)
    d = d(:).';
    if isscalar(d),         d = {d};
    elseif numel(d)==2,     d = {max(d),[d,max(d)],[d,max(d)]};
    else,                   d = {max(d),d,d};
    end
end
if numel(d)==1,     d{2} = [d{1},d{1},d{1}];    d{3} = d{2};    return,     end
if numel(d{2})==1,  d{2} = d{2}*ones(1,3);  elseif numel(d{2})==2,  d{2}(3) = max(d{2});  end
if numel(d)==2,     d{3} = d{2};    return,     end
if numel(d{3})==1,  d{3} = d{3}*ones(1,3);  elseif numel(d{3})==2,  d{3}(3) = max(d{3});  end
end
