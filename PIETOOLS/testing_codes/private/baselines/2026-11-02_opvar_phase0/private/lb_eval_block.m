function [lb__sys, lb__kind, lb__tg] = lb_eval_block(lb__code, lb__setter)
% [SYS, KIND, TG] = LB_EVAL_BLOCK(CODE, SETTER) runs one uncommented DDE
% or NDS/DDF library block, then the parameter setter, in a workspace of
% its own.
%
% INPUT
%   code:   char, the block as written by lb_blocks
%   setter: char, code run after the block ('' for none), e.g.
%           'tauN=1.2; DDE.tau(1)=tauN/2; DDE.tau(2)=tauN;'
%
% OUTPUT
%   sys:    the DDE, NDS or DDF struct the block defines
%   kind:   'DDE', 'NDS' or 'DDF'
%   tg:     struct of the executive toggles of the stock flow
%           (stability, stability_dual, Hinf_gain, ...), 0 if not set
%
% NOTES
%   - Local names carry the prefix lb__ so that block variables (n, i, h,
%     gamma, ...) cannot collide with them.

eval(lb__code);
if ~isempty(lb__setter)
    eval(lb__setter);
end
if exist('DDE','var')
    lb__sys = DDE;   lb__kind = 'DDE';
elseif exist('NDS','var')
    lb__sys = NDS;   lb__kind = 'NDS';
elseif exist('DDF','var')
    lb__sys = DDF;   lb__kind = 'DDF';
else
    error('lb_eval_block:nosys','The block defines no DDE, NDS or DDF struct.');
end
lb__names = lb_toggles();
lb__tg = struct();
for lb__k = 1:numel(lb__names)
    lb__val = 0;
    if exist(lb__names{lb__k},'var')
        % Statement form: R2025b crashed (access violation) in eval/evalin
        % calls that return a value.
        eval(['lb__val = double(',lb__names{lb__k},');']);
    end
    lb__tg.(lb__names{lb__k}) = lb__val;
end
end
