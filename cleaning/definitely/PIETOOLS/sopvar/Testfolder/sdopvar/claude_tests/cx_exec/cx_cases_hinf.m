function [cases,cases2d] = cx_cases_hinf(which)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [CASES,CASES2D] = CX_CASES_HINF(WHICH) returns the H-infinity family's
% comparison cases.
%
% CASES: 1-D cx_compare cases (kind 'obj': stock objective solve, stock
% program pinned at gamma and bisected, container bisected), on the baseline
% plants io1 (gain executives) and syn1 (synthesis executives), light
% settings, plus Hinf_gain at heavy. WHICH (optional) selects by id: a
% cellstr of ids, or 'gain' / 'syn' / 'heavy' for the three groups, so the
% set can be run in pieces under the per-run time budget.
%
% CASES2D: 2-D STRUCTURE-ONLY cases for cx_hinf_struct2d (never solved):
% the three 2-D gain executives on io2 at light, at the fixed gain
% 1.5 x the plant's closed-form L2 gain (bl_b_io2's formula, 0.1711).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = { % id               exec                        setname  plant
      'hinf_io1',        'Hinf_gain',                'light', {'io1'}
      'hinfco_io1',      'Hinf_gain_coercive',       'light', {'io1'}
      'hinfdu_io1',      'Hinf_gain_dual',           'light', {'io1'}
      'hinfduco_io1',    'Hinf_gain_dual_coercive',  'light', {'io1'}
      'ctrl_syn1',       'Hinf_control',             'light', {'syn1'}
      'est_syn1',        'Hinf_estimator',           'light', {'syn1'}
      'hinf_io1_heavy',  'Hinf_gain',                'heavy', {'io1'} };
grp = {'gain','gain','gain','gain','syn','syn','heavy'};
if nargin<1 || isempty(which),  sel = true(1,size(C,1));
elseif ischar(which),           sel = strcmp(grp,which);
else,                           sel = ismember(C(:,1)',which);
end
cases = struct('id',C(sel,1),'exec',C(sel,2),'setname',C(sel,3),'solver','mosek', ...
               'kind','obj','plant',C(sel,4),'thresh',false);

% Closed-form L2 gain of io2 (Ex_2D_ReactionDiffusion_DDDD, as in bl_b_io2).
nu = 1;     rr = 15;    Mm = 10;    Nn = 10;
mu_mn   = nu*pi^2*((2*(1:Mm)'-1).^2 + (2*(1:Nn)-1).^2);
mn_fact = (2*(1:Mm)'-1).*(2*(1:Nn)-1);
gex = (8/pi^2)*sqrt(sum((1./((mu_mn-rr).*mn_fact)).^2,'all'));
E2 = {'Hinf_gain_2D','Hinf_gain_2D_non_coercive','Hinf_gain_dual_2D'};
cases2d = struct('id',strcat(E2,'_io2'),'exec',E2,'setname','light','solver','sedumi', ...
                 'plant',{{'io2'}},'gam',1.5*gex);
end
