function row = pf_rows(name)
% ROW = PF_ROWS(NAME) candidate rows of the settings proposal (lengths; one
% cell entry per number of variables 1, 2, 3 of a space; len_R for R^q):
%   'Rlight'   1-D: d1 2, d2 2, d3 [2;2;2], all (today's light 1-D storage)
%              2-D: d1 [1 1], d2 [2 2], d3 [1 1;2 2;2 2], all
%              3-D: d1 0, d2 2, d3 [1;2;0], uniform (today's N-D entry)
%   'Rlight0'  as Rlight with d1 = 0 in 2-D and patterns separable (today's
%              default plus the mixed patterns)
%   'Rlight2'  as Rlight with d2 = [3 3], d3 [2 2;3 3;3 3] in 2-D, separable
%              (the storage that certified S18; d1 0)
%   'Rheavy'   1-D: d1 2, d2 2, d3 [3;2;2], all (today's heavy 1-D)
%              2-D: d1 [2 2], d2 [3 3], d3 [2 2;3 3;3 3], all
%              3-D: d1 0, d2 2, d3 [1;2;0], uniform
%   'Rdef'     today's DEG-omitted default in every entry (d1 0, uniform;
%              1-D entry a, b at d2 2, d3 [2;2;0])
% SCRATCH, MMP, 11/16/2026.
switch name
    case 'Rlight'
        row = struct('len_R',1,'len_inc',{{2,[1 1],[0 0 0]}},'len_com',{{2,[2 2],[2 2 2]}}, ...
                     'len_sos',{{[2;2;2],[1 1;2 2;2 2],[1 1 1;2 2 2;0 0 0]}},'patterns',{{'all','all','uniform'}});
    case 'Rlight0'
        row = struct('len_R',1,'len_inc',{{2,[0 0],[0 0 0]}},'len_com',{{2,[2 2],[2 2 2]}}, ...
                     'len_sos',{{[2;2;2],[1 1;2 2;0 0],[1 1 1;2 2 2;0 0 0]}},'patterns',{{'all','separable','uniform'}});
    case 'Rlight2'
        row = struct('len_R',1,'len_inc',{{2,[0 0],[0 0 0]}},'len_com',{{2,[3 3],[2 2 2]}}, ...
                     'len_sos',{{[2;2;2],[2 2;3 3;0 0],[1 1 1;2 2 2;0 0 0]}},'patterns',{{'all','separable','uniform'}});
    case 'Rheavy'
        row = struct('len_R',1,'len_inc',{{2,[2 2],[0 0 0]}},'len_com',{{2,[3 3],[2 2 2]}}, ...
                     'len_sos',{{[3;2;2],[2 2;3 3;3 3],[1 1 1;2 2 2;0 0 0]}},'patterns',{{'all','all','uniform'}});
    case 'Rlight_cap1'      % Rlight with the per-direction cap 1 in 2-D (theta_k + s_k <= 1: the mixed monomial of each direction dropped)
        row = pf_rows('Rlight');    row.cap = {Inf,[1 1],[Inf Inf Inf]};
    case 'Rlight_cap2'      % Rlight with the per-direction cap 2 in 2-D (no restriction at these lengths; the control)
        row = pf_rows('Rlight');    row.cap = {Inf,[2 2],[Inf Inf Inf]};
    case 'Rdef'
        row = struct('len_R',1,'len_inc',{{0,[0 0],[0 0 0]}},'len_com',{{2,[2 2],[2 2 2]}}, ...
                     'len_sos',{{[2;2;0],[1 1;2 2;0 0],[1 1 1;2 2 2;0 0 0]}},'patterns',{{'uniform','uniform','uniform'}});
    case 'Rheavy_u'         % today's heavy lengths as a row (08:30, D48): 1-D heavy (d3 [3;2;2]); 2-D d1 0, d2 2, d3 [2 2;2 2] (mult degree 1 -> length 2), uniform, no cap (the stock heavy block keeps 9 of these 16 pairs); 3-D as light
        row = struct('len_R',1,'len_inc',{{2,[0 0],[0 0 0]}},'len_com',{{2,[2 2],[2 2 2]}}, ...
                     'len_sos',{{[3;2;2],[2 2;2 2;0 0],[1 1 1;2 2 2;0 0 0]}},'patterns',{{'all','uniform','uniform'}});
    case 'Rheavy_s'         % as Rheavy_u with the separable patterns in 2-D (the mixed operators at the heavy lengths)
        row = pf_rows('Rheavy_u');  row.patterns{2} = 'separable';
    case 'Rlight_cap3'      % Rlight with the per-direction cap 1 in 3-D (27 of the 64 pairs of the complete-moment block; the stock joint cap 2 keeps 22)
        row = pf_rows('Rlight');    row.cap = {Inf,[Inf Inf],[1 1 1]};
    otherwise
        error('pf_rows:name','unknown row %s',name);
end
end
