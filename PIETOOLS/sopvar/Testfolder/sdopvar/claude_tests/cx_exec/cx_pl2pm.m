function deg = cx_pl2pm(dd,sp)
% DEG = CX_PL2PM(DD,SP) maps 1-D 'poslpivar' degrees DD = {d1,[d2],[d3]} to
% 'poscopvar' degrees, one entry per space in SP: an R^n space has the
% identity basis (degree 0 in the integration variable); an L_2 space gets
% poslpivar's Z1, Z2, Z3 with d{1} -> int, d{2}(1) -> int, d{2}(2) -> mult,
% d{2}(3) -> joint, likewise d{3} (poslpivar.m:321-346). This is the mapping
% test_poscopvar_vs_poslpivar verifies and test_copvar_kyp uses.
%
% Initial coding MMP, 09/25/2026
d1 = dd{1};     d2 = dd{2};     d3 = dd{3};
L2 = { struct('int',d1,'mult',0), ...
       struct('int',d2(1),'mult',d2(2),'joint',d2(3)), ...
       struct('int',d3(1),'mult',d3(2),'joint',d3(3)) };
deg = cell(1,numel(sp));
for k = 1:numel(sp)
    if isempty(sp{k}),  deg{k} = struct('int',0);
    else,               deg{k} = L2;
    end
end
end
