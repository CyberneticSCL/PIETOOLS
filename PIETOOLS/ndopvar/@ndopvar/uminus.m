function Pop = uminus(Pop)
% POP = UMINUS(POP) returns the 'ndopvar' object representing the scalar
% product (-1)*P for the PI operator P defined by POP.

for ii=1:numel(Pop.C)
    Pop.C{ii} = -Pop.C{ii};
end

end