classdef pdeparams
    properties
        n = struct('nx',0,'nw',0,'nu',0,'ny',0,'nz',0,'nv',0,'nr',0,'n_pde',[0]);
        dom = [0,1];
        PDE = struct('A',struct(),'Bpv',[],'Bpb',struct(),'Crp',struct(),'Drv',[],'Drb', struct());
        BC = struct('Ebb',struct(),'Ebp',struct(),'Ebv',[]);
        ODE = struct('Cv',[],'Dvw',[],'Dvu',[],'A',[],'Bxw',[],'Bxu',[],'Bxr',[],'Cz',[],'Dzw',[],'Dzu',[],'Dzr',[],'Cy',[],'Dyw',[],'Dyu',[],'Dyr',[]);
    end
end