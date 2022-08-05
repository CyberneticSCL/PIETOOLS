classdef pie_struct
    properties
        dom = [0,1];
        vars = [pvar('s'), pvar('theta')];
        dim = 1;
        T;
        Tw;
        Tu;
        A;
        B1;
        B2;
        C1
        D11;
        D12;
        C2;
        D21;
        D22;
    end
    properties (hidden)
        signal_dims = [0,0,0,0,0,0]; % store sizes of (x,X,z,y,w,u) here
    end
end