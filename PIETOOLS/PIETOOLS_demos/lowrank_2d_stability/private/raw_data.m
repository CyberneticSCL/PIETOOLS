function [Atf,bf,Ns,Kf] = raw_data(sos)
    Atf = [];   bf = [];
    for i = 1:sos.expr.num
        Atf = [Atf, sos.expr.At{i}];    %#ok<AGROW>
        bf  = [bf;  sos.expr.b{i}];     %#ok<AGROW>
    end
    Kf = sos.var.idx{1}-1;
    Ns = [];
    for i=1:sos.var.num
        if strcmp(sos.var.type{i},'sos')
            Ns(end+1) = sqrt(sos.var.idx{i+1}-sos.var.idx{i}); %#ok<AGROW>
        end
    end
    for i=1:sos.extravar.num
        Ns(end+1) = sqrt(sos.extravar.idx{i+1}-sos.extravar.idx{i}); %#ok<AGROW>
    end
end
