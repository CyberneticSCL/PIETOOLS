function probe4d_mmp(label)
% PROBE4D_MMP(LABEL) every tier on L2[s1,...,s_n], n = 3, 4, 5, through
% POSCOPVAR without DEG (the storage row): basis rows, Gram sizes, the
% number of decision variables and the time, or the error message. LABEL
% names the diary. The repository POSCOPVAR must resolve.
here = fileparts(mfilename('fullpath'));
diary(fullfile(here,[label '_diary.log']));
root = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
if isempty(which('poscopvar')),  run(fullfile(root,'pietools_path_update.m'));  end
w = strrep(which('poscopvar'),'\','/');
assert(strcmpi(w,[root '/lpi_programming_sopvar/poscopvar.m']),'poscopvar resolves to %s',w);
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
tiers = {'stripped','light','custom','heavy','veryheavy','extreme'};
fprintf('--- %s: %s (%s)\n',label,w,datestr(now,'yyyy-mm-dd HH:MM:SS'));
for n = [3 4 5]
    nm = arrayfun(@(k) sprintf('s%d',k),1:n,'UniformOutput',false);
    vn = polynomial(nm(:));     dn = repmat([0 1],n,1);     sp = {nm};
    for k = 1:numel(tiers)
        try
            t0 = tic;
            [p2,~,Qc,bl] = poscopvar(lpiprogram(vn,dn,tiers{k}),1,sp);
            g = zeros(size(bl,1),1);
            for b = 1:numel(g),     g(b) = size(Qc{b,b},1);     end
            fprintf('n=%d %-10s OK   %7.2f s  ndv %8d  rows %s\n',n,tiers{k},toc(t0), ...
                    size(p2.decvartable,1),mat2str([bl g]));
        catch ME
            fprintf('n=%d %-10s ERR  %s\n',n,tiers{k},ME.message);
        end
    end
end
diary off
end
