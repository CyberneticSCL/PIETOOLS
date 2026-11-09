% bench_run.m -- Benchmark of mtimes_tensopvar for the case specified in
% the JSON file named by the environment variable BENCH_SPEC (fields as for
% the input spec of bench_build.m, with optional additional fields 'fun',
% 'merge', 'native', 'build_only', 'profile'); the result is written (JSON)
% to the file named by the environment variable BENCH_OUT. Run in a fresh
% MATLAB process per case.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - intop.params are sparse matrices, so nnz counts are those of the
%   params.

sp = fileparts(mfilename('fullpath'));  addpath(sp);
spec = jsondecode(fileread(getenv('BENCH_SPEC')));
if isfield(spec,'fun'),     mt = str2func(spec.fun);    else,   mt = @mtimes_tensopvar;     end
if isfield(spec,'merge'),   setappdata(0,'P2_MERGE',logical(spec.merge));   end
if isfield(spec,'native'),  setappdata(0,'P3_NATIVE',logical(spec.native));  end
outfile = getenv('BENCH_OUT');
pr = System.Diagnostics.Process.GetCurrentProcess();
res = struct('spec',spec);
try
    % Warm-up on a tiny case (loads classes/functions)
    [A0,B0] = bench_build(struct('N',1,'d',1,'q',0,'seed',99));
    mt(A0,B0);    clear A0 B0
    % Build the case
    tb = tic;
    [Aop,Bop,info] = bench_build(spec);
    res.t_build = toc(tb);
    fprintf('BUILD_DONE t_build=%.2f s, kernelsA=%d, intmat sizes A: %s\n',res.t_build,numel(Aop.params),mat2str(cellfun(@(a) size(a,1),Aop.intmat)));
    nnzA = 0;   bytesA = 0;
    for k=1:numel(Aop.params)
        nnzA = nnzA + nnz(Aop.params{k});
        P = Aop.params{k};   w = whos('P');     bytesA = bytesA + w.bytes;
    end
    res.nkernelsA = numel(Aop.params);      res.nnzA = nnzA;      res.MB_A = bytesA/2^20;
    pr.Refresh();   ws0 = double(pr.WorkingSet64)/2^20;     peak0 = double(pr.PeakWorkingSet64)/2^20;
    res.ws_before = ws0;    res.peak_before = peak0;
    % Timed call
    fprintf('CALL_START\n');
    if isfield(spec,'build_only') && spec.build_only,   error('build_only');   end
    t1 = tic;
    Cop = mt(Aop,Bop);
    res.t_call = toc(t1);
    pr.Refresh();   res.peak_after = double(pr.PeakWorkingSet64)/2^20;     res.ws_after = double(pr.WorkingSet64)/2^20;
    % Output size
    nnzC = 0;   bytesC = 0;
    for k=1:numel(Cop.params)
        P = Cop.params{k};
        nnzC = nnzC + nnz(P);
        w = whos('P');      bytesC = bytesC + w.bytes;
    end
    res.nkernelsC = numel(Cop.params);      res.nnzC = nnzC;      res.MB_C = bytesC/2^20;
    res.nintmatC = cellfun(@(a) size(a,1),Cop.intmat);
    % Repeats for the time when cheap
    if res.t_call<3
        tt = res.t_call;
        for r=1:2
            t1 = tic;   mt(Aop,Bop);   tt(end+1) = toc(t1); %#ok<AGROW>
        end
        res.t_min = min(tt);
    else
        res.t_min = res.t_call;
    end
    % Profile (optional)
    if isfield(spec,'profile') && spec.profile
        profile clear;  profile on;
        mt(Aop,Bop);
        profile off;
        st = profile('info');
        ft = st.FunctionTable;
        selfT = zeros(numel(ft),1);
        for i=1:numel(ft)
            ch = 0;
            if ~isempty(ft(i).Children),  ch = sum([ft(i).Children.TotalTime]);  end
            selfT(i) = ft(i).TotalTime - ch;
        end
        [~,ord] = sort(selfT,'descend');
        top = struct('name',{},'self',{},'total',{},'calls',{},'parents',{});
        for i=1:min(25,numel(ord))
            j = ord(i);
            par = {};
            if isfield(ft,'Parents') && ~isempty(ft(j).Parents)
                par = arrayfun(@(p) ft(p.Index).FunctionName,ft(j).Parents,'UniformOutput',false);
                par = par(1:min(4,numel(par)));
            end
            top(end+1) = struct('name',ft(j).FunctionName,'self',selfT(j),'total',ft(j).TotalTime,'calls',ft(j).NumCalls,'parents',{par}); %#ok<AGROW>
        end
        res.profile_top = top;
        % Lines of mtimes_tensopvar (incl. local functions), by time
        lines = struct('fn',{},'line',{},'time',{},'calls',{});
        for j=1:numel(ft)
            if contains(ft(j).FunctionName,'mtimes_tensopvar') && ~contains(ft(j).FunctionName,{'_old','_alt','backup'})
                EL = ft(j).ExecutedLines;
                for r=1:size(EL,1)
                    lines(end+1) = struct('fn',ft(j).FunctionName,'line',EL(r,1),'time',EL(r,3),'calls',EL(r,2)); %#ok<AGROW>
                end
            end
        end
        [~,o2] = sort([lines.time],'descend');
        res.profile_lines = lines(o2(1:min(25,numel(o2))));
    end
    res.ok = true;
catch ME
    res.ok = false;     res.err = ME.message;
    res.stack = arrayfun(@(s) sprintf('%s:%d',s.name,s.line),ME.stack(1:min(6,numel(ME.stack))),'UniformOutput',false);
end
fid = fopen(outfile,'w');   fwrite(fid,jsonencode(res,'PrettyPrint',true));   fclose(fid);
if res.ok
    fprintf('BENCH t_call=%.3f t_min=%.3f peak_before=%.0f peak_after=%.0f MB nnzA=%d nnzC=%d\n',res.t_call,res.t_min,res.peak_before,res.peak_after,res.nnzA,res.nnzC);
else
    fprintf('BENCH FAILED: %s\n',res.err);
end
