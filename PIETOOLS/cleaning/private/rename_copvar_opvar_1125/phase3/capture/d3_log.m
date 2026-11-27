function d3_log(kind,names,vals,note)
% D3_LOG(KIND,NAMES,VALS,NOTE), 11/27/2026: appends one line to the file
% named by the environment variable D3_LOG: the row (D3_ROW), a counter
% within the row, KIND, NAME=digest for each value, and NOTE. Does nothing
% when D3_LOG is unset.
f = getenv('D3_LOG');
if isempty(f),  return,     end
row = getenv('D3_ROW');
persistent seq lastrow
if isempty(seq) || ~strcmp(lastrow,row),    seq = 0;    lastrow = row;  end
seq = seq + 1;
parts = cell(1,numel(vals));
for k = 1:numel(vals),  parts{k} = [names{k},'=',d3_hash(vals{k})];     end
fid = fopen(f,'a');
fprintf(fid,'%s\t%d\t%s\t%s\t%s\n',row,seq,kind,strjoin(parts,' '),note);
fclose(fid);
end
