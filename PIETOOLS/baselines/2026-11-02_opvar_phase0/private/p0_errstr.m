function s = p0_errstr(ME, root)
% S = P0_ERRSTR(ME,ROOT): the error message, its identifier and the first
% stack frame inside ROOT, on one line.
s = strtok(ME.message, newline);
if ~isempty(ME.identifier),     s = [ME.identifier,': ',s];     end
for k = 1:numel(ME.stack)
    if nargin<2 || startsWith(ME.stack(k).file, root)
        [~,f] = fileparts(ME.stack(k).file);
        s = sprintf('%s @ %s:%d', s, f, ME.stack(k).line);
        break
    end
end
s = strrep(s, newline, ' ');
end
