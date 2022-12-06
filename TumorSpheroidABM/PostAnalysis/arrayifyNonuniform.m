function out = arrayifyNonuniform(S,name)

% takes in a structure array S, a field name name and outputs an array of
% size [size(S),N]. N is the longest length of S(i).(name) for
% i=1...numel(S). For now, I will assume that S(i).(name) is always a
% vector (never a matrix with two dims of size greater than 1)

if isvector(S(1).(name))
    out = zeros([numel(S),0]);
    for i = 1:numel(S)
        nt = length(S(i).(name));
        if nt < size(out,2)
            out(i,1:nt) = S(i).(name);
            out(i,nt+1:end) = NaN;
        else
            out(:,end+1:nt) = NaN;
            out(i,:) = S(i).(name);
        end
    end

    out = reshape(out,[size(S),size(out,2)]);
else
    field_size = size(S(1).(name));
    out = zeros([numel(S),0,prod(field_size(2:end))]);

    for i = 1:numel(S)
        nt = size(S(i).(name),1);
        if nt < size(out,2)
            out(i,1:nt,:) = reshape(S(i).(name),nt,[]);
            out(i,nt+1:end,:) = NaN;
        else
            out(:,end+1:nt,:) = NaN;
            out(i,:,:) = reshape(S(i).(name),nt,[]);
        end
    end

    out = reshape(out,[size(S),size(out,2),field_size(2:end)]);

end