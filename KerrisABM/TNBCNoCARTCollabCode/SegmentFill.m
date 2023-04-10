function [AgentRe] = SegmentFill(AgentRe,N1,N2,num)
%SEGMENTFILL adds the points between N1 and N2 to AgentRe
%   Detailed explanation goes here

% N1 = [115,15,10];
% N2 = [113,14,25];

% AgentRe = zeros(500,500,500);

d = N1 - N2; %difference
dir = 0; %the direction to move
max = -1;
for ii = 1:3
    if abs(d(ii))>max
        dir = ii;
        max = abs(d(ii));
    end
end

dd = d./d(dir);
X = zeros(abs(d(dir))+1,3);
X(1,:) = N1;

for jj = 1:abs(d(dir))
    x = jj.*dd;
    x = N1+x;
    x = round(x);
    X(jj+1,:) = x;
end

%disp(X);

for kk = 1:length(X)
    AgentRe(X(kk,1),X(kk,2),X(kk,3)) = num;
    
end

end

