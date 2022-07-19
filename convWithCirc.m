function out = convWithCirc(sp, win1)
% function out = convWithCirc(sp, win1) - win1 convtrimmed with collumns of sp
% helper function: circular convolution
%
% Written by Andres Grosmark in 2021

edgeSize = ceil(length(win1)/2);
Qk = [zeros(edgeSize, 1); zeros(size(sp, 1), 1) + 1; zeros(edgeSize, 1)];

sp = [sp((end - (edgeSize - 1)):end, :); sp; sp(1:edgeSize, :)];
out = zeros(size(sp));

for I = 1:size(sp, 2)
    out(:, I) = convtrim(sp(:, I), win1)/sum(win1);
end

out = out(Qk == 1, :);

end

function [c] = convtrim(a,b)
% CONVTRIM trimmed convolution
% c = convtrim(a,b) convolves vectors A and B. The resulting
% vector is length LENGTH(a)
% 
% this function is a wrapper for conv - the only difference is the trimming



if (length(a) <= length(b))
  error('convtrim: the length of vector a must be larger than vector b');
end
  
tempC = conv(a,b);
FrontTrim = floor(length(b)/2);

if (mod(length(b),2) ~= 0)
  BackTrim = floor(length(b)/2);
else
  BackTrim = floor(length(b)/2)-1;
end

c = tempC(FrontTrim+1:end-BackTrim);
end