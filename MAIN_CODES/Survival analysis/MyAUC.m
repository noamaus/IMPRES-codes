function [ res ] = MyAUC( x,y )
%MYAUC Summary of this function goes here
%   Detailed explanation goes here
res = 0;
N = numel(x);
for a = 2:N
    res = res + y(a-1) * (x(a)-x(a-1));
end

end

