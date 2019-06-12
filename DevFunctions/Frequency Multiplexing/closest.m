function [ idx ] = closest( val, array )
%CLOSEST find index of closest value
%   
tmp = abs(array-val);
[idx, idx] = min(tmp); %index of closest value
end

