function rounded = roundNearest(num, res)
%ROUND NEAREST rounds to nearest number according to resolution specified
%   num : input array or number
%   res : resolution desired (can be an array but it must have same number
%   of elements as num
%
%   OUTPUT
%   rounded : rounded number

invres = 1./res;
rounded = round(num .* invres)./invres;

end

