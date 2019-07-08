function output = numPossibleTriplet(peaks,cap)
%numPossibleTriplet find number of pairs that have sum less than cap
%   peaks : array of peak positions
%   cap   : maximum sum
pair_comb = combnk(peaks,2);
sum_peak = sum(pair_comb,2);
output = sum(sum_peak <= cap);

end

