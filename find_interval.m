function index = find_interval(value, range)
% index = find_interval(value, range)
% 
% This function answers the question: given a VALUE and a set of RANGEs
% defined by each range's onsets, which INDEX corresponds to the RANGE that
% contains the VALUE?

index = find(range(2:end)>value & range(1:end-1)<=value);
if size(index,1)==0,
    index = -1;
end

end