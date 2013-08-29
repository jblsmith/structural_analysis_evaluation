function overlap = degree_of_overlap(a1, a2, b1, b2)
% function amount = degree_of_overlap(a1, a2, b1, b2)
% 
% Given two ranges determined by the scalars (A1, A2) and (B1, B2), what is
% the maximum amount of overlap?

% Four situations:

if a1<=b1 & a2<=b2,
    overlap = a2-b1;
elseif b1<=a1 & b2<=a2,
    overlap = b2-a1;
elseif b1<=a1 & a2<=b2,
    overlap = a2-a1;
elseif a1<=b1 & b2<=a2,
    overlap = b2-b1;
end

overlap = max(0,overlap);