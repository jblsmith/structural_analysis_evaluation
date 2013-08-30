function [m, d_ae] = directional_hamming_distance(e_onset, a_onset)

% Script to calculate directional hamming distance between two sets of boundaries.
% To get d_ea and f, just call the variables in reverse order.

% What we do is: for each segment of the annotation, we pick the segment of the
% estimation that overlaps it the most. The portion of the annotated segment that
% is not overlapped is added to the total distance, which is normalized by the 
% length of the song.


d_ae = zeros(length(e_onset-1),1);
for i=1:length(e_onset)-1,
    % Find the overlap of each estimated segment with the annotated one.
    tmp_ols = zeros(length(a_onset)-1,1);
    for j=1:length(tmp_ols),
        tmp_ols(j) = degree_of_overlap(e_onset(i),e_onset(i+1),a_onset(j),a_onset(j+1));
    end
    % Which overlap is the best?
    ol_amount = max(tmp_ols);
    % Calculate the remaining length of the annotated segment and add to
    % the result:
    d_ae(i) = e_onset(i+1) - e_onset(i) - ol_amount;
end

% Directional Hamming distance:
d_ae = sum(d_ae);
% Missed boundaries:
m = d_ae/max(e_onset);