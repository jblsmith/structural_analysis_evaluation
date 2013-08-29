function [e_onset, e_label, a_onset, a_label] = align_description_to_ann(e_onset, e_label, a_onset, a_label);

% Get rid of "END" labels:
if length(e_label)==length(e_onset),
    e_label = e_label(1:end-1);
end
if length(a_label)==length(a_onset),
    a_label = a_label(1:end-1);
end

% Pad the beginning:
if min(e_onset) < min(a_onset),
    % Delete appropriate amount from the start of the estimated description:
    cut_in_on = max(find(min(a_onset)>e_onset));
    e_label = e_label(cut_in_on:end);
    e_onset = e_onset(cut_in_on:end);
    e_onset(1) = min(a_onset);
elseif min(e_onset) > min(a_onset),
    % Add padding to the beginning of E:
    e_onset = [min(a_onset); e_onset];
    e_label = [{'unlabeled_begin_e'}; e_label];
end

% Pad the end:
if max(e_onset) > max(a_onset),
    % Delete appropriate amount from the end of the estimated description:
    cut_out_on = min(find(max(a_onset)<e_onset));
    e_label = e_label(1:cut_out_on-1);
    % e_label(end) = {'End'};
    e_onset = e_onset(1:cut_out_on);
    e_onset(end) = max(a_onset);
elseif max(e_onset) < max(a_onset),
    % Add padding to the end of E:
    e_onset = [e_onset; max(a_onset)];
    e_label = [e_label; {'unlabeled_end_e'}];
end

e_onset = e_onset - min(a_onset);
a_onset = a_onset - min(a_onset);