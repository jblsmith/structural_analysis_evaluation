function [e_onset, e_label, a_onset, a_label] = align_descriptions(e_onset, e_label, a_onset, a_label);

% Get rid of "END" labels:
if length(e_label)==length(e_onset),
    e_label = e_label(1:end-1);
end
if length(a_label)==length(a_onset),
    a_label = a_label(1:end-1);
end

% Pad the beginning:
if min(e_onset) < min(a_onset),
    % Add padding to the beginning of A:
    a_onset = [min(e_onset); a_onset];
    a_label = [{'unlabeled_begin_a'}; a_label];
elseif min(e_onset) > min(a_onset),
    % Add padding to the beginning of E:
    e_onset = [min(a_onset); e_onset];
    e_label = [{'unlabeled_begin_e'}; e_label];
end

% Pad the end:
if max(e_onset) > max(a_onset),
    % Add padding to the end of A:
    a_onset = [a_onset; max(e_onset)];
    a_label = [a_label; {'unlabeled_end_a'}];
elseif max(e_onset) < max(a_onset),
    % Add padding to the end of E:
    e_onset = [e_onset; max(a_onset)];
    e_label = [e_label; {'unlabeled_end_e'}];
end
