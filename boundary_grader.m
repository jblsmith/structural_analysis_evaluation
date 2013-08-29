function [grades which_were_found] = boundary_grader(estimated_tims, correct_tims, thr, warmup_time)

% function [grades which_were_found] = boundary_grader(e_onsets, a_onsets, thr, warmup_time)
%
%    Script to grade the estimated set of boundary times E_ONSETS against
%    the correct set of boundary times A_TIMES, using a threshold THR in seconds.
%
%    Optionally, specify a WARMUP_TIME (in seconds) which will ignore
%    any boundaries lying in the first and last WARMUP_TIME seconds
%    (this assumes that the last specified boundary is the end of the piece).
%    This is useful since predicting that there is a boundary in the first
%    and last 3 seconds of a piece is pretty trivial.
%
%    New to version 2 of this: the script is rewritten so that it
%    makes sure that two estimated boundaries that are both very close to
%    the same annotated boundary are NOT both counted as correct. (Very
%    important!)
%
%    New to version 3: the script now returns WHICH_WERE_FOUND, a vector
%    that indicates, for E_ONSETS(i), whether that boundary was a hit
%    or not. If it is a hit, it gives the time to the boundary. That way
%    you can do an evaluation with the maximum tolerance and then cut down
%    the tolerance when necessary. It also enables rather flexible post-analysis.
%
%    Notes on NaNs: when precision = recall = 0, then f-measure is 0, not NaN.
%    All output is NaN when E_ONSETS or A_ONSETS are empty.

if nargin<3,
    thr = 3;
end

% Only do the warm-up time thing if the parameter is specified. Otherwise, leave this all alone.
fake_mt2c = mean(correct_tims);
fake_mc2t = mean(estimated_tims);
if nargin>3,
    end_time = max([correct_tims; estimated_tims]);
    estimated_tims = estimated_tims(estimated_tims >= warmup_time  & estimated_tims <= end_time-warmup_time);
    correct_tims = correct_tims(correct_tims >= warmup_time  & correct_tims <= end_time-warmup_time);
end

if ~isempty(estimated_tims) & ~isempty(correct_tims),
    
    % Pseudocode:
    % # for each correct time in our list, 
    % #     check if there is an estimated time that matches it
    % #     if yes,
    % #         mark that the estimated time was correct
    % #         delete the correct estimated time from its list
    % #     end
    % #     delete the correct time from our list
    % # end

    ctcopy = correct_tims;
    etcopy = estimated_tims;

    ehits = [];
    xi = 1;
    while xi <= length(ctcopy)  % ~isempty(ctcopy),
        correct_time = ctcopy(xi);
        [y, ind] = min(abs(etcopy-correct_time));
        % This index gives the estimated boundary that is closest to the current annotated boundary, and the distance between them, y.
        if y<thr,
            % save the hit to our list
            ehits = [ehits; y etcopy(ind)];
            % but delete this boundary from our list (so that it does not hit two boundaries)
            etcopy = [etcopy(1:ind-1); etcopy(ind+1:end)];            
        end
        xi=xi+1;
        % ctcopy = ctcopy(2:end);
    end

    prec = size(ehits,1)/length(estimated_tims);
    reca = size(ehits,1)/length(correct_tims);
    
    % Make a vector of which were the correct tims:
    if ~isempty(ehits),
        found_times = ehits(:,2); % or, equivalently: setdiff(estimated_tims,etcopy);
        which_were_found = inf(length(estimated_tims),1);
        for xi=1:length(found_times),
            ind = find(estimated_tims==found_times(xi));
            which_were_found(ind) = ehits(xi,1);
        end
    else
        which_were_found = [];
    end
    
    % % % # # # Do the old method just so we can keep est_score and cor_score.

    % # Grade the result:
    
    % # Find the minimum true-to-guess distance for each true boundary:
    t2c = zeros(length(correct_tims),1);
    for xi=1:length(t2c),
        t2c(xi) = min(abs(estimated_tims-correct_tims(xi)));
    end
    mt2c = median(t2c);

    % # Find the minimum guess-to-true distance for each guess:
    c2t = zeros(length(estimated_tims),1);
    for xj=1:length(c2t),
        c2t(xj) = min(abs(correct_tims-estimated_tims(xj)));
    end
    mc2t = median(c2t);

    % # Precision: what fraction of the guesses are less than the threshold from the truth?
    % prec = length(find(c2t<=thr))/length(c2t);
    % # Recall: what fraction of the correct boundaries are less than the threshold from a guess?
    % reca = length(find(t2c<=thr))/length(t2c);
    % # F-measure
    fmea = 2 * prec * reca/(prec+reca);
    if isnan(fmea),
        fmea = 0;
    end
else
    prec = NaN;
    reca = NaN;
    fmea = NaN;
    mt2c = NaN;
    mc2t = NaN;
end

grades = [prec reca fmea mt2c mc2t];