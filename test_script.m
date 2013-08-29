% This is a script to manually test the functionality of the evaluation package.
% Note: this is not a diagnostic script. It is just a sandbox of code to test things.


% In comments, A = annotation, E = estimate.
ann_file = '/Users/jordan/Documents/structural_analysis_evaluation/struct_mrx_09_000001_gt.txt';
est_file = '/Users/jordan/Documents/structural_analysis_evaluation/struct_mrx_09_000001_pred.txt';

[a_t a_l] = load_annotation(ann_file,'two_column');
[e_t e_l] = load_annotation(est_file,'two_column');

