function [results vector_lab vector_seg] = compare_structures(e_onset, e_label, a_onset, a_label)

% function results = compare_structures(e_onset, e_label, a_onset, a_label)
% 
% This function takes in two sets of annotations. The onsets vector should be
% a single column. The label vector should be one element shorter than the
% onsets vector (so that the span between successive onsets is matched to a
% single label); if it is the same length, the final label will be assumed to
% be an 'end' tag and will be ignored.
% 
% This function interprets the second description as the annotation and
% the first as the estimation being graded.
% 
% Many evaluation metrics are calculated and all provided in the structure
% RESULTS as RESULTS.name, where 'name' is the name of the metric.

% Save the original values:
a = a_onset;
b = a_label;
c = e_onset;
d = e_label;

% PART 1: boundary precision, recall and f-measure; median claim2true and median true2claim.
res_1 = boundary_grader(e_onset, a_onset, 0.5, 5);
res_6 = boundary_grader(e_onset, a_onset, 3, 5);
b_p1 = res_1(1); b_r1 = res_1(2); b_f1 = res_1(3); mt2c = res_1(4); mc2t = res_1(5);
b_p6 = res_6(1); b_r6 = res_6(2); b_f6 = res_6(3);


%% Loading and data entry

% Pad the segmentation with "unlabeled_begin" and "unlabeled_end" to
% compensate for the fact that some annotations may refer to parts that lie
% outside of the other annotation.

% [e_onset, e_label, a_onset, a_label] = align_descriptions(e_onset, e_label, a_onset, a_label);
[e_onset, e_label, a_onset, a_label] = align_description_to_ann(e_onset, e_label, a_onset, a_label);

% OK, so after all the junky in/out jazz, we should have these things:
% A_ONSET, E_ONSET
% A_LABEL, E_LABEL
% Where each segment is given a particular onset time and label, matched up
% in the vectors. The onset vectors will be one cell longer than the label
% vectors, since the end of the song is supposed to be marked and not given
% a label.

% We can check that this is so:
check_sum = length(a_onset)==1+length(a_label) && length(e_onset)==1+length(e_label);
if check_sum~=1,
    fprintf('Annotations loaded incorrectly: %i\n',check_sum);
end

% Since they were altered, we include the updated onsets and labelling in
% the result:
results.a_onset = a_onset;
results.e_onset = e_onset;
results.a_label = a_label;
results.e_label = e_label;

%% Real work

% Create a common set of boundaries to partition the segments for efficient
% analysis:
bounds = unique([e_onset; a_onset]);
N = max(bounds);

% A/E_VOCAB are string vocabularies for each segmentation.
a_vocab = unique(a_label);
e_vocab = unique(e_label);
% N_A/E = number of states in the annotated and estimated segmentations.
N_a = length(a_vocab);
N_e = length(e_vocab);
% A/E_NVOCAB are corresponding numeric vocabularies for each segmentation.
% a_nvocab = 1:N_a;
% e_nvocab = 1:N_e;

% A_STATES has the numeric segment labels for each segment, expressed
% with reference to the vector A_NVOCAB, or expressed as indexes into the
% reference vector A_VOCAB.
a_states = zeros(length(bounds)-1,1);
e_states = zeros(length(bounds)-1,1);
for i=1:length(a_states),
    % What is the label of the segment represented?
    tmp_a = find_interval(bounds(i),a_onset);
    tmp_e = find_interval(bounds(i),e_onset);
    if tmp_a>0,
        a_states(i) = tmp_a;
    end
    if tmp_e>0,
        e_states(i) = tmp_e;
    end
end

% Let A/EB_LABELS represent the label of each segment with respect to the
% collected segmentation BOUNDS.
ab_labels = zeros(size(a_states));
eb_labels = zeros(size(e_states));
for i=1:length(ab_labels),
    if a_states(i)==0,
        ab_labels(i)=0;
    else
        ab_labels(i) = find_cell(a_label(a_states(i)), a_vocab);
    end
    if e_states(i)==0,
        eb_labels(i)=0;
    else
        eb_labels(i) = find_cell(e_label(e_states(i)), e_vocab);
    end
end

% N_IA and N_JE give the total number of frames that belong to the state I
% or J in the ground truth and estimated segmentations.
% N_IJ is a matrix thats N_A x N_E that gives at entry (i,j) the number of
% frames that correspond to both label i in the annotation and label j in
% the estimation.
n_ia = zeros(N_a,1);
n_je = zeros(N_e,1);
n_ij = zeros(N_a,N_e);
bounds_lengths = diff(bounds);
for i=1:length(bounds_lengths),
    n_ia(ab_labels(i)) = n_ia(ab_labels(i)) + bounds_lengths(i);
    n_je(eb_labels(i)) = n_je(eb_labels(i)) + bounds_lengths(i);
    n_ij(ab_labels(i),eb_labels(i)) = n_ij(ab_labels(i),eb_labels(i)) + bounds_lengths(i);
end
    
%% 2.1 Common notations

% Total number of frames:
results.N = N;
% Number of states in the annotated segmentation:
results.N_a = N_a;
% Number of states in the estimated segmentation:
results.N_e = N_e;
% Total number of frames belonging to the state i in the ground truth:
results.n_ia = n_ia;
% Total number of frames belong to the state j  in the estimation:
results.n_je = n_je;
% Number of frames that simultaneously belong to the state i in the
% annotated segmentation and to the state j in the estimated one.
results.n_ij = n_ij;

%% 2.2 Pairwise precision, recall, and F-measure

% We consider every pair of frames in the song, and ask whether they are
% given the same label in the annotation (belonging to the set MA) and the
% same label in the estimation (belonging to the set ME). For pairwise
% precision, recall, and F-measure, we consider the INTERSECTION of the
% sets MA and ME. (We also calculate the EXCLUSION set, and the two other
% sets, to ensure that we have weighted things properly.) Since we are
% making this calculation in a segment-wise fashion instead of pairwise, we
% need to calculate the WEIGHT of each segment as well, based on the length
% of the segments BOUNDS_LENGTHS. (We also calculate the EXCLUSION set, and
% the two other sets, to ensure that we have weighted things properly.)
Me = 0;
Ma = 0;
bounds_lengths = diff(bounds);
intersection = 0;
exclusion = 0;
Me_only =0;
Ma_only=0;

for i=1:length(bounds)-1,
    for j=1:length(bounds)-1,
        tmp_a = ab_labels(i)==ab_labels(j);
        tmp_e = eb_labels(i)==eb_labels(j);
        weight = bounds_lengths(i)*bounds_lengths(j);
        Ma = Ma + tmp_a*weight;
        Me = Me + tmp_e*weight;
        if tmp_a && tmp_e,
            intersection = intersection + weight;
        elseif ~tmp_a && ~tmp_e,
            exclusion = exclusion + weight;
        elseif ~tmp_a && tmp_e,
            Me_only = Me_only + weight;
        elseif tmp_a && ~tmp_e,
            Ma_only = Ma_only + weight;
        end                
    end
end

pw_p = intersection/(intersection+Me_only);
pw_r = intersection/(intersection+Ma_only);
pw_f = 2*pw_p*pw_r/(pw_p+pw_r);

% Pairwise precision:
results.pw_p = pw_p;
% Pairwise recall:
results.pw_r = pw_r;
% Pairwise f-measure:
results.pw_f = pw_f;


% BONUS MEASURE: The Rand Index
% The Rand Index is the ratio of A+B/(A+B+C+D) where all pairs of elements
% are classified under A if they have the same label in both descriptions
% and B if they have different labels in both descriptions, and under C or
% D if they have the same label in one description but not the other.
% http://en.wikipedia.org/wiki/Rand_index

rand = (exclusion+intersection)/(exclusion+intersection+Me_only+Ma_only);
results.rand = rand;

% 2.3 Purity concept

r_je = zeros(N_e,1);
for j=1:N_e,
    r_je(j) = sum(n_ij(:,j).^2)/n_je(j)^2;
end
r_je(isnan(r_je))=0;
acp = (1/N) * sum(r_je.*n_je);

r_ia = zeros(N_a,1);
for i=1:N_a,
    r_ia(i) = sum(n_ij(i,:).^2)/n_ia(i)^2;
end
r_ia(isnan(r_ia))=0;
asp = (1/N) * sum(r_ia.*n_ia);

K = sqrt(asp*acp);

% Purity:
results.r_je = r_je;
% Average cluster purity:
results.acp = acp;
% Speaker purity:
results.r_ia = r_ia;
% Average speaker purity:
results.asp = asp;
% Final purity score:
results.K = K;


% 2.4 Directional Hamming distance


[m, d_ae] = directional_hamming_distance(e_onset, a_onset);
[f, d_ea] = directional_hamming_distance(a_onset, e_onset);

% Directional Hamming distance:
results.d_ae = d_ae;
% Missed boundaries:
results.m = m;
% Inverse directional Hamming distance:
results.d_ea = d_ea;
% Segment fragmentation:
results.f = f;

%% 2.5 Mutual information

%% 3.0 Conditional entropy

% P_IJ is the joint distribution of labels (i) in the annotation and (j) in
% the estimation.
% P_IA and P_JE are the marginal distribution for the annotated and
% estimated segmentations, respectively.
p_ij = n_ij/sum(sum(n_ij));
p_ia = n_ia/sum(sum(n_ij));
p_je = n_je/sum(sum(n_ij));

% The conditional distributions are given respectively as:
p_ijae = zeros(size(n_ij));
p_jiea = zeros(size(n_ij));
for i=1:size(p_ijae,1),
    p_ijae(i,:) = n_ij(i,:)./transpose(n_je);
end
for j=1:size(p_jiea,2),
    p_jiea(:,j) = n_ij(:,j)./n_ia;
end

% We can now easily calculate the Mutual Information metric I_AE from
% section 2.5:
I_ij = 0;
for i=1:N_a,
    for j=1:N_e,
        I_ij(i,j) = p_ij(i,j)*log(p_ij(i,j)/(p_ia(i)*p_je(j)));
    end
end
I_ij(isnan(I_ij))=0;
I_AE = sum(sum(I_ij));

% The conditional entropies H_EA and H_AE can thus be written as:
tmpEA = p_jiea.*log2(p_jiea);
tmpAE = p_ijae.*log2(p_ijae);
tmpEA(isnan(tmpEA))=0;
tmpAE(isnan(tmpAE))=0;
H_EA = -sum(sum(tmpEA,2).*p_ia);
H_AE = -sum(sum(tmpAE,1).*transpose(p_je));

% The over- and under-segmentation scores then become:
S_o = 1 - H_EA/log2(N_e);
S_u = 1 - H_AE/log2(N_a);

% Mutual information:
results.I_ij = I_ij;
results.I_AE = I_AE;
% Conditional entropy measuring spurious information in estimation:
results.H_EA = H_EA;
% Conditional entropy measuring information missing from estimation:
results.H_AE = H_AE;
% Over-segmentation score:
results.S_o = S_o;
% Under-segmentation score:
results.S_u = S_u;

% Boundary median distance:
results.mt2c = mt2c;
results.mc2t = mc2t;
% Boundary evaluation at 0.5-second tolerance:
results.b_f1 = b_f1;
results.b_p1 = b_p1;
results.b_r1 = b_r1;
% Boundary evaluation at 3-second tolerance:
results.b_f6 = b_f6;
results.b_p6 = b_p6;
results.b_r6 = b_r6;

% Put it all in two vectors:
vector_lab = [pw_f, pw_p, pw_r, K, asp, acp, I_AE, H_EA, H_AE, S_o, S_u, rand];
vector_seg = [mt2c, mc2t, m, f, d_ae, d_ea, b_f1, b_p1, b_r1, b_f6, b_p6, b_r6];
