function [time_points text_labels] = load_annotation(filename, fileformat)

%   function [timepoints labels] = load_annotation(filename, fileformat)
%
%   This function loads the file with FILENAME and determines the
%   onset TIMEPOINTS of each section, as well as the LABEL assigned
%   to each section. The format can either be two-column or three-column
%   (AKA 'lab', the default). Assuming the sections run from no. 1 to
%   no. 0, this is what they look like:
%
%       'lab': [onset1 offset1 label1;
%               onset2 offset2 label2;
%               ...
%               onset0 offset0 label0]
%
%       'two_column': [onset1 label1;
%                      onset2 label2;
%                      ...
%                      onset0 label0
%                      offset0 'end']
%
%   Specifying the format is important because sometimes whitespace
%   delimitation is inconsistent, so that a two-column annotation can appear
%   three-column.

if nargin<2,
    fileformat = 'lab';
end

switch fileformat
    case 'two_column'
        fid = fopen(filename,'r');
        btimes = textscan(fid,'%f%[^\n]');
        % btimes now contains two columns: the first is the times, the second is the labels.
        fclose(fid);
        if isempty(btimes{1}),
            time_points = [];
            text_labels = [];
        else
            text_labels = btimes{2};
            time_points = btimes{1};
        end
    case 'lab'
        fid = fopen(filename,'r');
        btimes = textscan(fid,'%f%f%[^\n]');
        % btimes now contains three columns: the first two are the times, the third is the labels.
        fclose(fid);
        if isempty(btimes{1}),
            time_points = [];
            text_labels = [];
        else
            % Take all the onset times and add a final row.
            text_labels = [btimes{3}; 'end'];
            time_points = [btimes{1}; btimes{2}(end)];
        end
    end
end