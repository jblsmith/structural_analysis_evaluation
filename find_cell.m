function index = find_cell(cell1, cellarray)
% function index = find_cell(cell, cell_array)
%
% This function accepts a CELL and a CELL_ARRAY. If one of the entries in
% CELL_ARRAY is equal to the contents of the CELL, then the index of the
% entry is returned. Otherwise the result is 0.

index = 0;
for i=1:numel(cellarray),
    if compare_cells(cell1, cellarray(i)),
        index = i;
    end
end