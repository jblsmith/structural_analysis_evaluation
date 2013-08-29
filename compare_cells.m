function value = compare_cells(cell1, cell2)
% value = compare_cells(cell1, cell2)
% 
% This function compares two cells. If they are equal, the value is 1;
% otherwise the value is 0.


a = cell2mat(cell1);
b = cell2mat(cell2);

value = false;

if size(a)==size(b),
    if a==b,
        value = true;
    end
end

end