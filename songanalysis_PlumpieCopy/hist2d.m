function [Output] = hist2d(RowInput, ColInput, RowEdges, ColEdges)

for i = 2:length(RowEdges),
    Temp = (find((RowInput > RowEdges(i-1)) & (RowInput <= RowEdges(i))));
    for j = 2:length(ColEdges),
        Output(i-1,j-1) = length(find(((ColInput(Temp) > ColEdges(j-1)) & (ColInput(Temp) <= ColEdges(j)))));
    end
end
