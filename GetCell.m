%-------------------------------------------------------------------------%
%  Project       : Mitosis-Detection-Breast-C                             %
%  File          : GetCell.m                                              %
%  Description   : Export the cell sub-image                              %
%  Author        : Monachopoulos Konstantinos                              %
%-------------------------------------------------------------------------%

function [mitotic_cell]=GetCell(UnnormalizedImage,MitoticCoordinates)

% Normalize every Histological Image to make it independent of its colorant value ...
[Inorm1 H1 E1] = normalizeStaining(UnnormalizedImage);

% Threshold the normalized image ...
BW=thresh_by_iterations(Inorm1);

% Compare and hold only the set of pixels that represents the pattern 'octagon' (as it is close to the cell shape)

BWerode = imerode(BW, strel('octagon',3));

% Clear the border of the image at the edges..

BWfinal = imclearborder(BWerode, 6);

cell_cnt=0;
flag=0; % conditional variable that checks is the sub image does fit in a sub-square [51 51]
for csv_lines=1:size(MitoticCoordinates,1)
    
    % aligning of the square form and cut the sub-image cell, not selecting
    % the cells at the corners
    xy_mit=MitoticCoordinates(csv_lines,1:2);
    if (size(rgb2gray(imcrop(Inorm1,[xy_mit-30 50 50]))) == [51 51])
        flag=1;
        cell_cnt=cell_cnt+1;
        mitotic_cell(:,:,cell_cnt) = rgb2gray(imcrop(Inorm1,[xy_mit-30 50 50]));
    end
end

% If all of the mitotic or non mitotic cells in coordinates readed from the csv file
% are at the edge of the frame, then return all zeros sub-image ... if one
% of them fits in the sub-square i will return that one and not all zeros..
if (flag ==0)
    mitotic_cell=zeros(51,51);
end