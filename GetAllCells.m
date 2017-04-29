%-------------------------------------------------------------------------%
%  Project       : Mitosis-Detection-Breast-Cancer                        %
%  File          : GetAllCells.m                                          %
%  Description   : Export the cell sub-images and Coordinates             %
%  Author        : Monachopoulos Konstantinos                              %
%-------------------------------------------------------------------------%

function [test_cells test_cells_coor]=GetAllCells(UnnormalizedImage)

% Normalize every Histological Image to make it independent of its colorant value ...
[Inorm1 H1 E1] = normalizeStaining(UnnormalizedImage);

% Threshold the normalized image ...
BW=thresh_by_iterations(Inorm1);

% Compare and hold only the set of pixels that represents the pattern 'octagon' (as it is close to the cell shape)
BWerode = imerode(BW, strel('octagon',3));

% Clear the border of the image at the edges..
BWfinal = imclearborder(BWerode, 6);

% Label every area
[L, n] = bwlabel(BWfinal);

% Call regionprops extracting 'Area','Eccentricity','Extent','Orientation'
Region = regionprops(L,'Area','Eccentricity','Extent','Orientation');

% keep only the areas of gathered pixels that are more of 100 and less that 1500
% keep only the areas of gathered pixels that Eccentricity is less than 1 but cut these that look straight lines

idx=([Region.Area] > 100)&([Region.Area] < 1500)&([Region.Extent]<0.95)...
    &([Region.Eccentricity]<0.97);
BW2 = ismember(L,find(idx));

% Label every remaining area
[S, m] = bwlabel(BW2);
region2 = regionprops(S,'BoundingBox','centroid');
BW2 = ismember(S,find(idx));

% Square-plot every cell so the center of the cell will be in the midlle of the square, plot the results...
figure;
imshow(Inorm1);
for i=1:size(region2,1)
    hold on;
    region2(i).BoundingBox(:,1)=region2(i).BoundingBox(:,1)-...
        ((50-region2(i).BoundingBox(:,3))/2);
    region2(i).BoundingBox(:,2)=region2(i).BoundingBox(:,2)-...
        ((50-region2(i).BoundingBox(:,4))/2);
    region2(i).BoundingBox(:,3:4)=[50 50];
    rectangle('Position', [region2(i).BoundingBox],'edgecolor','r');
end
hold off;

% Cut every square-cell area image and extract it with the analogous coordinates .. Going for LBP ...

for i=1:size(region2,1)
    if (size(rgb2gray(imcrop(Inorm1, region2(i).BoundingBox))) == [51 51])
        test_cells(:,:,i) = rgb2gray(imcrop(Inorm1, region2(i).BoundingBox));
        test_cells_coor(i,1:2)=region2(i).Centroid;
    end
end
