%    This file is part of breast_density.

%    breast_density is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    breast_density is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with breast_density.  If not, see <http://www.gnu.org/licenses/>.

function [s_left,s_right,split_point, boundingBox] = split_left_right(s)
%Split the segmentation in a left and right part, assuming that there are
%at least two segmentation parts
%Parameter:
%s      ... Segmentation where both breasts are seperate blobs
% 
%Return:
%s_left      ... Segmentation of left breast
%s_right     ... Segmentation of right breast
%split_point ... Point halfway inbetween left and right segmentation
%bounding box
%boundingBox ... Overall bounding box of segmentation (s_full)

%Determine overall bounding box
%When there are more than one blop (e.g. noise), use the one with the largest area
p = regionprops(s,'basic');
[~,idx] = max([p.Area]);       
boundingBox = p(idx).BoundingBox;

%Determine two breast blobs
cc = bwconncomp(s);
p = regionprops(cc,'basic');

%Sort blobs by size
[~,idx] = sort([p.Area],'Descend');
%lm = labelmatrix(cc);

if (numel(idx) > 1 )
    if (p(idx(1)).Area >= 1000 && p(idx(2)).Area >= 1000)
        xmin(1) = round(p(idx(1)).BoundingBox(2)); %note: bw methods exchange x with y
        xmax(1) = xmin(1) + p(idx(1)).BoundingBox(5) - 1;

        xmin(2) = round(p(idx(2)).BoundingBox(2)); %note: bw methods exchange x with y
        xmax(2) = xmin(2) + p(idx(2)).BoundingBox(5) - 1;

        %Split point is the middle distance between bounding box borders
        if (xmin(2) > xmin(1))
            split_point = round((xmin(2) - xmax(1))/2) + xmax(1);
        else
            split_point = round((xmin(1) - xmax(2))/2) + xmax(2);
        end
        sz = size(s);
        s_right = zeros(sz);
        s_left = zeros(sz);
        s_right(1:split_point,:,:) = s(1:split_point,:,:);s_left(split_point+1:end,:,:) = s(split_point+1:end,:,:);
    else
        warning('Only one blob found with size >= 1000');
        s_left = 0;
        s_right = 0;
        split_point = 0;
    end
else
    warning('Only one blob found');
    s_left = 0;
    s_right = 0;
    split_point = 0;
end
end

