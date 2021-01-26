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

function [ s_left, s_right, splitting_point, boundingBox ] = segmentation_left_right( s )
%SEGMENTATION_LEFT_RIGHT Returns a segmentation of the breast into left and
%right breast. RAI is assumed as image coordinates order
%
% Input:
% ------
% s      ... matrix of segmentation (0 = background, 1 = foreground)
% 
% Return:
% -------
% s_left, s_right ... Segmentation of left and right breast
% splitting point ... x coordinate of splitting point
% boundingBox     ... Bounding box of segmentation
% Author: Wolf-Dieter Vogl
% August 2012

sz = size(s);
%Find approx. y startpoint of breast by analyzing segmentation slizes
%along y axis. 
p = regionprops(s,'basic');
[~,idx] = max([p.Area]);        %When there are more than one blop (e.g. noise), use the one with the largest area
p = p(idx);
xmin = round(p.BoundingBox(2)); %note: bw methods exchange x with y
xmax = xmin + p.BoundingBox(5) - 1;

ymin = round(p.BoundingBox(1));
ymax = ymin + p.BoundingBox(4) - 1;

zmin = round(p.BoundingBox(3));
zmax = zmin + p.BoundingBox(6) - 1;
boundingBox = p.BoundingBox;

%Calculate Distance transformation to segmentation
d = mat2gray(bwdist(s(xmin:xmax,ymin:ymax,zmin:zmax),'chessboard'));
zmean = round((zmax - zmin + 1)/2);
%Distance transformation has a 'w' shape along y axis. The center of the
%'w' shape is the center between the breasts. 

%Calculate 2nd derivatives to find "center" of w shape
%Find minimum of 2nd derivative in the first 10 slices along y axis
%Location of minimum of 2nd derivative is the x position splitting the two
%breasts.

%g = diff(d(:,1:10,zmean),2,1); 
%Use central difference method
g = (d(3:end,1:10,zmean)-d(1:end-2,1:10,zmean))/2;
g = (g(3:end,:) - g(1:end-2,:))/2;
[~,idx] = min (g(:));
[x,~] = ind2sub(size(g),idx);

x = x + xmin ;

s_right = zeros(sz);
s_left = zeros(sz);
s_right(1:x,:,:) = s(1:x,:,:);s_left(x+1:end,:,:) = s(x+1:end,:,:);    %reversed x direction is used for left-right, because images are assumed to be stored in RAI (right-left, anterior...);
splitting_point = x;
end


