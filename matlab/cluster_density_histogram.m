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

function [ density, assignment ] = cluster_density_histogram( f,w,s, s_left,s_right,pixdim,folder)
%CLUSTER_DENSITY_HISTOGRAM Calculate the 2d-histogram of fat and water
%image, determine the best dividing point, and calculate percentage of fat
%and water according to dividing point.
%
%Input
%-----
%f,w,s   Matrix of fat, water and binary segmentation image
%pixdim Dimension of a voxel in mm
%folder  (optional) folder where displayed diagrams are stored as .png
%
%Output
%------
%density structure
%   density.name             Name of the folder
%   density.volume           Volume of breasts (in cm³)
%   density.volume_left      Volume of left breast (in cm³)
%   density.volume_right     Volume of right breast (in cm³)
%   density.fat              Percentage of Fat in breasts
%   density.fat_left/right   Percentage of fat in left/right breast
%   density.water            Percentage of Water in breasts
%   density.water_left/right Percentage of Water in left/right breast
%   density.split_point      Splitting point (x coords) of left and right
%                            breast
%
%assignment ... Volume containing the assigned label for each voxel (1 = fat, 2 = water, 0 = background)

%Constants
BIN_SIZE=15;


density = struct('volume',0,'volume_left',0,'volume_right',0,'fat',0,...
                 'fat_left',0,'fat_right',0,'water',0,'water_left',0,...
             'water_right',0,'split_point',0);


if (nargin < 7)
    store_pix = false;
else
    store_pix = true;
end

fs=f(s>0);
ws=w(s>0);
X = double([fs ws]);
maxf = max(X(:,1)) + 1;
maxw = max(X(:,2)) + 1;
clear fs;
clear ws;

%Calculate histogram and binned histogram 
A = accumarray(X+1,1,[maxf maxw]);
N=hist3(double(X),[BIN_SIZE BIN_SIZE]);

%Left Breast
fs=f(s_left>0);
ws=w(s_left>0);
X_left = double([fs ws]);
clear fs;
clear ws;
A_left = accumarray(X_left+1,1,[maxf maxw]);

%Right Breast
fs=f(s_right>0);
ws=w(s_right>0);
X_right = double([fs ws]);
clear fs;
clear ws;
A_right = accumarray(X_right+1,1,[maxf maxw]);


%Find extrema in binned histogram
J = imextendedmax(N,500);

[L, numL] = bwlabel (J);
disp ('Morphological operations');
disp (['Number of peaks detected: ' num2str(max(L(:)))]);
stats = regionprops(L,'Centroid');
 
if (numL > 2)
    warning ('More than two peaks\nPeaks with highest histogram values chosen');
    value = [];
    for i = 1:numL        
        if (max(stats(i).Centroid) == 1)
            %Skip values close to (0,0)
            value = [value -inf];
        else
            r = round(stats(i).Centroid);
            value = [value N(r(2),r(1))];
        end
    end
    [~,idx] = sort(value,'descend');
    %Chose peaks with highest value, and peaks have to be on opposite site
    %If no further peak on opposite site is found, first peak is mirrored
    bincenter1 = stats(idx(1)).Centroid;
    bincenter2 = [bincenter1(2) bincenter1(1)];
    for i = 2:numL
        if (value(idx(i)) == -inf)
            continue;
        end
        t = stats(idx(i)).Centroid;
        if (bincenter1(1) >= bincenter1(2) && t(1) <= t(2) || ...
            bincenter1(1) <= bincenter1(2) && t(1) >= t(2))
            bincenter2 = t;
            break;
        end
    end
    
elseif (numL == 2)
    bincenter1 = stats(1).Centroid;
    bincenter2 = stats(2).Centroid;
    %Check if detected peaks makes sense 
    if (~(bincenter1(1) >= bincenter1(2) && bincenter2(1) <= bincenter2(2) || ...
        bincenter1(1) <= bincenter1(2) && bincenter2(1) >= bincenter2(2)))
        warning ('Cluster centers are not on different sides, peak with highest value is chosen and mirrored');
        value = [];
        for i = 1:numL 
            if (max(stats(i).Centroid) == 1)
                %Skip values close to (0,0)
                value = [value -inf];
            else
                r = round(stats(i).Centroid);
                value = [value N(r(2),r(1))];
            end
        end
        [~,idx] = sort(value,'descend');
        bincenter1 = stats(idx(1)).Centroid; 
        bincenter2 = [bincenter1(2) bincenter1(1)];
    end
elseif (numL == 1)
    bincenter1 = stats(1).Centroid;
    bincenter2 = [bincenter1(2) bincenter1(1)];  %No second peak, mirror center1 point along 45° axis by exchanging x and y
end

%Remove peaks close to (0,0)
if (max(bincenter1) == 1)
    warning ('One peak close to (0 0), peak is dropped and mirrored from other peak');
    bincenter1 = [bincenter2(2) bincenter2(1)];
end
if (max(bincenter2) == 1)
    warning ('One peak close to (0 0), peak is dropped and mirrored from other peak');
    bincenter2 = [bincenter1(2) bincenter1(1)];
end

h = double(max(X)) ./ BIN_SIZE;  %Size of bin intervall
    %Convert from binned coordinates (1..BIN_SIZE) to original histogram coordinates (center of interval (=h/2) is used)
    %x and y are exchanged due to results of regionprops are in image coordinates
center1 = ([bincenter1(2) bincenter1(1)] - 1) .* h + h/2; 
center2 = ([bincenter2(2) bincenter2(1)] - 1) .* h + h/2; 
center = (center1+center2)./2;

%Plot histogram
h1 = figure;
imagesc(log(A'));
set(gca,'YDir','normal');
%pcolor(xb,yb,n1);

ylabel('water');
xlabel('fat');

h2 = figure;
%pcolor(xb,yb,n1);
%imagesc(log(A'));


%Calculate percentage
%a=center(1)/center(2);
a = center(2)/center(1);
s1 = a*repmat(single(1:size(A,1)),[size(A,2),1])';
s2 = repmat(single(1:size(A,2)),[size(A,1),1]);
s_split = s1-s2;
s_split = s_split>0;
clear s1 s2
fw_sum = sum(A(:));
f_sum = sum(A(:).*s_split(:));
w_sum = fw_sum-f_sum;   
density.volume = fw_sum * pixdim(1)*pixdim(2)*pixdim(3) / 1000.0;  %Divide by 1000 to convert from mm³ to cm³
density.fat = f_sum / fw_sum * 100;
density.water = w_sum / fw_sum * 100;

%Calculate left breast percentage
fw_sum = sum(A_left(:));
f_sum = sum(A_left(:).*s_split(:));
w_sum = fw_sum-f_sum;   
density.volume_left = fw_sum * pixdim(1)*pixdim(2)*pixdim(3) / 1000.0;  %Divide by 1000 to convert from mm³ to cm³
density.fat_left = f_sum / fw_sum * 100;
density.water_left = w_sum / fw_sum * 100;

%Calculate right breast percentage
fw_sum = sum(A_right(:));
f_sum = sum(A_right(:).*s_split(:));
w_sum = fw_sum-f_sum;   
density.volume_right = fw_sum * pixdim(1)*pixdim(2)*pixdim(3) / 1000.0;  %Divide by 1000 to convert from mm³ to cm³
density.fat_right = f_sum / fw_sum * 100;
density.water_right = w_sum / fw_sum * 100;


%Plot splitting
imagesc(s_split' - mat2gray(A').*(s_split'.*2-1));
set(gca,'YDir','normal');
hold all;
plot (center1(1)-1,center1(2)-1,'ko','MarkerSize',12,'LineWidth',2); 
plot (center2(1)-1,center2(2)-1,'ko','MarkerSize',12,'LineWidth',2);

plot (center(1)-1,center(2)-1,'ko','MarkerSize',12,'LineWidth',2);
plot (center(1)-1,center(2)-1,'kx','MarkerSize',12,'LineWidth',2);
ylabel('water');
xlabel('fat');

idx = sub2indFast(size(A),X(:,1)+1,X(:,2)+1);
assignment = 1+ s_split(idx) ;                  %  + 1 to get label 1/2 (fat/water)
assignment = vectorToMask(assignment,s);        % convert vector back to volume

if (store_pix)
    saveas(h1,fullfile(folder,'hist'),'fig');
    saveas(h1,fullfile(folder,'hist'),'png');
    saveas(h2,fullfile(folder,'hist_div'),'fig');
    saveas(h2,fullfile(folder,'hist_div'),'png');
    
    close(h1);
    close(h2);
end

end

