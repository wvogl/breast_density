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

function [density] = measure_density( working_folder, output_folder)
%MEASURE_DENSITY Measure the density and volume of the breast. In particular
%percentage of fat and water is measured. For that reason breasts are
%segmented into parenchymal and fatty tissue. 
%
%Return:
%-------
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
%Notes:
%------
%The working_folder needs to have
%following files included:
%f.nii.gz Fat-weighted dixon image
%w.nii.gz Water-weighted dixon image
%s.nii.gz Binary Segmentation file of breast
%subject_wfdec.nii.gz combined f.nii.gz and w.nii.gz image using "ImageMath Decision" 
%All images needs to be oriented in RL AP IS 

%Orientation can be changed with fsl tool:
%fslswapdim filein.nii.gz RL AP IS fileout.nii.gz
%
%In the output folder following information images are stored:
%hist.png      2D-Histogram of fat-water distribution in images
%hist_div.png  Illustration of the division of the histogram
%div_left_right.png  Central axial slice of breast segmentation and division position into left/right breast. 
%assignment.png Central axial slice of breast segmentation into parenchymal/fatty tissue. 
%density.mat   Density measurement structure (as listed above). 

%Author:
%-------
%Wolf-Dieter Vogl

%Generate results folder

% Loading and preprocessing
% -------------------------

%Load water and fat image
w = load_untouch_nii_gzip(fullfile(working_folder, 'w.nii.gz'));
f = load_untouch_nii_gzip(fullfile(working_folder, 'f.nii.gz'));

%Load combined water/fat image 
wfdec = load_untouch_nii_gzip(fullfile(working_folder, 'subject_wfdec.nii.gz'));

%Load normalized fat + water image (don't load anymore, computed below)
%fw = load_untouch_nii_gzip(fullfile(patfolder,[prefix 'fw.nii.gz']));

%Convert to 16 bit integer, no matter what input format is like, 
%cut negative values. 
w.img = int16(round(w.img));
f.img = int16(round(f.img));
w.img(w.img < 0) = 0;
f.img(f.img < 0) = 0;

w.hdr.dime.datatype = 4;
f.hdr.dime.datatype = 4;
w.hdr.dime.bitpix = 16;
f.hdr.dime.bitpix = 16;

%Compute fat+water image 
fw = w;
fw.img = double(f.img + w.img);
%Normalize to [0,1]
fw.img = (fw.img - min(fw.img(:)))/(max(fw.img(:))-min(fw.img(:)));

%Load segmentation
s = load_untouch_nii_gzip(fullfile(working_folder,'s.nii.gz'));

dimension = w.hdr.dime.pixdim(2:4);

assignnifti = w;            %store nifty header for assignment image
assignnifti.img = [];

w = w.img;
f = f.img;

%s_without_middle only contains segmentations with label 1 (0 < x < 1.5),skipping other labels. 
%s contains both labels as binary mask. 
s = (s.img > 0 & s.img < 1.5 & fw.img > 0.05 & wfdec.img > 0.01); %wfdec is used to remove segmentation voxels outside of volume
s = imfill(s,'holes');
clear wfdec;
   
density = struct('volume',0,'volume_left',0,'volume_right',0,'fat',0,...
                 'fat_left',0,'fat_right',0,'water',0,'water_left',0,...
             'water_right',0,'split_point',0);

%Split into left/right breast
%----------------------------
         
%Split between left and right breast, asuming that the segmentations of
%these two breasts are not connected. 
[s_left,s_right,density.split_point, boundingBox] = split_left_right(s);
%Check if splitting worked
if (sum(s_left(:)) < 1000 || sum(s_right(:)) < 1000)
    %Splitting did not work, breasts are probably connected, use another algorithm based on distance transform. 
    [s_left,s_right,density.split_point, boundingBox] = segmentation_left_right(s);
end

%Bounding box of segmentation
split_point=density.split_point;
xmin = round(boundingBox(2)); %note: bw functions exchange x with y
xmax = xmin + boundingBox(5) - 1;

ymin = round(boundingBox(1));
ymax = ymin + boundingBox(4) - 1;

zmin = round(boundingBox(3));
zmax = zmin + boundingBox(6) - 1;
zmean = round((zmax - zmin + 1)/2);

%Calculate density
%-----------------
[density,assignnifti.img] = cluster_density_histogram(f,w,s,s_left,s_right,dimension,output_folder);

density.split_point=split_point;

disp(['Left Water Percentage: ',num2str(density.water_left),'%']);
disp(['Left Fat Percentage: ',num2str(density.fat_left),'%']);
disp(['Left Volume: ',num2str(density.volume_left)]);

disp(['Right Water Percentage: ',num2str(density.water_right),'%']);
disp(['Right Fat Percentage: ',num2str(density.fat_right),'%']);
disp(['Right Volume: ',num2str(density.volume_right)]);

disp(['Water Percentage: ',num2str(density.water),'%']);
disp(['Fat Percentage: ',num2str(density.fat),'%']);
disp(['Volume: ',num2str(density.volume)]);

%Save assignment nifti
save_untouch_nii_gzip(assignnifti,fullfile(output_folder, 'assignment.nii.gz'));
%clear assignnifti;

%Plot and save segmentations, save results
%---------------------

%Plot split point and breast segmentation of central axial slice
w_axi = mat2gray(w(:,:,zmin+zmean)');
w_axi_color(:,:,1) = w_axi;
w_axi_color(:,:,2) = w_axi;
w_axi_color(:,:,3) = double(s(:,:,zmin+zmean)');

h = figure;
imagesc(w_axi_color);
hold on;
plot ([density.split_point density.split_point],[1 size(w,2)],'r-','LineWidth',2);
saveas(h,fullfile(output_folder,'div_left_right'),'png');
close(h);

%Plot parenchymal/fatty tissue segmentation of central axial slice
h = figure;
w_axi = mat2gray(w(:,:,zmin+zmean)');
w_axi_color(:,:,1) = w_axi;
w_axi_color(:,:,2) = double(assignnifti.img(:,:,zmin+zmean)'==1);
w_axi_color(:,:,3) = double(assignnifti.img(:,:,zmin+zmean)'==2);

imagesc(w_axi_color);
hold on;
saveas(h,fullfile(output_folder,'assignment'),'png');
close(h);

save (fullfile(output_folder,'density.mat'),'density');


end


