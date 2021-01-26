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

function pipeline(water_dicom, fat_dicom, working_folder, output_folder, template_folder)
% PIPELINE Processing pipeline to determine breast density from MRI dixon sequence
% containing a water weighted and a fat weighted volumetric breast image
% Pipeline consists of three steps:
% 1. Pre-processing: Converting dicom to nifti, reorient images, do
% bias-field correction and combine images
% 2. breast segmentation: segment breast from background and other body
% parts using a template based registration approach
% 3. segment parenchymal and fatty tissue in breast and compute summary
% measures. 
%
% It is assumed that external programs from packages "ANTs, FSL, dcm2niix, NiftyReg" are
% added to system path. 
% Template folder should contain a subfolder for each template, and each
% subfolder needs to contain a "wfdec.nii.gz" and a "s.nii.gz" file, where
% first one is the combined fat/water image and second the corresponding
% template breast segmentation. 
%
%
% Parameters
% ----------
% water_dicom : string
%    Path to dicom directory containing water weighted MRI volume
% fat_dicom : string
%    Path to dicom directory containing fat weighted MRI volume
% working_folder : string
%    Path to folder where intermediate results are saved
% output_folder : string
%    Path to folder where final results are saved (breast segmentation,
%    parenchymal/fatty tissue and density measures)
% template_folder : string
%    Path where templates for breast segmentation are located


%mkdir working and output folders
mkdir(working_folder)
mkdir(output_folder)

%-------
%Pre-processing
%-------
pre_processing (water_dicom, fat_dicom, working_folder)

%---------
%Segment breast
%---------
template_segment_breasts (working_folder, template_folder)

%---------
%Segment fatty/parenchymal tissue
%---------
measure_density (working_folder, output_folder)
