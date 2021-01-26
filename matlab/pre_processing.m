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

function pre_processing(water_dicom, fat_dicom, working_folder)
% PRE_PROCESSING Converting dicom to nifti, reorient images, do
% bias-field correction and combine images.
% Images are swaped to be in RL AP IS orientation. Fat weighted image is
% resampled such that it has same pixel spacing as water weighted image. 
%
% Following programs needs to be present in path:
% "fslswapdim (from FSL tools), antsApplyTransforms, N4BiasFieldCorrection, ImageMath (from ANTs), dcm2niix"
% Parameters
% ----------
% water_dicom : string
%    Path to dicom directory containing water weighted MRI volume
% fat_dicom : string
%    Path to dicom directory containing fat weighted MRI volume
% working_folder : string
%    Path to folder where intermediate results are saved (bias field
%    correction

%Filenames
fnii = ['"' fullfile(working_folder,'f.nii.gz') '"'];
wnii = ['"' fullfile(working_folder,'w.nii.gz') '"'];
fnii2 = [fullfile(working_folder,'f.nii.gz')];
wnii2 = [fullfile(working_folder,'w.nii.gz')];


fbnii = ['"' fullfile(working_folder,'fb.nii.gz') '"'];
wbnii = ['"' fullfile(working_folder,'wb.nii.gz') '"'];
fw = ['"' fullfile(working_folder,'fw.nii.gz') '"'];
fn = ['"' fullfile(working_folder,'fn.nii.gz') '"'];
wn = ['"' fullfile(working_folder,'wn.nii.gz') '"'];
subject_wfdec = ['"' fullfile(working_folder,'subject_wfdec.nii.gz') '"'];        
      
disp ('Convert from DICOM to NIFTI');

%If f.nii or w.nii already exists, remove them, as dcm2niix do not
%overwrite files
if (exist(fnii2,'file')==2) 
    delete(fnii2);
end
if (exist(wnii2,'file')==2)
    delete(wnii2);
end

[status output] = system (['dcm2niix -9 -f "f" -o "' working_folder '" -v 1 -z y "' fat_dicom '/"']);
[status output] = system (['dcm2niix -9 -f "w" -o "' working_folder '" -v 1 -z y "' water_dicom '/"']);

disp ('Reorient images');
[status output] = system (['fslswapdim ' wnii ' RL AP IS ' wnii]);
[status output] = system (['fslswapdim ' fnii ' RL AP IS ' fnii]);

disp ('Resample fat weighted image');
[status output] = system (['antsApplyTransforms -d 3 -n NearestNeighbor -i ' fnii ' -o ' fnii ' -r ' wnii ]);

disp('BiasField Correction');
[status output] = system (['N4BiasFieldCorrection -d 3 -i ' wnii ' -o ' wbnii ' -c [150x100,0.001] -b [20] -s 2']);
[status output] = system (['N4BiasFieldCorrection -d 3 -i ' fnii ' -o ' fbnii ' -c [150x100,0.001] -b [20] -s 2']);

disp ('Calculate f+w');
[status output] = system (['ImageMath 3 ' fw ' + ' wbnii ' ' fbnii ]);
[status output] = system (['ImageMath 3 ' fw ' Normalize ' fw ]);

disp ('Calculate wf decision (combine w and f)');
%First normalize bias-field corrected images
[status output] = system (['ImageMath 3 ' fn ' Normalize ' fbnii]);
[status output] = system (['ImageMath 3 ' wn ' Normalize ' wbnii]);

%Combine using ImageMath Decision function from ANTs
%Combined image is used to register to template
[status output] = system (['ImageMath 3 ' subject_wfdec ' Decision ' wn ' ' fn ]);
[status output] = system (['ImageMath 3 ' subject_wfdec ' Normalize ' subject_wfdec  ]);

end