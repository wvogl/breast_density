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
%    along with breast_segmentation.  If not, see <http://www.gnu.org/licenses/>.


function  template_segment_breasts( working_folder, template_folder )
% TEMPLATE_SEGMENT_BREASTS Breasts segmentation: segment breast from background and other body
% parts using a template based registration approach.
% In a first step the best matching template breast is determined by
% affinely registering it to subject breast and then check overlap by dice
% score. The selected template is then registered non-rigid using NiftyReg
% to refine contours of template breast to subject breast. 
%
% Parameters
% ----------
% working_folder : string
%    Path to folder where intermediate results are saved
% template_folder : string
%    Path where templates for breast segmentation are located
%
 
    subject_wfdec = fullfile(working_folder,'subject_wfdec.nii.gz');
    
    %Filenames
    template_affine_wfdec = fullfile(working_folder,'template_affine_wfdec.nii.gz'); 
    affine_rigid_reg = fullfile(working_folder,'affine_rigid_reg.nii.gz');
    subject_segmented = fullfile(working_folder,'s.nii.gz');
    template_affine_seg = fullfile(working_folder,'s_affine_template.nii.gz');
    regwfdec = fullfile(working_folder,'regwfdec');
    regwfdecmat = fullfile(working_folder,'regwfdec0GenericAffine.mat');
    
    disp('Determining best matching template')
    
    subject=load_untouch_nii_gzip(subject_wfdec);

    templateFolders = getSubFoldersOfFolder(template_folder);
    numberOfTemplates = length(templateFolders);
    dice = zeros(numberOfTemplates,1);

    %Iterate over all templates and register each one affinely to subject
    for i = 1:numberOfTemplates
        %Path to template
        template_wfdec = fullfile(templateFolders{i},'wfdec.nii.gz');

        disp (['Registering template ' templateFolders{i}])
        [status output] = system(['antsRegistration -d 3 -o "' regwfdec '" -c 10000x10000x1 -s 2x1x0 -f 4x2x1 -m MI["' subject_wfdec '","' template_wfdec '",1,32] -t Affine[0.1]']);

        disp(['ApplyTransformation ' num2str(i)]);
        [status output] = system(['antsApplyTransforms -d 3 -i "' template_wfdec '" -o "' template_affine_wfdec '" -r "' subject_wfdec '" -t "' regwfdecmat '"']);

        
        template = load_untouch_nii_gzip(template_affine_wfdec);
        dice(i) = LabelOverlapMeasures(subject.img>0.1,template.img>0.1);
        disp(['Dice: ' num2str(dice(i))]);
    end

    [maxDice, iMax]=max(dice);
    
    disp(['Template ' templateFolders{iMax} ' chosen. DICE: ' num2str(maxDice)])
    template_wfdec = fullfile(templateFolders{iMax},'wfdec.nii.gz');
    template_seg = fullfile(templateFolders{iMax},'s.nii.gz');
    
    %Do final registration
    %Two steps, first affine registration with ants and then non-rigid
    %using nifty-reg as it is much faster when using GPU
   	disp ('Final registration')
    disp ('Affine registration')
    [status output] = system(['antsRegistration -d 3 -o "' regwfdec '" -c 10000x10000x1 -s 2x1x0 -f 4x2x1 -m MI["' subject_wfdec '","' template_wfdec '",1,32] -t Affine[0.1]']);
    
    %Warp combined w/f image and template segmentation. 
	[status output] = system(['antsApplyTransforms -d 3 -i "' template_wfdec '" -o "' template_affine_wfdec '" -r "' subject_wfdec '" -t "' regwfdecmat '"']);
    [status output] = system(['antsApplyTransforms -d 3 -i "' template_seg '" -o "' template_affine_seg '" -r "' subject_wfdec '" -t "' regwfdecmat '"']);
    
    %Non-rigid registration of combined affine registered w/f image
    disp('Non-rigid registration');    
	[status output] = system(['reg_f3d -flo "' template_affine_wfdec '" -ref "' subject_wfdec '" -gpu -maxit 10000 -cpp "' affine_rigid_reg '"']); %remove -gpu parameter when it is not supported by reg_f3d
	
	disp('Warp template segmentation');
	[status output] = system(['reg_resample -flo "' template_affine_seg '" -ref "' subject_wfdec '" -cpp "' affine_rigid_reg '" -res "' subject_segmented '" -NN']);

end

