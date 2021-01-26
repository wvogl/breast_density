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

function [MO, UO, TO, VS, FP, FN, PR] = LabelOverlapMeasures(gt, seg)
%LABELOVERLAPMEASURES This function provides several measures for the
%similarity of a ground truth labeling and a segmentation labeling.
%
%Parameter:
%----------
% gt  ... Ground Truth (logical or 0/1 array)
% seg ... Segmentation (must have same size as ground truth)
%
% Return:
%--------
% Note: G = ground truth, S = segmentation label
% MO ...   Mean Overlap (Dice Similarity Coefficient) MO = 2*|G n S|/(|G| + |S|)
% UO ...   Union Overlap (Jaccard Similarity Coefficient) JO = |G n S|/|G u S|
% TO ...   Target Overlap (Sensitivity, Recall) TO = |G n S| / |G|
% VS ...   Volume Similarity VS = 2 * (|S| - |G|) / (|S| + |G|)
% FP ...   False Positive Error FP = |S\G| / |S|
% FN ...   False Negative Error FN = |G\S| / |G|
% PR ...   Precision PR = |G n S| / |S|
% Detailed description of overlap measures given in:
% Klein et al. "Evaluation of 14 nonlinear deformation algorithms applied 
% to human brain MRI registration" (2009)
%
%Author: Wolf-Dieter Vogl
%Todo: Support more than one label

%Ensure that gt and seg are vectorized
%
gt = logical(gt(:));
seg = logical(seg(:));

%Support variables
SUM_GT = sum(gt == 1);
SUM_SEG = sum(seg == 1);
UNION = sum(gt == 1 | seg==1);
INTERSECTION = sum(gt == 1 & seg == 1);
MEANVOLUME = SUM_GT + SUM_SEG;

%Overlap measures
MO = 2 * INTERSECTION / MEANVOLUME;
UO = INTERSECTION / UNION;
TO = INTERSECTION / SUM_GT;
VS = 2 * (SUM_SEG - SUM_GT) / MEANVOLUME;
FP = sum(seg == 1 & gt ~= 1) / SUM_SEG;
if (isnan(FP))
    FP = 1;
end
FN = sum(gt == 1 & seg ~= 1) / SUM_GT;
PR = INTERSECTION / SUM_SEG;
end