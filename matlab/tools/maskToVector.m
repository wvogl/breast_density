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

function [ vec ] = maskToVector( volume, mask )
%MASKTOVECTOR Convert the data of a volume into a vector. Only
%data inside the specified mask is converted.
%
% [ vec ] = maskToVector( volume, mask )
%
% Parameter:
% ----------
%
% volume ... Volume, which should be vectorized
% mask   ... Logical Mask. Must have same size and dimension than volume.
%            <= 0 is assumed as background and > 0 is foreground
%
% Return:
% ------
% vec    ... 1-dimensional column vector containing data under mask
%
% See also:
% VECTORTOMASK

if (size(mask(:)) ~= size(volume(:)))
    error ('Size of mask and volume must be the same.');
end
if (ndims(mask) ~= ndims(volume))
    error ('Dimension of mask and volume must be the same.');
end 

vec = volume (mask > 0);

end

