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

function [ volume ] = vectorToMask( vec, mask )
%VECTOR2MASK Reshapes the data of a vector to a volume. Only
%data inside the specified mask is converted.
%
% [ volume ] = vectorToMask( vec, mask )
%
% Parameter:
% ----------
% vec    ... 1-dimensional column vector containing data under mask
% mask   ... Logical Mask. Number of foreground elements must be the same
%            as elements in the vector
%            <= 0 is assumed as background and > 0 is foreground
%
% Return:
% ------
% volume ... Reshaped volume, which has the same size as mask. Vector data
%            is copied into the foreground region, background region is set
%            to zero
% 
% See also:
% MASKTOVECTOR

%Columnize vector
iptcheckinput(vec,{'logical','numeric'},...
                    {'vector','nonempty'},...
                    mfilename, 'VEC', 1);

if (ndims(vec) > 0)

volume = zeros (size(mask));

volume (mask > 0) = vec;

end

