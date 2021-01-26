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
%    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

function [ files ] = getFilesOfFolder( folder, wildcard, sortnat, addfolder )
%GETFILESOFFOLDER Returns a cell array with all filenames of a specific
%folder. Subdirectories and '.','..' are not in the filelist.
%
%Parameter
%---------
%folder   ... folder from where files should be listed
%wildcard ... wildcard to filter the filelist (eg '*.zip') (default: all files)
%sortnat  ... flag, do a natural order sort of the filelist (eg. 'a1.zip','a2.zip','a10.zip') (default: 0)
%addfolder ... flag, add folder to filename (default: 1)
%
% Author: Wolf-Dieter Vogl
% Mail: wolf-dieter.vogl at meduniwien.ac.at
% Date: March 2013

if (nargin < 4 || isempty(addfolder))
    addfolder = 1;
end

if (nargin < 3 || isempty(sortnat))
    sortnat = 0;
end
if (nargin < 2 || isempty(wildcard))
    wildcard = '*';
end

mask_files = dir(fullfile(folder,wildcard));

isfile = ~logical(cat(1,mask_files.isdir));

mask_files = mask_files(isfile);

if (addfolder ~= 0)
    files = cellfun(@(x)fullfile(folder,x),{mask_files.name},'UniformOutput',false);
else
    files = cellfun(@(x)x,{mask_files.name},'UniformOutput',false);
end

if (sortnat ~= 0)
    files = sort_nat(files);
end

end

