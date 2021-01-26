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

function [ subfolder ] = getSubFoldersOfFolder( folder, wildcard, sortnat, addfolder  )
%GETFILESOFFOLDER Returns a cell array with all filenames of a specific
%folder. Subdirectories and '.','..' are not in the filelist.
%
%Parameter
%---------
%folder   ... folder from where subfolders should be listed
%wildcard ... wildcard to filter the folderlist (eg '*.zip') (default: all folders)
%sortnat  ... flag, do a natural order sort of the filelist (eg. 'a1.zip','a2.zip','a10.zip') (default: 0)
%addfolder ... flag, add parent folder to filename (default: 1)
if (nargin < 4)
    addfolder = 1;
end

if (nargin < 3)
    sortnat = 0;
end
if (nargin < 2)
    wildcard = '*';
end

mask_files = dir(fullfile(folder,wildcard));
%Remove . and ..
mask_files(structfind (mask_files,'name','.'))=[];
mask_files(structfind (mask_files,'name','..'))=[];

isdir = logical(cat(1,mask_files.isdir));

mask_files = mask_files(isdir);

if (addfolder ~= 0)
    subfolder = cellfun(@(x)fullfile(folder,x),{mask_files.name},'UniformOutput',false);
else
    subfolder = cellfun(@(x)x,{mask_files.name},'UniformOutput',false);
end

if (sortnat ~= 0)
    subfolder = sort_nat(subfolder);
end

end

