% Script to extract anterior and posterior coordinates for each data set
close all
clear 
%------------------------Set Path Specs, ID Vars------------------------%
FolderPath = 'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\';
project = 'eve7stripes_inf_2018_04_28'; %Project Identifier
% folders
fig_path = ['../../fig/experimental_system/' project '/preprocessing/'];
out_path = ['../../dat/' project '/'];
% cleaning params
keyword = '_30uW_550V'; % Keyword to ensure only sets from current project are pulled

% store set names
dirinfo = dir(FolderPath);
dirinfo(~[dirinfo.isdir]) = []; %remove non-directories
ap_filenames = {}; % ap info
set_nums = [];
for d = 1 : length(dirinfo)
    thisdir = dirinfo(d).name;
    % Skip files lacking project keyword 
    if isempty(strfind(thisdir,keyword)) 
        continue
    end    
    % append file paths    
    ap_filenames = [ap_filenames {[thisdir '/APDetection.mat']}];        
end

coord_set = NaN(length(ap_filenames),6);
header = {'set_id', 'a_x','a_y', 'p_x', 'p_y','ap_angle'};
for i = 1:length(ap_filenames) % Loop through filenames    
    % read in raw files
    load([FolderPath ap_filenames{i}]) % AP Info   
    coord_set(i,1) = i;
    coord_set(i,2:3) = coordAZoom;
    coord_set(i,4:5) = coordPZoom;
    coord_set(i,6) = atan2((coordPZoom(2)-coordAZoom(2)),(coordPZoom(1)-coordAZoom(1)));        
end

csvwrite_with_headers([out_path '\AP_XY_coordinates.csv'], coord_set, header,9); 