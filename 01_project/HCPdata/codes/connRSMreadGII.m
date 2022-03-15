%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code by Selen Atasoy
% 05/12/2013
% 
% reads the .gii file of the Human Connectome Project
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vertices, faces] = connRSMreadGII(file_path, subject_name, surf_type)

if strcmp(subject_name,'')
    file_name.left  = ['L.', surf_type, '.32k_fs_LR.surf.gii'];
    file_name.right = ['R.', surf_type, '.32k_fs_LR.surf.gii'];
else
    file_name.left  = [subject_name, '.L.', surf_type, '.32k_fs_LR.surf.gii'];
    file_name.right = [subject_name, '.R.', surf_type, '.32k_fs_LR.surf.gii'];
end

temp = gifti(fullfile(file_path, file_name.left));
vertices.left   = temp.vertices;
faces.left      = temp.faces;
clear temp;

temp = gifti(fullfile(file_path, file_name.right));
vertices.right   = temp.vertices;
faces.right      = temp.faces;

vertices.all    = [vertices.left; vertices.right];
faces.all       = [faces.left; faces.right+size(vertices.left,1)];

