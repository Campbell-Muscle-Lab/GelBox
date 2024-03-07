function GelBox_old

version_string = 'GelBox Version 2.0';

folders = {'box','interface','update','utilities'};

for i = 1 : numel(folders)
addpath(genpath(folders{i}));
end

gui = create_interface(version_string);
