function gelbox

version_string = 'GelBox Version 2.0';

addpath(genpath('../MATLAB_Utilities'));
addpath('box','interface','update','utilities');

gui = create_interface(version_string);
