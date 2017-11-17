% ----------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and           %
% simulation. See http://opensim.stanford.edu and the NOTICE file         %
% for more information. OpenSim is developed at Stanford University       %
% and supported by the US National Institutes of Health (U54 GM072970,    %
% R24 HD065690) and by DARPA through the Warrior Web program.             %
%                                                                         %   
% Copyright (c) 2005-2012 Stanford University and the Authors             %
% Author(s): Dan Lichtwark                                                %
%                                                                         %
% Licensed under the Apache License, Version 2.0 (the "License");         %
% you may not use this file except in compliance with the License.        %
% You may obtain a copy of the License at                                 %
% http://www.apache.org/licenses/LICENSE-2.0.                             %
%                                                                         % 
% Unless required by applicable law or agreed to in writing, software     %
% distributed under the License is distributed on an "AS IS" BASIS,       %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         %
% implied. See the License for the specific language governing            %
% permissions and limitations under the License.                          %
% ----------------------------------------------------------------------- %

% strengthScaler.m                                                        
% Author: Dan Lichtwark

function strengthScaler_Neck(scaleFactor1,scaleFactor2, Model_In, Model_Out)
% OSIMstrength_scaler(scaleFactor, Model_In, Model_Out)
% Test program to load muscles and change strength of muscles and re-save
% model
%
% Inputs - scaleFactor (double) - amount to scale all muscle forces
%          Model_In (string) - existing model path and file name 
%          Model_Out (string) - new model path and file name 
%
% eg. strengthScaler(2)
% eg. strengthScaler(2, 'mySimpleBlockModel.osim')
% eg. strengthScaler(2, 'mySimpleBlockModel.osim', 'myStrongerBlockModel.osim')
%
% Author: Glen Lichtwark (The University of Queensland)
% with invaluable assistance from Ayman Habib (Stanford University)
% Initial code replicating the muscleStrengthScaler.cpp file developed by
% Edith Arnold and Ajay Seth

import org.opensim.modeling.*

error(nargchk(1, 4, nargin));

if nargin < 3
    [Model_In, path] = uigetfile('.osim');
    fileoutpath = [Model_In(1:end-5),'_MuscleScaled.osim'];    
    filepath = [path Model_In];
elseif nargin < 4
    fileoutpath = [Model_In(1:end-5),'_MuscleScaled.osim'];
    filepath = Model_In;
else
    filepath = Model_In;
    fileoutpath = Model_Out;
end

%Create the Original OpenSim model from a .osim file
Model1 = Model(filepath);
Model1.initSystem;

% Create a copy of the original OpenSim model for the Modified Model
Model2 = Model(Model1);
Model2.initSystem;

% Rename the modified Model so that it comes up with a different name in
% the GUI navigator
Model2.setName('modelModified');

% Get the set of muscles that are in the original model
Muscles1 = Model1.getMuscles(); 

%Count the muscles
nMuscles = Muscles1.getSize();

disp(['Number of muscles in orginal model: ' num2str(nMuscles)]);

% Get the set of forces that are in the scaled model
% (Should be the same as the original at this point.)
Muscles2 = Model2.getMuscles();

flexion_group={ 'cleid_mast' 'cleid_mast_l' 'cleid_occ' 'cleid_occ_l' 'long_cap_sklc4' 'long_cap_sklc4_l' 'long_col_c1c5'...
    'long_col_c1c5_l' 'long_col_c1thx' 'long_col_c1thx_l' 'long_col_c5thx' 'long_col_c5thx_l' 'scalenus_ant'...
    'scalenus_ant_l' 'stern_mast' 'stern_mast_l'};

extension_group={ 'deepmult-C4/5-C2' 'deepmult-C4/5-C2_l' 'deepmult-C5/6-C3' 'deepmult-C5/6-C3_l' 'deepmult-C6/7-C4' 'deepmult-C6/7-C4_l' ...
    'deepmult-T1-C5' 'deepmult-T1-C5_l' 'deepmult-T1-C6' 'deepmult-T1-C6_l' 'deepmult-T2-C7' 'deepmult-T2-C7_l' 'deepmult-T2-T1' ...
    'deepmult-T2-T1_l' 'iliocost_cerv_c5rib' 'iliocost_cerv_c5rib_l' 'levator_scap' 'levator_scap_l' 'longissi_cap_sklc6' ...
    'longissi_cap_sklc6_l' 'longissi_cerv_c4thx' 'longissi_cerv_c4thx_l' 'obl_cap_inf' 'obl_cap_inf_l' 'obl_cap_sup' 'obl_cap_sup_l' ...
    'rectcap_post_maj' 'rectcap_post_maj_l' 'rectcap_post_min' 'rectcap_post_min_l' 'scalenus_med' 'scalenus_med_l' 'scalenus_post' ...
    'scalenus_post_l' 'semi_cap_sklc5' 'semi_cap_sklc5_l' 'semi_cap_sklthx' 'semi_cap_sklthx_l' 'semi_cerv_c3thx' 'semi_cerv_c3thx_l' ...
    'splen_cap_sklc6' 'splen_cap_sklc6_l' 'splen_cap_sklthx' 'splen_cap_sklthx_l' 'splen_cerv_c3thx' 'splen_cerv_c3thx_l' 'supmult-C4/5-C2' ...
    'supmult-C4/5-C2_l' 'supmult-C5/6-C2' 'supmult-C5/6-C2_l' 'supmult-C6/7-C2' 'supmult-C6/7-C2_l' 'supmult-T1-C4' 'supmult-T1-C4_l' 'supmult-T1-C5'...
    'supmult-T1-C5_l' 'supmult-T2-C6' 'supmult-T2-C6_l' 'trap_acr trap_acr_l' 'trap_cl trap_cl_l'};

% loop through forces and scale muscle Fmax accordingly (index starts at 0)
for i = 0:nMuscles-1
        
    %get the muscle that the original muscle set points to
    %to read the muscle type and the max isometric force
    currentMuscle = Muscles1.get(i);
    
    %define the muscle in the modified model for changing
    newMuscle = Muscles2.get(i);

    %define the new muscle force by multiplying current muscle max
    %force by the scale factor
    
    % Separated scale factors in relation to the muscle group
    for n=1:length(extension_group),
    if strcmp(newMuscle,extension_group(n))==1
    
       newMuscle.setMaxIsometricForce(currentMuscle.getMaxIsometricForce()*scaleFactor1);
    break
    else
        newMuscle.setMaxIsometricForce(currentMuscle.getMaxIsometricForce()*scaleFactor2);
   
    end    
    end
    

end
 
% save the updated model to an OSIM xml file
Model2.print(fileoutpath)
disp(['The new model has been saved at ' fileoutpath]);

end
