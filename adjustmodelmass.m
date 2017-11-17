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

function adjustmodelmass(scaleFactor, Model_In, Model_Out,RAA_Log,SCALE_Log)
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

error(nargchk(1, 5, nargin));

% if nargin < 2
%     [Model_In, path] = uigetfile('.osim');
%     fileoutpath = [path Model_In(1:end-14),'_mass_adjusted.osim'];    
%     filepath = [path Model_In];
% elseif nargin < 3
%     fileoutpath = [Model_In(1:end-14),'_mass_adjusted.osim'];
%     filepath = Model_In;
% else
    filepath = Model_In;
    fileoutpath = Model_Out;
%end

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

% Get the set of bodies that are in the original model
Bodies1 = Model1.getBodySet(); 

%Count the bodies
nBodies = Bodies1.getSize();
disp(['Number of bodies in orginal model: ' num2str(nBodies)]);

% Get the set of bodies that are in the scaled model
% (Should be the same as the original at this point.)
Bodies2 = Model2.getBodySet();

% Read the adjusted mass suggestions from the RRA log file
% and read the pelvis width scale factor from the Scale log file
segment_mass = read_rra_mass(RAA_Log);
%wrap_scale_factor = read_pelvis_scale(SCALE_Log);

% loop through forces and scale muscle Fmax accordingly (index starts at 0)
for i = 0:nMuscles-1
        
    %get the muscle that the original muscle set points to
    %to read the muscle type and the max isometric force
    currentMuscle = Muscles1.get(i);
    
    %define the muscle in the modified model for changing
    newMuscle = Muscles2.get(i);

    %define the new muscle force by multiplying current muscle max
    %force by the scale factor
    newMuscle.setMaxIsometricForce(currentMuscle.getMaxIsometricForce()*scaleFactor);

end

% loop through bodies and scale body mass accordingly (index starts at 0)
for i = 0:nBodies-1
        
    %get the body that the original body set points to
    %to read the body type and the body mass
    currentBody = Bodies1.get(i);
    
    %define the body in the modified model for changing
    newBody = Bodies2.get(i);

    %define the new body mass by multiplying current body mass
    %by the scale factor
    %newBody.setMass(currentBody.getMass()*scaleFactor);
    newBody.setMass(segment_mass(i+1));
           
end

% sclae the domensions of the wrap object at the back joint relative to
% pelvis scale

%    wrapElip = WrapEllipsoid();
%    BodyRef = Bodies2.get('torso');
%    wrap = BodyRef.getWrapObject('BACK');
%    %wrap_back = wrapObjSet('BACK');
%    size = wrap.setDimensions();
%    new_wrap = wrap_back.setDimensions(size*wrap_scale_factor);
%    
%    wrapObjSet.cloneAndAppend(new_wrap);

% cyl = new WrapCylinder(); cyl.setRadius(.5); cyl.setLength(1.0);
% wrapCyl=modeling.WrapCylinder()  #Create a new wrap object
% wrapCyl.setAllPropertiesUseDefault(True)
% wrapCyl.setName('NewWrap')
% wrapCyl.setRadius(1)
% wrapCyl.setLength(1)
% wrapCyl.setQuadrantName('x')
% 
% BodyRef=model.getBodySet().get('NameOfBodyToAddTo')
% wrapObjSet=BodyRef.getWrapObjectSet()
% wrapObjSet.cloneAndAppend(wrapCyl)

% add five makers per foot as per required to run the IAA plugin
% set the offset
location_r = [0.01,0.00,-0.05; 0.00,0.00,0.03; 0.135,0.00,0.07; 0.205,0.00,-0.05; 0.275,0.00,0.02];
location_l = [0.01,0.00,0.05; 0.00,0.00,-0.03; 0.135,0.00,-0.07; 0.205,0.00,0.05; 0.275,0.00,-0.02];
% add the maker name convention
nfootmarker = ['A';'B';'C';'D';'E'];
% get the current marker set
markers = Model2.getMarkerSet();
% loop through and write the five markers to each foot
for i = 1:length(location_r)
newmarkers = markers.addMarker(['RIGHT_' nfootmarker(i)],[location_r(i,:)],Bodies2.get(5));
newmarkers = markers.addMarker(['LEFT_' nfootmarker(i)],[location_l(i,:)],Bodies2.get(10));
end


% save the updated model to an OSIM xml file
Model2.print(fileoutpath)
disp(['The new model has been saved at ' fileoutpath]);

end
