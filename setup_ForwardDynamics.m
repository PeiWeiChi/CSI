function FD = setup_ForwardDynamics(varargin)

% function FD = setup_ForwardDynamics(data, varagin)
%
% This function will create a Setup_FD.XML file for use with OpenSIM
% based on the data from the C3D file and the names of the model files
% which will be used to scale the model.
% 
% Input - 'ModelFile' - string which is the filename (including path) of the
%                OSIM file (model file)      
%         'MOTFile' - desired kinematics file MOT file for ID
%         'GRFFile' - filename string of the XML file containing GRF
%               information 
%
%         OPTIONAL RECOMMENDED parameters
%         'data' - structure containing data from C3D file after processing
%                with C3D2TRC.m and any other processing (e.g. IK)
%         'CMCTaskFile' - filename string of the Tasks XML file
%         'CMCForceFile' - filename string of the Actuator XML file
%         'CMCConstraintsFile' - File containing the constraints on the
%               controls.
%         'RRAControlsFile' - File containing the controls output by RRA. 
%               These can be used to place constraints on the residuals during CMC.
%
%         OPTIONAL parameters
%         'OutputPrecision' - number between 1-50 (default 8)
%         'LowPassFilterFreq' - number between 1 and 60 for low pass filter 
%                       of kinematics (default -1 = none)
%         'DirectoryName' - string which is the directory name to be made
%                           for the results (default 'CMCResults'
%         'ReplaceForceSet' - 'true' or 'false' to deteremine whether 
%                          the model actuator set is replaced with those
%                          from actuator file (default 'true' -replaced)
%         'InitialTime' - initial time for ID to run from (defaults to
%                       start of the MOT file time)
%         'FinalTime' - final time for ID to run until (defaults to
%                       time at end of the MOT file)
%         'AuxillaryStates' - whether or not to compute equilibrium values for
%       		    states other than the coordinates or speeds (default
%       		    - 'false')
%         'MaxIntegratorSteps' - maximum number of integrator steps
%                   (default=20000)
%         'MaxStepSize' - maximum integrator step size (default=1)
%         'MinStepSize' - maximum integrator step size (default=1e-6)
%         'ErrorTol' - inegrator error tolerance (default=1e-5)
%         'OptimizerAlgorithm' - 'ipopt' or 'cfsqp' (default - 'ipopt')
%         'OptimizerDx' - optimizer derivative dx value (default = 1e-4)
%         'OptimConvergCrit' - optmiser convergence criterion value
%               (default = 1e-6)
%         'OptimMaxIterations' - Maximum number of iterations for the optimize
%               (default - 30000)
%         'CMCTimeWindow' - Time window over which the desired actuator
%               forces are achieved (default - 0.01);
%         'UseCurvatureFilter' - Flag (true or false) indicating whether or 
%               not to use the curvature filter. Setting this flag to true
%               can reduce oscillations in the computed muscle excitations
%               (default - 'false')	    
%         'FastTarget' - Flag (true or false) indicating whether to    
%		    use the fast CMC optimization target. The fast target requires 
%		    the desired accelerations to be met (default - 'false')
%
% Output - CMC (optional) - This is the structure which is used to make the
%                         XML output file
%
% E.g.  CMC = setup_CMC('data',data, 'ModelFile',ModelFile,'MOTFile',MOTFile,'GRFFile',GRFFile,...
%        'CMCTaskFile',CMCTaskFile,'CMCActuatorFile',CMCActuatorFile,'CMCControlFile',CMCControlFile,...
%        'AdjCOMRes','true','OptimMaxIter',20000,'LowPassFilterFreqLoad',-1,...
%        'LowPassFilterFreq',-1);
% 
% Written by Glen Lichtwark (The University of Queensland)
%
% Inspirations from Tim Dorn's Gait Extract Toolbox -writeXML.m (University of
% Melbourne)

% setup default input files to emtpy
ModelFile = [];

ReplaceForceSet = 'false';
FDForceFile = [];

DirectoryName = 'FDResults';
OutputPrecision = 12;
InitialTime = [];
FinalTime = [];

AuxillaryStates = 'true';
MaxIntegratorSteps = 20000;
MaxStepSize = 0.00025;
MinStepSize = 1e-6;
ErrorTol = 1e-5;
Specified_dt = 'true';

GRFFile = [];
MOTFile = [];

FDForceFile = [];
FDStatesFile = [];
%FDConstraintsFile = [];
FDControlsFile = ''; 

kinematics_on = 'true';
kinematics_step = 10;
kinematics_in_degrees = 'true';

actuation_on = 'true';
actuation_step = 10;
actuation_in_degrees = 'true';

bodykinematics_on = 'true';
bodykinematics_step = 10;
bodykinematics_in_degrees = 'true';

%actuator_list = ;
enable_controller = 'ALL';

LowPassFilterFreq = -1;

FDTimeWindow = 0.01;
UseCurvatureFilter = 'true';
FastTarget = 'true';

OptimizerAlgorithm = 'ipopt';
OptimizerDx = 1e-4;
OptimConvergCrit = 1e-6;
OptimMaxIterations = 30000;

data = [];
CMCdata = [];

if ~isempty(varargin)
    if rem(length(varargin),2)
        error('Incorrect input arguments - must specify property and input')
    end
    for i = 1:2:length(varargin)
       n = varargin{i+1};
       eval([varargin{i} '= n;']); 
    end    
end

% define the initial and final times for forward dynamics from the data structure 
% if this is passed to function and these aren't prescribed otherwise
if ~isempty(data)
    if isfield(data,'time')
        if isempty(InitialTime)
            InitialTime = CMCdata.CMC.States.time(1);
        end
        if isempty(FinalTime)
            FinalTime = CMCdata.CMC.States.time(end);
        end
    end
end

%setup root and RRATool
root = 'OpenSimDocument';

V.ATTRIBUTE.Version = '30000';
% define some names for outputting... use the data in the data structure to
% limit the filename size to important parts if data structure is passed
if ~isempty(data)
    Model = data.Name;
    [~,trial, ~] = fileparts(data.TRC_Filename);
else [~, Model, ~] = fileparts(ModelFile);
    [~, trial, ~] = fileparts(MOTFile);
end

V.ForwardTool.ATTRIBUTE.name = [Model '_' trial '_FD'];

% define the model file
if ~isempty(ModelFile)
    V.ForwardTool.model_file = ModelFile;
else error('Please specify a model file')
end

% define the force set and determine whether this is to replace or
% append to current actuator set
V.ForwardTool.force_set_files = FDForceFile;
V.ForwardTool.replace_force_set = ReplaceForceSet;

% Define times to perform analysis over 
V.ForwardTool.initial_time = num2str(InitialTime,12);
V.ForwardTool.final_time = num2str(FinalTime,12);


% define integrator settings
V.ForwardTool.maximum_number_of_integrator_steps = MaxIntegratorSteps;
V.ForwardTool.maximum_integrator_step_size = MaxStepSize;
V.ForwardTool.minimum_integrator_step_size = MinStepSize;
V.ForwardTool.integrator_error_tolerance = ErrorTol;
V.ForwardTool.use_specified_dt = Specified_dt;

% define results directory and precision
V.ForwardTool.results_directory = DirectoryName;
V.ForwardTool.output_precision = OutputPrecision;

V.ForwardTool.AnalysisSet.ATTRIBUTE.name = 'Analyses';
V.ForwardTool.AnalysisSet.objects.Kinematics.ATTRIBUTE.name = 'Kinematics';
V.ForwardTool.AnalysisSet.objects.Kinematics.on = kinematics_on;
V.ForwardTool.AnalysisSet.objects.Kinematics.step_interval = kinematics_step;
V.ForwardTool.AnalysisSet.objects.Kinematics.in_degrees = kinematics_in_degrees;

V.ForwardTool.AnalysisSet.objects.Actuation.ATTRIBUTE.name = 'Actuation';
V.ForwardTool.AnalysisSet.objects.Actuation.on = actuation_on;
V.ForwardTool.AnalysisSet.objects.Actuation.step_interval = actuation_step;
V.ForwardTool.AnalysisSet.objects.Actuation.in_degrees = actuation_in_degrees;

V.ForwardTool.AnalysisSet.objects.BodyKinematics.ATTRIBUTE.name = 'BodyKinematics';
V.ForwardTool.AnalysisSet.objects.BodyKinematics.on = actuation_on;
V.ForwardTool.AnalysisSet.objects.BodyKinematics.step_interval = actuation_step;
V.ForwardTool.AnalysisSet.objects.BodyKinematics.in_degrees = actuation_in_degrees;

% define input and output files for FD
V.ForwardTool.states_file = FDStatesFile;

V.ForwardTool.ControllerSet.ATTRIBUTE.name = 'Controllers';
V.ForwardTool.ControllerSet.objects.ControlSetController.ATTRIBUTE.name = '';
%V.ForwardTool.ControllerSet.objects.ControlsetController.actuator_list;
V.ForwardTool.ControllerSet.objects.ControlSetController.enable_controller = enable_controller;
V.ForwardTool.ControllerSet.objects.ControlSetController.controls_file = FDControlsFile;

% External Loads
V.ForwardTool.external_loads_file = GRFFile;

V.ForwardTool.optimizer_print_level = 0;
V.ForwardTool.use_verbose_printing = 'false';

fileout = [Model '_Setup_ForwardDynamics.xml'];

Pref.StructItem = false;

xml_write(fileout, V, root,Pref);

FD = V;
        