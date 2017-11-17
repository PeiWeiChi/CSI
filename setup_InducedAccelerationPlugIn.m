function IndAccPI = setup_InducedAccelerationPlugIn(varargin)

% function IndAccPI = setup_InducedAccelerationPlugIn(data, varagin)
%
% Inspirations from Tim Dorn's Gait Extract Toolbox -writeXML.m (University of
% Melbourne)

% setup default input files to emtpy
ModelFile = [];
ReplaceForceSet = 'false';
ForceSetFiles = [];

%DirectoryName = 'IndAccPI_Results';
DirectoryName = [];
OutputPrecision = 20;
InitialTime = [];
FinalTime = [];

CoordinatesFile = [];
LowPassFilterFreq = [];

ExternalLoadsFile = [];
ExternalLoadsModelKinematicsFile = [];
LowPassFilterFreq = -1;

StepInterval = 1; 
KineticsFile = [];
ForcesFile = [];
ForceThreshold = 15;
FootpointMarkers = ['RIGHT_A ' 'RIGHT_B ' 'RIGHT_C ' 'RIGHT_D ' 'RIGHT_E ' 'LEFT_A ' 'LEFT_B ' 'LEFT_C ' 'LEFT_D ' 'LEFT_E '];

data = [];

if ~isempty(varargin)
    if rem(length(varargin),2)
        error('Incorrect input arguments - must specify property and input')
    end
    for i = 1:2:length(varargin)
       n = varargin{i+1};
       eval([varargin{i} '= n;']); 
    end    
end

% define the initial and final times for inverse dynamics from the data structure 
% if this is passed to function and these aren't prescribed otherwise
if ~isempty(data)
    if isfield(data,'time')
        if isempty(InitialTime)
            InitialTime = data.time(1);
        end
        if isempty(FinalTime)
            FinalTime = data.time(end);
        end
    end
end

if isempty(CoordinatesFile)
    error('CoordinatesFile must be included for RRA analysis')
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

V.AnalyzeTool.ATTRIBUTE.name = [Model '_' trial '_IndAccPI'];

% define the model file
if ~isempty(ModelFile)
    V.AnalyzeTool.model_file = ModelFile;
else error('Please specify a model file')
end

% define the force set and determine whether this is to replace or
% append to current actuator set
V.AnalyzeTool.replace_force_set = ReplaceForceSet;
V.AnalyzeTool.force_set_files = ForceSetFiles;

% define results directory and precision
%V.AnalyzeTool.results_directory = ['./' DirectoryName];
V.AnalyzeTool.results_directory = DirectoryName;
V.AnalyzeTool.output_precision = OutputPrecision;

% Define times to perform analysis over 
V.AnalyzeTool.initial_time = num2str(InitialTime,6);
V.AnalyzeTool.final_time = num2str(FinalTime,6);

% Kinematics and kinetics foles
V.AnalyzeTool.coordinates_file = CoordinatesFile;
V.AnalyzeTool.lowpass_cutoff_frequency_for_coordinates = LowPassFilterFreq;
V.AnalyzeTool.external_loads_file = ExternalLoadsFile;
V.AnalyzeTool.external_loads_model_kinematics_file = ExternalLoadsModelKinematicsFile;
V.AnalyzeTool.lowpass_cutoff_frequency_for_load_kinematics = LowPassFilterFreq;

V.AnalyzeTool.AnalysisSet.ATTRIBUTE.name = 'Analyses';
V.AnalyzeTool.AnalysisSet.objects.IndAccPI.ATTRIBUTE.name = 'IndAccPI';

V.AnalyzeTool.AnalysisSet.objects.IndAccPI.on = 'true';
V.AnalyzeTool.AnalysisSet.objects.IndAccPI.step_interval = StepInterval;
V.AnalyzeTool.AnalysisSet.objects.IndAccPI.in_degrees = 'true';
V.AnalyzeTool.AnalysisSet.objects.IndAccPI.kinetics_file = KineticsFile;
V.AnalyzeTool.AnalysisSet.objects.IndAccPI.forces_file = ForcesFile;
V.AnalyzeTool.AnalysisSet.objects.IndAccPI.force_threshold = ForceThreshold;
V.AnalyzeTool.AnalysisSet.objects.IndAccPI.weights = ['10000 ' '100 ' '1'];
V.AnalyzeTool.AnalysisSet.objects.IndAccPI.footpoint_markers = FootpointMarkers;
V.AnalyzeTool.AnalysisSet.objects.IndAccPI.accleration_match = 'true';
V.AnalyzeTool.AnalysisSet.objects.IndAccPI.output_grf_contribution = 2;

fileout = [data.Name '_Setup_InducedAccelerationPlugIn.xml'];

Pref.StructItem = false;

xml_write(fileout, V, root,Pref);

IndAccPI = V;