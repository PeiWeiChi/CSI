function data = btk_c3d2trc_DC(varargin)
% function btk_c3d2trc(file) OR
% function btk_c3d2trc(data)
%
% Function to convert data from a C3D file into the TRC and MOT file
% formats for OpenSim
%
% INPUT -   file - the C3D file path that you wish to load (leave blank to
%               choose from a dialog box) OR
%           data - structure containing fields from from previously loaded
%               C3D file using btk_loadc3d.m
%           anim - animate 'on' or 'off' (default - 'on')
%
% OUTPUT -  data - structure containing the relevant data from the c3dfile
%                  Creates the TRC file and _grf.MOT file for OpenSim
%
% example - data = btk_c3dtrc('filein.c3d','off');
%           data = btk_c3dtrc(data,'on');
%
% Written by Glen Lichtwark (University of Queensland)
% Updated September 2012

%% load data
if nargin > 0
    if ~isstruct(varargin{1})
    % load C3d file
        file = varargin{1};
        if isempty(fileparts(file))
            pname = cd;
            if ispc
                pname = [pname '\'];
            else pname = [pname '/'];
            end
            fname = file;
        else [pname, name, ext] = fileparts(file);
            fname = [name ext];
        end
        % load the c3dfile
        [data.marker_data, data.analog_data, data.fp_data, data.sub_info] = btk_loadc3d([pname, fname], 10);
        
    else
        data = varargin{1};
        if ~isfield(data,'marker_data') 
            error('Please ensure that the following field is included in structure - marker_data. Please use btk_loadc3d for correct outputs');
        end
        if isfield(data,'sub_info')
          
            
 if isfield(data.sub_info, 'Filename')
[pname, name, ext] = fileparts(data.sub_info.Filename);
else
   
[pname, name, ext] = fileparts(data.sub_info.Name); 
end
            if ispc
                pname = [pname '\'];
            else pname = [pname '/'];
            end
            fname = [name ext];
        else fname = data.marker_data.Filename;
        end
        
    end
    if length(varargin) < 2
        anim = 'on';
    else anim = varargin{2};
    end
else [fname, pname] = uigetfile('*.c3d', 'Select C3D file');
    % load the c3dfile
    [data.marker_data, data.analog_data, data.fp_data, data.sub_info] = btk_loadc3d([pname, fname], 10);
    anim = 'on';
end

%%
% try and find the mass and height of the subject from the MP file and
% store into the data structure
if ~isfield(data.sub_info,'Bodymass')
    data.Mass = 75; % default to 75kg
end

if ~isfield(data.sub_info,'Height')
    data.Height = 1800; % default to 1800mm
end

if ~isfield(data.sub_info,'Name')
    data.Name = 'NoName';
else data.Name = data.sub_info.Name;
end

%% define the start and end frame for analysis as first and last frame unless 
% this has already been done to change the analysed frames
if ~isfield(data,'Start_Frame')
    data.Start_Frame = 1;
    data.End_Frame = data.marker_data.Info.NumFrames;
end

%%
% define some parameters 
nrows = data.End_Frame-data.Start_Frame+1;
nmarkers = length(fieldnames(data.marker_data.Markers));

 if isfield(data.marker_data, 'Info')
data.time = (1/data.marker_data.Info.frequency:1/data.marker_data.Info.frequency:(data.End_Frame-data.Start_Frame+1)/data.marker_data.Info.frequency)';

 else
     
data.marker_data.Info=data.Info;     
data.time = (1/data.marker_data.Info.frequency:1/data.marker_data.Info.frequency:(data.End_Frame-data.Start_Frame+1)/data.marker_data.Info.frequency)';

end

nframe = 1:nrows;

% anim the trial if animation = on
if strcmp(anim,'on')
    data.marker_data.First_Frame = data.Start_Frame;
    data.marker_data.Last_Frame = data.End_Frame;
    if isfield(data,'fp_data')
        btk_animate_markers(data.marker_data, data.fp_data, 5)
    else btk_animate_markers(data.marker_data)
    end
end

%%
% we need to reorder the lab coordinate system to match that of the OpenSim
% system --> SKIP THIS STEP IF LAB COORDINATE SYSTEM IS SAME AS MODEL
% SYSTEM
markers = fieldnames(data.marker_data.Markers); % get markers names

if strcmp(data.marker_data.Info.units.ALLMARKERS,'mm')
    p_sc = 1000;
    data.marker_data.Info.units.ALLMARKERS = 'm';
else p_sc = 1;
end

% go through each marker field and re-order from X Y Z to Z X Y
for i = 1:nmarkers
% % % % %     data.marker_data.Markers.(markers{i}) =  [-1*(data.marker_data.Markers.(markers{i})(:,2))...
% % % % %         data.marker_data.Markers.(markers{i})(:,3) -1*(data.marker_data.Markers.(markers{i})(:,1))];
    data.marker_data.Markers.(markers{i}) =  [(data.marker_data.Markers.(markers{i})(:,2))...
        data.marker_data.Markers.(markers{i})(:,3) (data.marker_data.Markers.(markers{i})(:,1))]/p_sc;
end

%%
% now we need to make the headers for the column headings for the TRC file
% which are made up of the marker names and the XYZ for each marker

% first initialise the header with a column for the Frame # and the Time
% also initialise the format for the columns of data to be written to file
dataheader1 = 'Frame#\tTime\t';
dataheader2 = '\t\t';
format_text = '%i\t%2.4f\t';
% initialise the matrix that contains the data as a frame number and time row
data_out = [nframe; data.time'];

%%
% Write motion file containing GRFs

disp('Writing grf.mot file...')

if isfield(data,'fp_data')
F = data.fp_data.Info(1).frequency/data.marker_data.Info.frequency; % assume that all force plates are collected at the same frequency!!!

fp_time = 0:1/data.fp_data.Info(1).frequency:((F*data.End_Frame)-(F*data.Start_Frame))/data.fp_data.Info(1).frequency;

% initialise force data matrix with the time array and column header
force_data_out = fp_time';
force_header = 'time\t';
force_format = '%20.6f\t';

% go through each marker field and re-order from X Y Z to Y Z X and place
% into data array and add data to the force data matrix --> also need to
% divide by 1000 to convert to mm from m if necessary


for i = 1:length(data.fp_data.FP_data) 
    % reoder data so lab coordinate system to match that of the OpenSim
    % system
    data.fp_data.GRF_data(i).P =  [data.fp_data.GRF_data(i).P(:,1)/p_sc ...
       data.fp_data.GRF_data(i).P(:,2)/p_sc data.fp_data.GRF_data(i).P(:,3)/p_sc];
   data.fp_data.GRF_data(i).F =  [data.fp_data.GRF_data(i).F(:,1) ...
       data.fp_data.GRF_data(i).F(:,2) data.fp_data.GRF_data(i).F(:,3)];
%    data.fp_data.GRF_data(i).M =  [-1*(data.fp_data.GRF_data(i).M(:,2)) ...
%        data.fp_data.GRF_data(i).M(:,3) -1*(data.fp_data.GRF_data(i).M(:,1))];
%    data.fp_data.GRF_data(i).M =  [-0.001*(data.fp_data.GRF_data(i).M(:,2)) ...
%        0.001*(data.fp_data.GRF_data(i).M(:,3)) -0.001*(data.fp_data.GRF_data(i).M(:,1))];
    data.fp_data.GRF_data(i).M =  [data.fp_data.GRF_data(i).M(:,1)/p_sc ...
        data.fp_data.GRF_data(i).M(:,2)/p_sc data.fp_data.GRF_data(i).M(:,3)/p_sc]; 

   %*** GT addition to try and filter force data
   data.fp_data.GRF_data(i).P=Filtmat(0.001,12,2,data.fp_data.GRF_data(i).P);
   data.fp_data.GRF_data(i).F=Filtmat(0.001,12,2,data.fp_data.GRF_data(i).F);
   data.fp_data.GRF_data(i).M=Filtmat(0.001,12,2,data.fp_data.GRF_data(i).M);
   
   % do some cleaning of the COP before and after contact
   b = find(abs(diff(data.fp_data.GRF_data(i).P(:,3)))>0);
   if ~isempty(b)
       for j = 1:3
            data.fp_data.GRF_data(i).P(1:b(1),j) = data.fp_data.GRF_data(i).P(b(1)+1,j);
            data.fp_data.GRF_data(i).P(b(end):end,j) = data.fp_data.GRF_data(i).P(b(end)-1,j);
       end
   end
   
   % define the period which we are analysing
   K = (F*data.Start_Frame):1:(F*data.End_Frame);
   
   % add the force, COP and moment data for current plate to the force matrix 
   force_data_out = [force_data_out data.fp_data.GRF_data(i).F(K,:) data.fp_data.GRF_data(i).P(K,:) data.fp_data.GRF_data(i).M(K,:)];
   % define the header and formats
   force_header = [force_header num2str(i) '_ground_force_vx\t' num2str(i) '_ground_force_vy\t' num2str(i) '_ground_force_vz\t'...
       num2str(i) '_ground_force_px\t' num2str(i) '_ground_force_py\t' num2str(i) '_ground_force_pz\t' ...
       num2str(i) '_ground_torque_x\t' num2str(i) '_ground_torque_y\t' num2str(i) '_ground_torque_z\t'];

   force_format = [force_format '%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t'];
   
end

force_header = [force_header(1:end-2) '\n'];
force_format = [force_format(1:end-2) '\n'];

% assign a value of zero to any NaNs
force_data_out(logical(isnan(force_data_out))) = 0;

newfilename = [fname(1:end-4) '_grf.mot'];

data.GRF_Filename = [pname newfilename];

fid_2 = fopen([pname newfilename],'w');

% write the header information
fprintf(fid_2,'name %s\n',newfilename);
fprintf(fid_2,'datacolumns %d\n', size(force_data_out,2));  % total # of datacolumns
fprintf(fid_2,'datarows %d\n',length(fp_time)); % number of datarows
fprintf(fid_2,'range %f %f\n',data.marker_data.Time(1,1),data.marker_data.Time(end,1)); % range of time data
fprintf(fid_2,'endheader\n');
fprintf(fid_2,force_header);

% write the data
fprintf(fid_2,force_format,force_data_out');

fclose(fid_2);

disp('Done.')
else disp('No force plate information available.')
end
