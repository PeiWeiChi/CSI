function data = rugby_scrum_pipeline_v1(file)

% This function will create the appropriate files to run an OpenSim
% simulation sequence for rugby scrum data
% Input - file - c3d file to process
% Originally written by Glen Lichtwark (The University of Queensland)
% Copyright 2011
% Contributions by David Graham and Chris Carty (Griffith University) 2012
%
% Dario Cazzola implemented rugby_scrum version
%  DC Included:
%  -) Scrum forces CoP calculation via pressure sensors;
%  -) Joint center calculation (for static)


% First, import the classes from the jar file so that these can be called
% directly
import org.opensim.modeling.*

%% Define path and model

osim_path = 'C:\Users\install\Desktop\CSI\CSI_Scrum_2015\ISBS2015 - Workshop\OpenSim_Files_Workshop_TBE\';

model = 'Rugby_Model';


global subdir_current;
global trialdir_current;


if nargin < 1
    [fname, pname] = uigetfile('*.c3d', 'Select C3D file');
else
    if isempty(fileparts(file))
        pname = cd;
        pname = [pname '\'];
        fname = file;
    else [pname, name, ext] = fileparts(file);
        fname = [name ext];
    end
end

cd(pname);

%% Load the c3dfile

% load the c3d file using BTK
data = btk_loadc3d([pname, fname], 5);

%% Naming in case of Static cal files
if ~isempty(strfind(lower(fname),'cal'))
    data.sub_info.Name = fname(1:end-14);
    data.Name = fname(1:end-14);
else
    data.sub_info.Name = fname(7:end-4);
    data.Name = fname(7:end-4);
end

%% BTK uses First/Last so convert to fit existing routine so Start/End
data.marker_data.Last_Frame = data.marker_data.Last_Frame - data.marker_data.First_Frame;
data.marker_data.First_Frame = 1;

data.Start_Frame = data.marker_data.First_Frame;
data.End_Frame = data.marker_data.Last_Frame;

% % % % % %% JOINT CENTERS  %%
% % % % % 
% % % % % if ~isempty(strfind(lower(fname),'cal'))
% % % % %     
% % % % %     %     ------ HIP ---------
% % % % %     GTs_dist_l=([data.marker_data.Markers.GT_R(:,1) data.marker_data.Markers.GT_R(:,2) data.marker_data.Markers.GT_R(:,3)]...
% % % % %         - [data.marker_data.Markers.GT_L(:,1) data.marker_data.Markers.GT_L(:,2) data.marker_data.Markers.GT_L(:,3)])/4;
% % % % %     GTs_dist_r=([data.marker_data.Markers.GT_L(:,1) data.marker_data.Markers.GT_L(:,2) data.marker_data.Markers.GT_L(:,3)]...
% % % % %         - [data.marker_data.Markers.GT_R(:,1) data.marker_data.Markers.GT_R(:,2) data.marker_data.Markers.GT_R(:,3)])/4;
% % % % %     
% % % % %     HJC_R= ([data.marker_data.Markers.GT_R(:,1) data.marker_data.Markers.GT_R(:,2) data.marker_data.Markers.GT_R(:,3)]...
% % % % %         + [GTs_dist_r(:,1) GTs_dist_r(:,2) GTs_dist_r(:,3)]);
% % % % %     
% % % % %     HJC_L= ([data.marker_data.Markers.GT_L(:,1) data.marker_data.Markers.GT_L(:,2) data.marker_data.Markers.GT_L(:,3)]...
% % % % %         + [GTs_dist_l(:,1) GTs_dist_l(:,2) GTs_dist_l(:,3)]);
% % % % %     
% % % % %     GT_PSIS_L_dir= ([HJC_L(:,1) HJC_L(:,2) HJC_L(:,3)]...
% % % % %         - [data.marker_data.Markers.PSIS_L(:,1) data.marker_data.Markers.PSIS_L(:,2) data.marker_data.Markers.PSIS_L(:,3)]);
% % % % %     
% % % % %     GT_PSIS_R_dir= ([HJC_R(:,1) HJC_R(:,2) HJC_R(:,3)]...
% % % % %         - [data.marker_data.Markers.PSIS_R(:,1) data.marker_data.Markers.PSIS_R(:,2) data.marker_data.Markers.PSIS_R(:,3)]);
% % % % %     
% % % % %     data.marker_data.Markers.HJC_R=[];
% % % % %     data.marker_data.Markers.HJC_R=HJC_R;
% % % % %     
% % % % %     data.marker_data.Markers.HJC_L=[];
% % % % %     data.marker_data.Markers.HJC_L=HJC_L;
% % % % %     
% % % % %     save ('GT_PSIS_L_dir.mat','GT_PSIS_L_dir');
% % % % %     save ('GT_PSIS_R_dir.mat','GT_PSIS_R_dir');
% % % % %     
% % % % %     
% % % % %     % ------ KNEE ---------
% % % % %     
% % % % %     %%% Create new Marker in the model to represent the joint center
% % % % %     
% % % % %     if isfield(data.marker_data.Markers,'MC_R'),
% % % % %         KJC_r=([data.marker_data.Markers.MC_R(:,1) data.marker_data.Markers.MC_R(:,2) data.marker_data.Markers.MC_R(:,3)]...
% % % % %             + [data.marker_data.Markers.LC_R(:,1) data.marker_data.Markers.LC_R(:,2) data.marker_data.Markers.LC_R(:,3)])/2;
% % % % %         
% % % % %         KJC_r_dir= ([data.marker_data.Markers.MC_R(:,1) data.marker_data.Markers.MC_R(:,2) data.marker_data.Markers.MC_R(:,3)]...
% % % % %             - [data.marker_data.Markers.LC_R(:,1) data.marker_data.Markers.LC_R(:,2) data.marker_data.Markers.LC_R(:,3)])/2;
% % % % %         
% % % % %         save ('KJC_r_dir.mat','KJC_r_dir');
% % % % %         
% % % % %         data.marker_data.Markers.KJC_R=[];
% % % % %         data.marker_data.Markers.KJC_R=KJC_r;
% % % % %         
% % % % %         % Left
% % % % %         KJC_l=([data.marker_data.Markers.MC_L(:,1) data.marker_data.Markers.MC_L(:,2) data.marker_data.Markers.MC_L(:,3)]...
% % % % %             + [data.marker_data.Markers.LC_L(:,1) data.marker_data.Markers.LC_L(:,2) data.marker_data.Markers.LC_L(:,3)])/2;
% % % % %         
% % % % %         KJC_l_dir= ([data.marker_data.Markers.MC_L(:,1) data.marker_data.Markers.MC_L(:,2) data.marker_data.Markers.MC_L(:,3)]...
% % % % %             - [data.marker_data.Markers.LC_L(:,1) data.marker_data.Markers.LC_L(:,2) data.marker_data.Markers.LC_L(:,3)])/2;
% % % % %         
% % % % %         save ('KJC_l_dir.mat','KJC_l_dir');
% % % % %         
% % % % %         data.marker_data.Markers.KJC_L=[];
% % % % %         data.marker_data.Markers.KJC_L=KJC_l;
% % % % %     else % MC and MM markers missing in static
% % % % %         KJC_l=data.marker_data.Markers.LC_L(:,:);
% % % % %         KJC_r=data.marker_data.Markers.LC_R(:,:);
% % % % %         KJC_l_dir=zeros(length(data.marker_data.Markers.LC_L(:,1)),3);
% % % % %         KJC_r_dir=zeros(length(data.marker_data.Markers.LC_R(:,1)),3);
% % % % %         data.marker_data.Markers.KJC_R=[];
% % % % %         data.marker_data.Markers.KJC_R=KJC_r;
% % % % %         data.marker_data.Markers.KJC_L=[];
% % % % %         data.marker_data.Markers.KJC_L=KJC_l;
% % % % %         save ('KJC_r_dir.mat','KJC_r_dir');
% % % % %         save ('KJC_l_dir.mat','KJC_l_dir');
% % % % %     end
% % % % %     % ------ ANKLE --------- %
% % % % %     if isfield(data.marker_data.Markers,'MM_R'),
% % % % %         % Right
% % % % %         AJC_r=([data.marker_data.Markers.MM_R(:,1) data.marker_data.Markers.MM_R(:,2) data.marker_data.Markers.MM_R(:,3)]...
% % % % %             + [data.marker_data.Markers.LM_R(:,1) data.marker_data.Markers.LM_R(:,2) data.marker_data.Markers.LM_R(:,3)])/2;
% % % % %         
% % % % %         AJC_r_dir=([data.marker_data.Markers.MM_R(:,1) data.marker_data.Markers.MM_R(:,2) data.marker_data.Markers.MM_R(:,3)]...
% % % % %             - [data.marker_data.Markers.LM_R(:,1) data.marker_data.Markers.LM_R(:,2) data.marker_data.Markers.LM_R(:,3)])/2;
% % % % %         
% % % % %         save ('AJC_r_dir.mat','AJC_r_dir');
% % % % %         
% % % % %         data.marker_data.Markers.AJC_R=[];
% % % % %         data.marker_data.Markers.AJC_R=AJC_r;
% % % % %         
% % % % %         % Left
% % % % %         AJC_l=([data.marker_data.Markers.MM_L(:,1) data.marker_data.Markers.MM_L(:,2) data.marker_data.Markers.MM_L(:,3)]...
% % % % %             + [data.marker_data.Markers.LM_L(:,1) data.marker_data.Markers.LM_L(:,2) data.marker_data.Markers.LM_L(:,3)])/2;
% % % % %         
% % % % %         AJC_l_dir=([data.marker_data.Markers.MM_L(:,1) data.marker_data.Markers.MM_L(:,2) data.marker_data.Markers.MM_L(:,3)]...
% % % % %             - [data.marker_data.Markers.LM_L(:,1) data.marker_data.Markers.LM_L(:,2) data.marker_data.Markers.LM_L(:,3)])/2;
% % % % %         
% % % % %         save ('AJC_l_dir.mat','AJC_l_dir');
% % % % %         
% % % % %         data.marker_data.Markers.AJC_L=[];
% % % % %         data.marker_data.Markers.AJC_L=AJC_l;
% % % % %     else
% % % % %         AJC_l=data.marker_data.Markers.LM_L(:,:);
% % % % %         AJC_r=data.marker_data.Markers.LM_R(:,:);
% % % % %         AJC_l_dir=zeros(length(data.marker_data.Markers.LM_L(:,1)),3);
% % % % %         AJC_r_dir=zeros(length(data.marker_data.Markers.LM_R(:,1)),3);
% % % % %         data.marker_data.Markers.AJC_R=[];
% % % % %         data.marker_data.Markers.AJC_R=AJC_r;
% % % % %         data.marker_data.Markers.AJC_L=[];
% % % % %         data.marker_data.Markers.AJC_L=AJC_l;
% % % % %         save ('AJC_r_dir.mat','AJC_r_dir');
% % % % %         save ('AJC_l_dir.mat','AJC_l_dir');
% % % % %     end
% % % % %     
% % % % %     % ------ FOOT --------- %
% % % % %     % Right - Middle point betwwen CM1 and CM5
% % % % %     if isfield(data.marker_data.Markers,'C1M_R'),
% % % % %         CM5_1_R= ([data.marker_data.Markers.C1M_R(:,1) data.marker_data.Markers.C1M_R(:,2) data.marker_data.Markers.C1M_R(:,3)]...
% % % % %             + [data.marker_data.Markers.C5M_R(:,1) data.marker_data.Markers.C5M_R(:,2) data.marker_data.Markers.C5M_R(:,3)])/2;
% % % % %         
% % % % %         CM5_1_R_dir= ([data.marker_data.Markers.C1M_R(:,1) data.marker_data.Markers.C1M_R(:,2) data.marker_data.Markers.C1M_R(:,3)]...
% % % % %             - [data.marker_data.Markers.C5M_R(:,1) data.marker_data.Markers.C5M_R(:,2) data.marker_data.Markers.C5M_R(:,3)])/2;
% % % % %         
% % % % %         data.marker_data.Markers.CM5_1_R=[];
% % % % %         data.marker_data.Markers.CM5_1_R=CM5_1_R;
% % % % %         
% % % % %         save ('CM5_1_R_dir.mat','CM5_1_R_dir');
% % % % %         % Left - Middle point betwwen CM1 and CM5
% % % % %         
% % % % %         CM5_1_L= ([data.marker_data.Markers.C1M_L(:,1) data.marker_data.Markers.C1M_L(:,2) data.marker_data.Markers.C1M_L(:,3)]...
% % % % %             + [data.marker_data.Markers.C5M_L(:,1) data.marker_data.Markers.C5M_L(:,2) data.marker_data.Markers.C5M_L(:,3)])/2;
% % % % %         
% % % % %         CM5_1_L_dir= ([data.marker_data.Markers.C1M_L(:,1) data.marker_data.Markers.C1M_L(:,2) data.marker_data.Markers.C1M_L(:,3)]...
% % % % %             - [data.marker_data.Markers.C5M_L(:,1) data.marker_data.Markers.C5M_L(:,2) data.marker_data.Markers.C5M_L(:,3)])/2;
% % % % %         
% % % % %         data.marker_data.Markers.CM5_1_L=[];
% % % % %         data.marker_data.Markers.CM5_1_L=CM5_1_L;
% % % % %         
% % % % %         save ('CM5_1_L_dir.mat','CM5_1_L_dir');
% % % % %     else
% % % % %         CM5_1_l=data.marker_data.Markers.C5M_L(:,:);
% % % % %         CM5_1_r=data.marker_data.Markers.C5M_R(:,:);
% % % % %         CM5_1_L_dir=zeros(length(data.marker_data.Markers.C5M_L(:,1)),3);
% % % % %         CM5_1_R_dir=zeros(length(data.marker_data.Markers.C5M_R(:,1)),3);
% % % % %         data.marker_data.Markers.CM5_1_R=[];
% % % % %         data.marker_data.Markers.CM5_1_R=CM5_1_r;
% % % % %         data.marker_data.Markers.CM5_1_L=[];
% % % % %         data.marker_data.Markers.CM5_1_L=CM5_1_l;
% % % % %         save ('CM5_1_R_dir.mat','CM5_1_R_dir');
% % % % %         save ('CM5_1_L_dir.mat','CM5_1_L_dir');
% % % % %     end
% % % % %     
% % % % %     
% % % % %     
% % % % %     %% WeldConstraints %%
% % % % %     if isfield(data.marker_data.Markers,'C1M_R'),
% % % % %         dist_r=[data.marker_data.Markers.C1M_R(:,2)- data.marker_data.Markers.C5M_R(:,2),data.marker_data.Markers.C1M_R(:,3)- data.marker_data.Markers.C5M_R(:,3),data.marker_data.Markers.C1M_R(:,1)- data.marker_data.Markers.C5M_R(:,1)];
% % % % %         dist_l=[data.marker_data.Markers.C1M_L(:,2)- data.marker_data.Markers.C5M_L(:,2),data.marker_data.Markers.C1M_L(:,3)- data.marker_data.Markers.C5M_L(:,3),data.marker_data.Markers.C1M_L(:,1)- data.marker_data.Markers.C5M_L(:,1)];
% % % % %     else
% % % % %         dist_r=zeros(length(data.marker_data.Markers.C5M_L(:,1)),3);
% % % % %         dist_l=zeros(length(data.marker_data.Markers.C5M_L(:,1)),3);
% % % % %     end
% % % % %     
% % % % %     
% % % % %     save ('dist_r.mat','dist_r');
% % % % %     save ('dist_l.mat','dist_l');
% % % % %     
% % % % % 
% % % % %     
% % % % %     
% % % % % else
% % % % %     
% % % % %     % ------ HIP --------- %
% % % % %     load('GT_PSIS_L_dir.mat');
% % % % %     load('GT_PSIS_R_dir.mat');
% % % % %     
% % % % %     HJC_r=[(data.marker_data.Markers.PSIS_R(:,1)+ mean(GT_PSIS_R_dir(:,1))) (data.marker_data.Markers.PSIS_R(:,2)+ mean(GT_PSIS_R_dir(:,2))) ...
% % % % %         (data.marker_data.Markers.PSIS_R(:,3)+ mean(GT_PSIS_R_dir(:,3)))];
% % % % %     
% % % % %     data.marker_data.Markers.HJC_R=[];
% % % % %     data.marker_data.Markers.HJC_R=HJC_r;
% % % % %     
% % % % %     HJC_l=[(data.marker_data.Markers.PSIS_L(:,1)+ mean(GT_PSIS_L_dir(:,1))) (data.marker_data.Markers.PSIS_L(:,2)+ mean(GT_PSIS_L_dir(:,2))) ...
% % % % %         (data.marker_data.Markers.PSIS_L(:,3)+ mean(GT_PSIS_L_dir(:,3)))];
% % % % %     
% % % % %     data.marker_data.Markers.HJC_L=[];
% % % % %     data.marker_data.Markers.HJC_L=HJC_l;
% % % % %     
% % % % %     
% % % % %     % ------ KNEE --------- %
% % % % %     
% % % % %     %%% Create new Marker in the model to represent the joint center
% % % % %     % Right
% % % % %     load('KJC_r_dir.mat');
% % % % %     
% % % % %     KJC_r=[(data.marker_data.Markers.LC_R(:,1)+ mean(KJC_r_dir(:,1))) (data.marker_data.Markers.LC_R(:,2)+ mean(KJC_r_dir(:,2))) ...
% % % % %         (data.marker_data.Markers.LC_R(:,3)+ mean(KJC_r_dir(:,3)))];
% % % % %     
% % % % %     data.marker_data.Markers.KJC_R=[];
% % % % %     data.marker_data.Markers.KJC_R=KJC_r;
% % % % %     
% % % % %     % Left
% % % % %     load ('KJC_l_dir.mat');
% % % % %     KJC_l=[(data.marker_data.Markers.LC_L(:,1)+ mean(KJC_l_dir(:,1))) (data.marker_data.Markers.LC_L(:,2)+ mean(KJC_l_dir(:,2))) ...
% % % % %         (data.marker_data.Markers.LC_L(:,3)+ mean(KJC_l_dir(:,3)))];
% % % % %     
% % % % %     data.marker_data.Markers.KJC_L=[];
% % % % %     data.marker_data.Markers.KJC_L=KJC_l;
% % % % %     
% % % % %     % ------ ANKLE ---------
% % % % %     % Right
% % % % %     load ('AJC_r_dir.mat');
% % % % %     AJC_r=[(data.marker_data.Markers.LM_R(:,1)+ mean(AJC_r_dir(:,1))) (data.marker_data.Markers.LM_R(:,2)+ mean(AJC_r_dir(:,2))) ...
% % % % %         (data.marker_data.Markers.LM_R(:,3)+ mean(AJC_r_dir(:,3)))];
% % % % %     
% % % % %     data.marker_data.Markers.AJC_R=[];
% % % % %     data.marker_data.Markers.AJC_R=AJC_r;
% % % % %     
% % % % %     % Left
% % % % %     load ('AJC_l_dir.mat');
% % % % %     AJC_l=[(data.marker_data.Markers.LM_L(:,1)+ mean(AJC_l_dir(:,1))) (data.marker_data.Markers.LM_L(:,2)+ mean(AJC_l_dir(:,2))) ...
% % % % %         (data.marker_data.Markers.LM_L(:,3)+ mean(AJC_l_dir(:,3)))];
% % % % %     
% % % % %     data.marker_data.Markers.AJC_L=[];
% % % % %     data.marker_data.Markers.AJC_L=AJC_l;
% % % % %     
% % % % %     % ------ FOOT --------- %
% % % % %     % Right
% % % % %     load('CM5_1_R_dir.mat');
% % % % %     
% % % % %     CM5_1_R=[(data.marker_data.Markers.C5M_R(:,1)+ mean(CM5_1_R_dir(:,1))) (data.marker_data.Markers.C5M_R(:,2)+ mean(CM5_1_R_dir(:,2))) ...
% % % % %         (data.marker_data.Markers.C5M_R(:,3)+ mean(CM5_1_R_dir(:,3)))];
% % % % %     
% % % % %     data.marker_data.Markers.CM5_1_R=[];
% % % % %     data.marker_data.Markers.CM5_1_R=CM5_1_R;
% % % % %     
% % % % %     % Left
% % % % %     load('CM5_1_L_dir.mat');
% % % % %     
% % % % %     CM5_1_L=[(data.marker_data.Markers.C5M_L(:,1)+ mean(CM5_1_L_dir(:,1))) (data.marker_data.Markers.C5M_L(:,2)+ mean(CM5_1_L_dir(:,2))) ...
% % % % %         (data.marker_data.Markers.C5M_L(:,3)+ mean(CM5_1_L_dir(:,3)))];
% % % % %     
% % % % %     data.marker_data.Markers.CM5_1_L=[];
% % % % %     data.marker_data.Markers.CM5_1_L=CM5_1_L;
% % % % %     
% % % % %     
% % % % % end


%%

% sort the C3D file so we know what is Marker data and what is calculated
% data - Remember that AC_L, AC_R, SC_L and SC_R are missing during dynamic
% trials.

if ~isempty(strfind(lower(fname),'cal'))
    marker_names = {'HD_1';'HD_2';'HD_3';... %Head
        'C7';'T8';'L4';... % Spine
        'SC_L1';'SC_L2';'SC_L3';... % Scap Left
        'SC_R1';'SC_R2';'SC_R3';... % Scap Right
        'AC_L';'SC_L';... % Clavl L
        'AC_R';'SC_R'; 'ClavM';... % Clav
        'PSIS_L';'PSIS_R';... % Pelvis
        'GT_L'; 'LC_L'; 'MC_L'; 'LM_L'; 'MM_L';'IC_L';'TH1_L';'TH2_L';'TH3_L';'TH4_L'; 'SK1_L';'SK2_L';'SK3_L';'SK4_L';...% Left Leg
        'GT_R'; 'LC_R'; 'MC_R'; 'LM_R'; 'MM_R';'IC_R';'TH1_R';'TH2_R';'TH3_R';'TH4_R'; 'SK1_R';'SK2_R';'SK3_R';'SK4_R';...% Right Leg
        'B2_1';'B2_2';'B2_3';'B2_4';'B3_1';'B3_2';'B3_3';'B3_4';'ARM2';'ARM3';... % Pads markers
% % % % %         'HJC_R';'HJC_L';'KJC_R';'AJC_R';'KJC_L';'AJC_L'; ...% Joint Centers Markers
        'EL_L';'EL_R';'H1_L';'H2_L';'H3_L';'H1_R';'H2_R';'H3_R'; ... % Arms
        'HL_L';'HL_R';'C5M_R'; 'C1M_R';'C5M_L'; 'C1M_L';'CM5_1_L'; 'CM5_1_R';}; % Foot
    
else
    marker_names = {'HD_1';'HD_2';'HD_3';... %Head
        'C7';'T8';'L4';... % Spine
        'ClavM';'ClavL';'ClavR';... % Clavicle
        'SC_L1';'SC_L2';'SC_L3';... % Scap Left
        'SC_R1';'SC_R2';'SC_R3';... % Scap Right
        'PSIS_L';'PSIS_R';... % Pelvis
        'IC_L'; 'LC_L';  'LM_L'; 'TH1_L';'TH2_L';'TH3_L';'TH4_L'; 'SK1_L';'SK2_L';'SK3_L';'SK4_L';...% Left Leg
        'IC_R'; 'LC_R';  'LM_R'; 'TH1_R';'TH2_R';'TH3_R';'TH4_R'; 'SK1_R';'SK2_R';'SK3_R';'SK4_R';...% Right Leg
        'B2_1';'B2_2';'B2_3';'B2_4';'B3_1';'B3_2';'B3_3';'B3_4';'ARM2';'ARM3';... % Pads markers}
% % % % %         'HJC_R';'HJC_L';'KJC_R';'AJC_R';'KJC_L';'AJC_L'; ... % Joint Centers Markers
        'EL_L';'EL_R';'H1_L';'H2_L';'H3_L';'H1_R';'H2_R';'H3_R'; ... % Arms
        'HL_L';'HL_R';'C5M_L';'C5M_R';'CM5_1_L'; 'CM5_1_R';};  % Foot
end



data.marker_data = btk_sortc3d(data.marker_data,marker_names); % creates .trc file


if isempty(strfind(lower(fname),'cal'))
    %% DEFINE EVENTS E1 E2 E3 E4
    
    E(1) = 5; %the time of "set" call
    
    E(2) = E(1)-0.8; %Start just before set call
    
    E(3) = E(1)+1.2; % Change 2.00 if you want a longer trial
    
    %% Define start and end frames
    % Define start and end frame from the events to write the appropriate TRC
    % and MOT files for the OpenSim simulations
    
    data.Start_Frame = round(E(2)/(1/data.marker_data.Info.frequency));
    data.End_Frame = round(E(3)/(1/data.marker_data.Info.frequency));
    
    %% Define file name prefix structure
    % Added a nameing structure to avoid confusion throughout the script,
    % participant code refers to the actual participant number (eg '066')
    % where as trial code refers to the trial (eg '20_4')
    
    %** Hard coding mass and height at the moment  **
    mass=[90.5 109.1 125.5 116.5 100.8 115.5 86.7 94.1 82.3];
    
    switch data.Name(1)
        case 'A'
            data.Mass= mass(1);
            data.Height=1693;
        case 'B'
            data.Mass=mass(2);
            data.Height=1766;
        case 'C'
            data.Mass=mass(3);
            data.Height=1855;
        case 'D'
            data.Mass=mass(4);
            data.Height=1761;
        case 'E'
            data.Mass=mass(5);
            data.Height=1792;
        case 'F'
            data.Mass=mass(6);
            data.Height=1930;
        case 'G'
            data.Mass=mass(7);
            data.Height=1752;
        case 'H'
            data.Mass=mass(8);
            data.Height=1865;
        case 'I'
            data.Mass=mass(9);
            data.Height=1830;
    end
    
    participant_code = data.Name;
    
    trial_code = fname(1:end-4);
    %% Open and read Scrum Machine forces - MRF file needs to be in the same folder of .c3D files (xx_cal and trial_xx)
    if ~isempty(ls ('*MRF.mat'))
        SMfilename = [participant_code 'MRF.mat'];
        load(SMfilename)
        % Resample Scrum Machine and Pressure sensors signals from 500 hz up to
        % 2500 hz
        
        SMF=[fR.Rx(~isnan(fR.Rx(:,2)),2:3) fR.Ry(~isnan(fR.Ry(:,2)),2:3) fR.Rz(~isnan(fR.Rz(:,2)),2:3)];
        SMFraw=[zeros( length(fR.Rx(isnan(fR.Rx(:,2)),2:3)),6); SMF];
        SMFr = resample(SMFraw(:,:),2500,500);
        [SMRr_peaky,SMRr_fpeaky]=FINDMAX(fR.Ry(:,2:3));
        [SMRr_peakx,SMRr_fpeakx]=FINDMAX(fR.Rx(:,2:3));
        [SMRr_peakz,SMRr_fpeakz]=FINDMAX(fR.Rz(:,2:3));
        time=fR.t+5; % to make it consistent with OpenSim
        t_peakSMF=time(round(mean(SMRr_fpeaky(1:2))));
        t_rENG=time(round(mean(fR.f_rENG(2:3))));
        SMRF=[SMRr_peakx,SMRr_peaky,SMRr_peakz];
    else
        [fR]=INITIALISE_MRF_500Hz_CSI('*.csv',pname,abs((data.marker_data.Markers.ARM2(:,2)+data.marker_data.Markers.ARM3(:,2))/2- data.marker_data.Markers.ClavM(:,2)));
        
        % Resample Scrum Machine and Pressure sensors signals from 500 hz up to
        % 2500 hz
        
        SMF=[fR.Rx(~isnan(fR.Rx(:,2)),2:3) fR.Ry(~isnan(fR.Ry(:,2)),2:3) fR.Rz(~isnan(fR.Rz(:,2)),2:3)];
        SMFraw=[zeros( length(fR.Rx(isnan(fR.Rx(:,2)),2:3)),6); SMF];
        SMFr = resample(SMFraw(:,:),2500,500);
        [SMRr_peaky,SMRr_fpeaky]=FINDMAX(fR.Ry(:,2:3));
        [SMRr_peakx,SMRr_fpeakx]=FINDMAX(fR.Rx(:,2:3));
        [SMRr_peakz,SMRr_fpeakz]=FINDMAX(fR.Rz(:,2:3));
        time=fR.t+5; % to make it consistent with OpenSim
        t_peakSMF=time(round(mean(SMRr_fpeaky(1:2))));
        t_rENG=time(round(mean(fR.f_rENG(2:3))));
        SMRF=[SMRr_peakx,SMRr_peaky,SMRr_peakz];
    end
    
    %% Pressure
    % 0.508 cm row and col spacing
    % row 78 - cols 31
    % Save scrum machine markers before sortc3D function. These markers will be
    % used to reconstruct scrum machine force COP. This procedure has to be
    % run with both static and dynamic trials, because of different markers positions on Scrum
    % machine.
    
    load('press_B2_1_s.mat');
    load('press_B2_2_s.mat');
    load('press_B2_3_s.mat');
    load('press_B2_4_s.mat');
    
    load('press_B3_1_s.mat');
    load('press_B3_2_s.mat'); % no tracked
    load('press_B3_3_s.mat');
    load('press_B3_4_s.mat');
    
    load('DM204TekR.mat');
    load('DM204TekL.mat');
    
    % Resample Tekscan CoP 500 to 250 Hz
    DM204TekR=[resample(DM204R(:,2),250,500) resample(DM204R(:,1),250,500) zeros(2500,1)];
    DM204TekL=[resample(DM204L(:,2),250,500) resample(DM204L(:,1),250,500) zeros(2500,1)];
    % % %     % Get CoP period of interest
    % % %     DM204TekR_event=DM204TekR(data.Start_Frame:data.End_Frame,:);
    
    
    
    % Position of the markers placed on the rear part of the pad
    press_B2_1_d=data.marker_data.Markers.B2_1;
    press_B2_2_d=data.marker_data.Markers.B2_2;
    press_B2_3_d=data.marker_data.Markers.B2_3;
    press_B2_4_d=data.marker_data.Markers.B2_4;
    
    
    press_B3_1_d=data.marker_data.Markers.B3_1;
    press_B3_2_d=data.marker_data.Markers.B3_2;
    press_B3_3_d=data.marker_data.Markers.B3_3;
    press_B3_4_d=data.marker_data.Markers.B3_4;
    
    for i=1:length(press_B2_1_d)
        plane_Pad2{i} = createPlane(press_B2_1_d(i,:), press_B2_2_d(i,:), press_B2_4_d(i,:));
        plane_Pad3{i} = createPlane(press_B3_1_d(i,:), press_B3_2_d(i,:), press_B3_3_d(i,:));
        % Create reference system on pads
        angle_rot_Pad2(i,:)=createRS(press_B2_1_d(i,:),press_B2_2_d(i,:),press_B2_4_d(i,:));
        angle_rot_Pad3(i,:)=createRS(press_B3_1_d(i,:),press_B3_2_d(i,:),press_B3_3_d(i,:));
        %
    end
    % unwrapping
    angle_rot_Pad2= unwrap(angle_rot_Pad2);
    angle_rot_Pad3= unwrap(angle_rot_Pad3);
    
    for i=1:length(press_B2_1_d)
        DM204TekR_IRS(i,:)=DM204TekR(i,:)*eulerAngles2rotationMatrix(angle_rot_Pad2(i,:))'*1000;
        DM204TekL_IRS(i,:)=DM204TekL(i,:)*eulerAngles2rotationMatrix(angle_rot_Pad3(i,:))'*1000;
        
    end
    
    
    % CoP shoulder force = Pos Marker rear pad + (mediaCoP with respect to
    % Pos Marker in the back)/3
    % Tekscan)-(media
    press_B2_1(:,1)=press_B2_1_d(:,1) + (mean(press_B2_1_s(:,3))*1000-mean(press_B2_1_d(1:10,1)))+DM204TekR_IRS(1:length(press_B2_1_d),1); %z - mediolateral
    press_B2_1(:,2)=press_B2_1_d(:,2) + (mean(press_B2_1_s(:,1))*1000-mean(press_B2_1_d(1:10,2)))/3; %x - anteroposterior
    press_B2_1(:,3)=press_B2_1_d(:,3) + (mean(press_B2_1_s(:,2))*1000-mean(press_B2_1_d(1:10,3)))+DM204TekR_IRS(1:length(press_B2_1_d),3); %y - vertical
    
    
    press_B3_1(:,1)=press_B3_1_d(:,1) + (mean(press_B3_1_s(:,3))*1000-mean(press_B3_1_d(1:10,1)))+DM204TekL_IRS(1:length(press_B3_1_d),1); %z - mediolateral
    press_B3_1(:,2)=press_B3_1_d(:,2) + (mean(press_B3_1_s(:,1))*1000-mean(press_B3_1_d(1:10,2)))/3;%x - anteroposterior
    press_B3_1(:,3)=press_B3_1_d(:,3) + (mean(press_B3_1_s(:,2))*1000-mean(press_B3_1_d(1:10,3)))+DM204TekL_IRS(1:length(press_B3_1_d),3);%y - vertical
    
    
    % Still x (mediolateral), y (antero-posterior) and z (vertical)
    SM_B2_x=press_B2_1(:,1);
    SM_B2_y=press_B2_1(:,2);
    SM_B2_z=press_B2_1(:,3);
    SM_B3_x=press_B3_1(:,1);
    SM_B3_y=press_B3_1(:,2);
    SM_B3_z=press_B3_1(:,3);
    
    % Force is zero before engagement
    frametrENG=find(time==t_rENG)*5;
    SMFr(1:frametrENG,2)=zeros(frametrENG,1);
    SMFr(1:frametrENG,4)=zeros(frametrENG,1);
    SMFr(1:frametrENG,6)=zeros(frametrENG,1);
    SMFr(1:frametrENG,1)=zeros(frametrENG,1);
    SMFr(1:frametrENG,3)=zeros(frametrENG,1);
    SMFr(1:frametrENG,5)=zeros(frametrENG,1);
    
    % Rotate forces
    angle_rot_Pad2_2500=resample(angle_rot_Pad2,2500,250);
    angle_rot_Pad3_2500=resample(angle_rot_Pad3,2500,250);
    
    angle_rot_pad2_avg=mean(angle_rot_Pad2_2500(length(angle_rot_Pad2_2500)-100:length(angle_rot_Pad2_2500)-50,:));
    angle_rot_pad3_avg=mean(angle_rot_Pad2_2500(length(angle_rot_Pad3_2500)-100:length(angle_rot_Pad3_2500)-50,:));
    
    
    diff_2=25000-length(angle_rot_Pad2_2500);
    array_diff_2=[ones(diff_2,1)*angle_rot_pad2_avg(1,1)' ones(diff_2,1)*angle_rot_pad2_avg(1,2)' ones(diff_2,1)*angle_rot_pad2_avg(1,3)'];
    angle_rot_Pad2_force=[angle_rot_Pad2_2500;array_diff_2];
    
    diff_3=25000-length(angle_rot_Pad3_2500);
    array_diff_3=[ones(diff_3,1)*angle_rot_pad3_avg(1,1)' ones(diff_3,1)*angle_rot_pad3_avg(1,2)' ones(diff_3,1)*angle_rot_pad3_avg(1,3)'];
    angle_rot_Pad3_force=[angle_rot_Pad3_2500;array_diff_3];
    
    for i=1:length(angle_rot_Pad2_2500)
        
        switch data.Name(1)
            case 'H'
                SMFr_rotB3(i,:)=[SMFr(i,2) -SMFr(i,4) SMFr(i,6)];
                SMFr_rotB2(i,:)=[SMFr(i,1) -SMFr(i,3) SMFr(i,5)];
            case 'I'
                SMFr_rotB3(i,:)=[SMFr(i,2) -SMFr(i,4) SMFr(i,6)];
                SMFr_rotB2(i,:)=[SMFr(i,1) -SMFr(i,3) SMFr(i,5)];
            case 'C'
                SMFr_rotB3(i,:)=[SMFr(i,2) SMFr(i,4) SMFr(i,6)]*eulerAngles2rotationMatrix(angle_rot_Pad2_force(i,:))';
                SMFr_rotB2(i,:)=[SMFr(i,1) SMFr(i,3) SMFr(i,5)]*eulerAngles2rotationMatrix(angle_rot_Pad3_force(i,:))';
            case 'D'
                SMFr_rotB3(i,:)=[SMFr(i,2) SMFr(i,4) SMFr(i,6)]*eulerAngles2rotationMatrix(angle_rot_Pad2_force(i,:))';
                SMFr_rotB2(i,:)=[SMFr(i,1) SMFr(i,3) SMFr(i,5)]*eulerAngles2rotationMatrix(angle_rot_Pad3_force(i,:))';
            case 'I'
                SMFr_rotB3(i,:)=[SMFr(i,2) SMFr(i,4) SMFr(i,6)]*eulerAngles2rotationMatrix(angle_rot_Pad2_force(i,:))';
                SMFr_rotB2(i,:)=[SMFr(i,1) SMFr(i,3) SMFr(i,5)]*eulerAngles2rotationMatrix(angle_rot_Pad3_force(i,:))';
                
        end
        
    end
    
    data.fp_data.FP_data(3,1).Beam3=[-SMFr_rotB3(:,1) SMFr_rotB3(:,2) -SMFr_rotB3(:,3)]; % Change sign on Z due to different ref system between scrum machine forces and pad ref system
    data.fp_data.FP_data(4,1).Beam2=[-SMFr_rotB2(:,1) SMFr_rotB2(:,2) -SMFr_rotB2(:,3)]; %
    
    data.fp_data.GRF_data(3,1).F=[-SMFr_rotB3(:,1) SMFr_rotB3(:,2) -SMFr_rotB3(:,3)];
    data.fp_data.GRF_data(4,1).F=[-SMFr_rotB2(:,1) SMFr_rotB2(:,2) -SMFr_rotB2(:,3)];
    
    data.fp_data.GRF_data(3,1).M=zeros(length(SMFr),3);
    data.fp_data.GRF_data(4,1).M=zeros(length(SMFr),3);
    
    data.fp_data.GRF_data(3,1).P=[resample(SM_B3_x,2500,250) resample(SM_B3_y,2500,250) resample(SM_B3_z,2500,250)];
    data.fp_data.GRF_data(4,1).P=[resample(SM_B2_x,2500,250) resample(SM_B2_y,2500,250) resample(SM_B2_z,2500,250)];
    
    
    %% C3D TO TRC
    %'v5' of this script is modified from previous ones since person is now
    %walking in the 'forward' direction of the treadmill (towards the wall) and
    %so coordinate system transformations have changed.
    
    data = btk_c3d2trc_v5(data,'off'); % animation off
        
    
        %% INVERSE KINEMATIC ANALYSIS 
        % Define the standard OpenSim files for IK
        ModelFile = [participant_code '_SCALED_MANUAL.osim'];
        ResultsDirectory = cd;
        IKTasksFile = [osim_path 'Rugby_Model_FULL_IK_Tasks.xml'];
    
        % Set up the XML file
        setup_InverseKinematics('data',data,'ModelFile',ModelFile,...
            'IKTasksFile',IKTasksFile,'ResultsDirectory',...
            ResultsDirectory,'Accuracy',0.00005);
    
        % Call the IK tool
        com = ['ik -S ' participant_code '_Setup_InverseKinematics.xml'];
        system(com)
    
        % Load the data from MOT file and tie kinematic data to data structure
        D = load_sto_file([trial_code '_ik.mot']);   % Use D in case you want to plot or do any processing on joint angles 
    
        clear a
        
    %% BODY KINEMATIC ANALYSIS
    % SetUp BodyKinematics
    % Define the standard OpenSim files for BodyKinematics
    
    % DEFINE STANDARD OPENSIM FILES AND PATHS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CoordinatesFile =[pname trial_code '_ik.mot'];
    MOTFile = [pname trial_code '_ik.mot'];
    % % % % ForceSetFile        = [osim_path model 'GT_CMC_Actuators.xml'];
    % % % % ExternalLoadsFile   = GRFFileXML;
    MADirectoryName     = [pname 'BK_Results'];
    ModelFile = [pname participant_code '_SCALED.osim'];
    
    % Set up the XML file
    setup_kinematicsDG( 'data',data,...
        'ModelFile',ModelFile,...
        'CoordinatesFile',CoordinatesFile,...
        'MOTFILE',MOTFile,...
        'AnalysisOn', true,...
        'UseFLV', false,...
        'StartTime',data.Start_Frame,...
        'SolveForEquilibrium',false,...
        'UseModelForceSet',false,...
        'EndTime',data.End_Frame,...
        'LowPassFilterFreq', 6,...
        'DirectoryName', MADirectoryName);
    
    % Call the Body Kinematics Tool
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    com=['analyze -S ', participant_code '_Setup_BodyKinematics.xml'];
    system(com);

    % Load the data from Body Kinematics MOT file and tie body kinematic data to data structure for
    % later use
    B = load_sto_file([pname 'BK_Results\' participant_code '_' trial_code '_BK_BodyKinematics_pos_global.sto']);
    data.BodyPos = B;
    clear B
    
    % % % % %% Change CoP Scrum Forces
    Scap_CoP_r=[resample(data.BodyPos.rscapula_X,2500,250) resample(data.BodyPos.rscapula_Y,2500,250) resample(data.BodyPos.rscapula_Z,2500,250)];
    Scap_CoP_l=[resample(data.BodyPos.lscapula_X,2500,250) resample(data.BodyPos.lscapula_Y,2500,250) resample(data.BodyPos.lscapula_Z,2500,250)];
    
    Clav_CoP_r=[resample(data.BodyPos.rclavicle_X,2500,250) resample(data.BodyPos.rclavicle_Y,2500,250) resample(data.BodyPos.rclavicle_Z,2500,250)];
    Clav_CoP_l=[resample(data.BodyPos.lclavicle_X,2500,250) resample(data.BodyPos.lclavicle_Y,2500,250) resample(data.BodyPos.lclavicle_Z,2500,250)];
    
    Mid_point_r=(Scap_CoP_r+Clav_CoP_r)/2;
    Mid_point_l=(Scap_CoP_l+Clav_CoP_l)/2;
    
    
    radius=0.03;
    
    % create Sphere from 4 points R
    for i=1:length(Scap_CoP_r)
        
        p1= Mid_point_r(i,:);
        p2= [Mid_point_r(i,1)+radius Mid_point_r(i,2) Mid_point_r(i,3)];
        p3= [Mid_point_r(i,1) Mid_point_r(i,2)+radius Mid_point_r(i,3)];
        p4= [Mid_point_r(i,1) Mid_point_r(i,2) Mid_point_r(i,3)+radius];
        
        sphereR(i,:) = createSphere(p1,p2,p3,p4);
        
        CoP_r(i,:)=data.fp_data.GRF_data(3,1).P(data.Start_Frame*10+i-1,:);
        
        distR(i) = distancePoints3d(p1,CoP_r(i,:));
        
        if distR>radius
            
            LineSpeherePoints= intersectLineSphere([p1 (CoP_r(i,:)-p1)/norm(p1-CoP_r(i,:))], [p1 radius]);
            if LineSpeherePoints(1,2)>LineSpeherePoints(2,2)
                CoP_newR(i,:) =LineSpeherePoints(1,:);
            else
                CoP_newR(i,:) =LineSpeherePoints(2,:);
            end
            
        else
            CoP_newR(i,:)=data.fp_data.GRF_data(3,1).P(data.Start_Frame*10+i-1,:);
        end
        
    end
    % create Sphere from 4 points L
    for i=1:length(Scap_CoP_l)
        
        p1= Mid_point_l(i,:);
        p2= [Mid_point_l(i,1)+radius Mid_point_l(i,2) Mid_point_l(i,3)];
        p3= [Mid_point_l(i,1) Mid_point_l(i,2)+radius Mid_point_l(i,3)];
        p4= [Mid_point_l(i,1) Mid_point_l(i,2) Mid_point_l(i,3)+radius];
        sphereL(i,:) = createSphere(p1,p2,p3,p4);
        
        CoP_l(i,:)=data.fp_data.GRF_data(4,1).P(data.Start_Frame*10+i-1,:);
        
        distL(i) = distancePoints3d(p1,CoP_l(i,:));
        
        if distL>radius
            
            LineSpeherePoints= intersectLineSphere([p1 (CoP_l(i,:)-p1)/norm(p1-CoP_l(i,:))], [p1 radius]);
            
            if LineSpeherePoints(1,2)>LineSpeherePoints(2,2)
                
                CoP_newL(i,:) =LineSpeherePoints(1,:);
            else
                CoP_newL(i,:) =LineSpeherePoints(2,:);
            end
            
            
            
        else
            CoP_newL(i,:)=data.fp_data.GRF_data(4,1).P(data.Start_Frame*10+i-1,:);
        end
        
    end
    
    figure;
    plot(distL);hold on
    plot(distR,'-.o');
    
    FPA_r=data.fp_data.GRF_data(3,1).P;
    FPA_l=data.fp_data.GRF_data(4,1).P;
    
    FPA_r(data.Start_Frame*10:data.End_Frame*10,:)=CoP_newR(1:length(data.fp_data.GRF_data(3,1).P(data.Start_Frame*10:data.End_Frame*10,:)),:);
    FPA_l(data.Start_Frame*10:data.End_Frame*10,:)=CoP_newL(1:length(data.fp_data.GRF_data(4,1).P(data.Start_Frame*10:data.End_Frame*10,:)),:);
    
    figure;
    plot(FPA_l);hold on
    plot(data.fp_data.GRF_data(4,1).P,'-.');
    
    figure;
    plot(FPA_r);hold on
    plot(data.fp_data.GRF_data(3,1).P,'-.');
    
    data.fp_data.GRF_data(3,1).P(data.Start_Frame*10:data.End_Frame*10,:)=CoP_newR(1:length(data.fp_data.GRF_data(3,1).P(data.Start_Frame*10:data.End_Frame*10,:)),:);
    data.fp_data.GRF_data(4,1).P(data.Start_Frame*10:data.End_Frame*10,:)=CoP_newL(1:length(data.fp_data.GRF_data(4,1).P(data.Start_Frame*10:data.End_Frame*10,:)),:);
    
    % % data.fp_data.GRF_data(3,1).P(data.Start_Frame*10:data.End_Frame*10,:)=Scap_CoP_r(1:length(data.fp_data.GRF_data(3,1).P(data.Start_Frame*10:data.End_Frame*10,:)),:);
    % % data.fp_data.GRF_data(4,1).P(data.Start_Frame*10:data.End_Frame*10,:)=Scap_CoP_l(1:length(data.fp_data.GRF_data(4,1).P(data.Start_Frame*10:data.End_Frame*10,:)),:);
    
    
    %% rewrite grf
    data = btk_c3d2trc_DC(data,'off'); % animation off
    
    
    %% INVERSE DYNAMIC ANALYSIS
    
    % % Define the standard OpenSim files for ID
    MOTFile = [pname trial_code '_ik.mot'];
    GRFFile = [pname trial_code '_grf.mot'];
    GRFFileXML = [pname trial_code '_grf.xml'];
    ModelFile = [pname participant_code '_SCALED_MANUAL.osim'];
    OutputFile = [trial_code '_ID.sto'];
    ResultsDirectory = [pname 'ID\'];
    
    
    % Write GRF XML file
    grf2xml(data,'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2','Scrum_Machine_1','Scrum_Machine_2'},...
        'AppliedToBodies',{'toes_r','toes_l','rscapula','lscapula'},...
        'ForceExpressedinBody',{'ground','ground','ground','ground'},...
        'PointExpressedinBody',{'ground','ground','ground','ground'},...
        'GRFFile',GRFFile,'MOTFile',...
        MOTFile,'LowPassFilterForKinematics',6,'OutputFile',GRFFileXML);
    
    % Set up the XML file
    setup_InverseDynamics('data',data, 'ModelFile',ModelFile,'MOTFile',MOTFile,...
        'GRFFile',GRFFileXML,'LowPassFilterFreq',6,'OutputFile', OutputFile, 'ResultsDirectory',  pname, 'InputDirectory', ResultsDirectory);
    
    % Call the ID tool
    com = ['id -S ' participant_code '_Setup_InverseDynamics.xml'];
    
    system(com);
    
    clear D
    
    % Load the data from ID MOT file and tie kinetic data to data structure
    % if you want to perform any analysis or plot data;
    D = load_sto_file([pname trial_code '_ID.sto']);
    data.Kinetics = D;

    clear D
    %
    
    
    %% RESIDUAL REDUCTION ANALYSIS
    
    
    % Set up RRA - define the standard OpenSim files for RRA
    %(MOT and GRF_MOT defined above along with Model file)
    
    MOTFile = [pname trial_code '_ik.mot'];
    GRFFile = [pname trial_code '_grf.mot'];
    GRFFileXML = [pname trial_code '_grf.xml'];
    
    % %  Needed to change 'osim_path' to  '../../../RugbyModel/' to make it
    % working. Don't know the reason why it is needed.
    
    RRATaskFile = [' ../../../RugbyModel/'  model '_FULL_RRA_Tasks.xml'];%[data.Name '_RRA_Tasks.xml'];
    RRAForceFile = [' ../../../RugbyModel/' model '_FULL_RRA_Actuators.xml'];%[data.Name '_RRA_Actuators.xml'];
    %RRAConstraintsFile = [osim_path model '_RRA_ControlConstraints.xml'];%[data.Name '_RRA_ControlConstraints.xml'];
    RRADirectoryName = [pname 'RRAResults\'];
    ModelFile = [pname participant_code '_SCALED_MANUAL.osim'];
    
    % Write GRF XML file
    grf2xml(data,'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2','Scrum_Machine_1','Scrum_Machine_2'},...
        'AppliedToBodies',{'toes_r','toes_l','rscapula','lscapula'},...
        'ForceExpressedinBody',{'ground','ground','ground','ground'},...
        'PointExpressedinBody',{'ground','ground','ground','ground'},...
        'GRFFile',GRFFile,'MOTFile',...
        MOTFile,'LowPassFilterForKinematics',6,'OutputFile',GRFFileXML);
    
    
    % Set up the XML file
    RRA = setup_ReduceResiduals('data',data,...
        'ModelFile',ModelFile,...
        'MOTFile',MOTFile,...
        'GRFFile',GRFFileXML,...
        'RRATaskFile',RRATaskFile,...
        'RRAForceFile',RRAForceFile,...
        'AdjCOMRes','true',...
        'LowPassFilterFreq',6,...
        'DirectoryName',RRADirectoryName,...
        'InitialTime',data.time(1),...
        'FinalTime',data.time(end));
    %'RRAConstraintsFile',RRAConstraintsFile,...%'PointKinematics',PointKinematics,...
    % CALL THE RRA TOOL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,log_mes]=dos(['RRA -S ', participant_code '_setup_ReduceResiduals.xml'],'-echo');
    
    % SAVE THE WORKSPACE AND PRINT A LOG FILE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen([RRADirectoryName,data.Name,'.log'],'w+');
    fprintf(fid,'%s\n', log_mes);
    fclose(fid);
    
    
    % Call the RRA tool
    clear com
    com = ['RRA -S ' participant_code '_setup_ReduceResiduals.xml'];
    system(com);
    
    clear D
    % load the data from RRA
    D = load_sto_file([RRADirectoryName  participant_code '_' trial_code '_RRA_Actuation_force.sto']);
    data.RRA.Kinetics = D;
    
    clear D
    % load the data from RRA
    D = load_sto_file([RRADirectoryName  participant_code '_' trial_code '_RRA_states.sto']);
    data.RRA.States = D;
    %     save([(subdir_current) '.mat']);
    clear D
    % load the data from RRA
    D = load_sto_file([RRADirectoryName  participant_code '_' trial_code '_RRA_pErr.sto']);
    data.RRA.pErr = D;
    
    %% Make plots for quick data check of residuals and pErr
    figure
    
    subplot(2,1,1), plot(data.Kinetics.time,[data.Kinetics.pelvis_tx_force data.Kinetics.pelvis_ty_force data.Kinetics.pelvis_tz_force]); hold on
    subplot(2,1,1), plot(data.RRA.Kinetics.time,[data.RRA.Kinetics.FX data.RRA.Kinetics.FY data.RRA.Kinetics.FZ],'LineWidth',2);
    title('Residuals - PELVIS')
    xlabel('time [s]')
    ylabel('Force [N]')
    legend('ID Fx','ID Fy','ID Fz','RRA Fx','RRA Fy','RRA Fz'); hold off
    subplot(2,1,2), plot(data.Kinetics.time,[data.Kinetics.gndroll_moment data.Kinetics.gndpitch_moment data.Kinetics.gndyaw_moment]); hold on
    subplot(2,1,2), plot(data.RRA.Kinetics.time,[-data.RRA.Kinetics.MY data.RRA.Kinetics.MZ data.RRA.Kinetics.MX],'LineWidth',2);
    title('Residuals - SPINE')
    xlabel('time [s]')
    ylabel('Moment [Nm]')
    legend('ID roll','ID pitch','ID yaw','RRA MY','RRA MZ','RRA MX'); hold off
    
    figure
    
    subplot(2,1,1), plot(data.RRA.pErr.time,[data.RRA.pErr.pelvis_tx data.RRA.pErr.pelvis_ty data.RRA.pErr.pelvis_tz]); hold on
    
    title('pErr Pelvis Trasl')
    xlabel('time [s]')
    ylabel('pErr [m]')
    legend('err_tx','err_ty','err_tz'); hold off
    subplot(2,1,2), plot(data.RRA.pErr.time,[data.RRA.pErr.spine_list*(180/pi) data.RRA.pErr.spine_rotation*(180/pi) data.RRA.pErr.spine_tilt*(180/pi)]); hold on
    subplot(2,1,2), plot(data.RRA.pErr.time,[data.RRA.pErr.gndroll*(180/pi) data.RRA.pErr.gndpitch*(180/pi) data.RRA.pErr.gndyaw*(180/pi)],'LineWidth',2);
    title('pErr - SPINE and PELVIS Angular Kinematics')
    xlabel('time [s]')
    ylabel('pErr [degree]')
    legend('Spine List','Spine rotation','Spine tilt','Pelvis roll','Pelvis Pitch','Pelvis Yaw'); hold off
    
    figure
    
    subplot(3,2,1), plot(data.RRA.pErr.time,[data.RRA.pErr.hip_flexion_r*(180/pi) data.RRA.pErr.hip_adduction_r*(180/pi) data.RRA.pErr.hip_rotation_r*(180/pi)]); hold on
    subplot(3,2,1), plot(data.RRA.pErr.time,[data.RRA.pErr.hip_flexion_l*(180/pi) data.RRA.pErr.hip_adduction_l*(180/pi) data.RRA.pErr.hip_rotation_l*(180/pi)],'LineWidth',2); hold on
    title('pErr HIP')
    xlabel('time [s]')
    ylabel('pErr [degree]')
    legend('err hip flex R','err hip adduction R','err hip rotation R','err hip flex L','err hip adduction L','err hip rotation L'); hold off
    subplot(3,2,2), plot(data.RRA.pErr.time,[data.RRA.pErr.knee_angle_r*(180/pi) data.RRA.pErr.knee_angle_l*(180/pi) data.RRA.pErr.ankle_angle_r*(180/pi) data.RRA.pErr.ankle_angle_l*(180/pi)],'LineWidth',2); hold on
    title('pErr - KNEE and ANKLE')
    xlabel('time [s]')
    ylabel('pErr [degree]')
    legend('Knee R','Knee L','Ankle R','Ankle L'); hold off
    subplot(3,2,3), plot(data.RRA.pErr.time,[data.RRA.pErr.pitch1*(180/pi) data.RRA.pErr.roll1*(180/pi) data.RRA.pErr.yaw1*(180/pi)],'LineWidth',2); hold on
    title('pErr Head')
    xlabel('time [s]')
    ylabel('pErr [degree]')
    legend('Head flex/ext','Head Lat Bend','Head Axial rot'); hold off
    subplot(3,2,4), plot(data.RRA.pErr.time,[data.RRA.pErr.arm_flex_r*(180/pi) data.RRA.pErr.arm_flex_l*(180/pi) data.RRA.pErr.arm_add_r*(180/pi) data.RRA.pErr.arm_add_l*(180/pi) data.RRA.pErr.arm_rot_r*(180/pi) data.RRA.pErr.arm_rot_l*(180/pi)],'LineWidth',2); hold on
    title('pErr - ARM')
    xlabel('time [s]')
    ylabel('pErr [degree]')
    legend('Arm flex R','Arm flex L','Arm add R','Arm add L','Arm rot R','Arm rot L'); hold off
    subplot(3,2,5), plot(data.RRA.pErr.time,[data.RRA.pErr.spine_list*(180/pi) data.RRA.pErr.spine_rotation*(180/pi) data.RRA.pErr.spine_tilt*(180/pi)],'LineWidth',2); hold on
    title('pErr - SPINE')
    xlabel('time [s]')
    ylabel('pErr [degree]')
    legend('Spine List','Spine Rotation','Spine tilt'); hold off
    subplot(3,2,6), plot(data.RRA.pErr.time,[data.RRA.pErr.pitch2*(180/pi) data.RRA.pErr.roll2*(180/pi) data.RRA.pErr.yaw2*(180/pi)],'LineWidth',2); hold on
    title('pErr - NECK')
    xlabel('time [s]')
    ylabel('pErr [degree]')
    legend('Neck flex/ext','Neck Lat Bend','Neck Axial rot'); hold off
    
    figure
    
    subplot(3,2,1), plot(data.Kinetics.time,[data.Kinetics.pitch2_moment data.Kinetics.roll2_moment data.Kinetics.yaw2_moment],'--'); hold on
    subplot(3,2,1), plot(data.RRA.Kinetics.time,[data.RRA.Kinetics.pitch2 data.RRA.Kinetics.roll2 data.RRA.Kinetics.yaw2],'LineWidth',2);
    title('Residuals - NECK')
    xlabel('time [s]')
    ylabel('Moments [Nm]')
    legend('ID Flex-Ext','ID Lat Bend','ID Axial Rot','RRA FLex-Ext','RRA Lat Bend','RRA Axial Rot'); hold off
    subplot(3,2,2), plot(data.Kinetics.time,[data.Kinetics.pitch1_moment data.Kinetics.roll1_moment data.Kinetics.yaw1_moment],'--'); hold on
    subplot(3,2,2), plot(data.RRA.Kinetics.time,[data.RRA.Kinetics.pitch1 data.RRA.Kinetics.roll1 data.RRA.Kinetics.yaw1],'LineWidth',2);
    title('Residuals - HEAD')
    xlabel('time [s]')
    ylabel('Moment [Nm]')
    legend('ID Flex-Ext','ID Lat Bend','ID Axial Rot','RRA FLex-Ext','RRA LAt Bend','RRA Axial Rot'); hold off
    subplot(3,2,3), plot(data.Kinetics.time,[data.Kinetics.hip_flexion_r_moment data.Kinetics.hip_adduction_r_moment data.Kinetics.hip_rotation_r_moment data.Kinetics.hip_flexion_l_moment data.Kinetics.hip_adduction_l_moment data.Kinetics.hip_rotation_l_moment],'--'); hold on
    subplot(3,2,3), plot(data.RRA.Kinetics.time,[data.RRA.Kinetics.hip_flexion_r data.RRA.Kinetics.hip_adduction_r data.RRA.Kinetics.hip_rotation_r data.RRA.Kinetics.hip_flexion_l data.RRA.Kinetics.hip_adduction_l data.RRA.Kinetics.hip_rotation_l],'LineWidth',2);
    title('Residuals - HIP')
    xlabel('time [s]')
    ylabel('Moment [Nm]')
    legend('ID Hip R Flex','ID Hip R Add','ID Hip R Rot','ID Hip L Flex','ID Hip L Add','ID Hip L Rot','RRA Hip R Flex','RRA Hip R Add','RRA Hip R Rot','RRA Hip L Flex','RRA Hip L Add','RRA Hip L Rot'); hold off
    subplot(3,2,4), plot(data.Kinetics.time,[data.Kinetics.knee_angle_r_moment data.Kinetics.ankle_angle_r_moment data.Kinetics.knee_angle_l_moment data.Kinetics.ankle_angle_l_moment],'--'); hold on
    subplot(3,2,4), plot(data.RRA.Kinetics.time,[data.RRA.Kinetics.knee_angle_r data.RRA.Kinetics.ankle_angle_r data.RRA.Kinetics.knee_angle_l data.RRA.Kinetics.ankle_angle_l],'LineWidth',2);
    title('Residuals - KNEE and ANKLE')
    xlabel('time [s]')
    ylabel('Moment [Nm]')
    legend('ID Knee R Flex','ID Ankle R Plantaflex','ID Knee L Flex','ID Ankle L Plantaflex','RRA Knee R Flex','RRA Ankle R Plantaflex','RRA Knee L Flex','RRA Ankle L Plantaflex'); hold off
    subplot(3,2,5), plot(data.Kinetics.time,[data.Kinetics.arm_flex_r_moment data.Kinetics.arm_add_r_moment data.Kinetics.arm_rot_r_moment data.Kinetics.arm_flex_l_moment data.Kinetics.arm_add_l_moment data.Kinetics.arm_rot_l_moment],'--'); hold on
    subplot(3,2,5), plot(data.RRA.Kinetics.time,[data.RRA.Kinetics.arm_flex_r data.RRA.Kinetics.arm_add_r data.RRA.Kinetics.arm_rot_r data.RRA.Kinetics.arm_flex_l data.RRA.Kinetics.arm_add_l data.RRA.Kinetics.arm_rot_l],'LineWidth',2);
    title('Residuals - ARMS')
    xlabel('time [s]')
    ylabel('Moment [Nm]')
    legend('ID Arm R Flex','ID Arm R Add','ID Arm R Rot','ID Arm L Flex','ID Arm L Add','ID Arm L Rot','RRA Arm R Flex','RRA Arm R Add','RRA Arm R Rot','RRA Arm L Flex','RRA Arm L Add','RRA Arm L Rot'); hold off
    subplot(3,2,6), plot(data.Kinetics.time,[data.Kinetics.spine_list_moment data.Kinetics.spine_rotation_moment data.Kinetics.spine_tilt_moment],'--'); hold on
    subplot(3,2,6), plot(data.RRA.Kinetics.time,[data.RRA.Kinetics.spine_list data.RRA.Kinetics.spine_rotation data.RRA.Kinetics.spine_tilt],'LineWidth',2);
    title('Residuals - SPINE')
    xlabel('time [s]')
    ylabel('Moment [Nm]')
    legend('ID spine list','ID spine rot','ID spine tilt','RRA spine list','RRA spine rot','RRA spine tilt'); hold off
    
    
    %% Define file name prefix structure
    % Added a nameing structure to avoid confusion throughout the script,
    % particpant code refers to the actual participant number (eg '066')
    % where as trial code refers to the trial (eg '20_4')
    
    participant_code = data.Name;
    
    trial_code = fname(1:end-4);
    
    
    %% ADJUST MODEL PROPERTIES ROUTINE
    
    % Using new API functions
    InputModel = [participant_code '_RRA_adjusted.osim'];
    OutputModel = [participant_code '_RRA_ADJ_Mass_.osim'];
    
    adjustmodelmass(2,InputModel,OutputModel,log_mes);
    
    
    
else % Scaling process for the static trial
    
    data = btk_c3d2trc_v5(data, 'off');
    
    %  %** Hard coding mass and height at the moment AND again below in dynamic trials **
    data.Mass=116.5;
    data.Height=1761;
    
    press_B2_1_s=[data.marker_data.Markers.B2_1(:,1) data.marker_data.Markers.B2_1(:,2) data.marker_data.Markers.B2_1(:,3)];
    press_B2_2_s=[data.marker_data.Markers.B2_2(:,1) data.marker_data.Markers.B2_2(:,2) data.marker_data.Markers.B2_2(:,3)];
    press_B2_3_s=[data.marker_data.Markers.B2_3(:,1) data.marker_data.Markers.B2_3(:,2) data.marker_data.Markers.B2_3(:,3)];
    press_B2_4_s=[data.marker_data.Markers.B2_4(:,1) data.marker_data.Markers.B2_4(:,2) data.marker_data.Markers.B2_4(:,3)];
    
    press_B3_1_s=[data.marker_data.Markers.B3_1(:,1) data.marker_data.Markers.B3_1(:,2) data.marker_data.Markers.B3_1(:,3)];
    press_B3_2_s=[data.marker_data.Markers.B3_2(:,1) data.marker_data.Markers.B3_2(:,2) data.marker_data.Markers.B3_2(:,3)];
    press_B3_3_s=[data.marker_data.Markers.B3_3(:,1) data.marker_data.Markers.B3_3(:,2) data.marker_data.Markers.B3_3(:,3)];
    press_B3_4_s=[data.marker_data.Markers.B3_4(:,1) data.marker_data.Markers.B3_4(:,2) data.marker_data.Markers.B3_4(:,3)];
    
    save('press_B2_1_s','press_B2_1_s');
    save('press_B2_2_s','press_B2_2_s');
    save('press_B2_3_s','press_B2_3_s');
    save('press_B2_4_s','press_B2_4_s');
    
    save('press_B3_1_s','press_B3_1_s');
    save('press_B3_2_s','press_B3_2_s');
    save('press_B3_3_s','press_B3_3_s');
    save('press_B3_4_s','press_B3_4_s');
    
    % Define the starndard OpenSim files for scaling
    ModelFile = [osim_path model '.osim'];
    MeasurementSetFile = [osim_path model '_FULL_Scale_MeasurementSet.xml'];
    ScaleTasksFile = [osim_path model '_FULL_Scale_Tasks.xml'];
    
    if isfield(data.marker_data.Markers,'C1M_R'),
        loc_r=[(data.marker_data.Markers.C5M_R(1,:)+data.marker_data.Markers.C1M_R(1,:))/2];
        loc_l=[(data.marker_data.Markers.C5M_L(1,:)+data.marker_data.Markers.C1M_L(1,:))/2];
    else
        loc_r=data.marker_data.Markers.C5M_R(1,:);
        loc_l=data.marker_data.Markers.C5M_L(1,:);
    end
    save ('loc_r.mat','loc_r');
    save ('loc_l.mat','loc_l');
    
    % % % %     RugbyModel=Model(ModelFile);
    
    
    setup_scale('data',data,'ModelFile',ModelFile,'ScaleTasksFile',ScaleTasksFile,...
        'MeasurementSetFile',MeasurementSetFile, 'PreserveMass','true','InitialTime',data.time(1),...
        'FinalTime',data.time(10));
    
    com = ['scale -S ' data.Name '_Setup_Scale.xml'];
    system(com);
    
    %save('Scale.mat');
    
end



