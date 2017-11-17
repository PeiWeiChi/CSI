function [f_ONP] = ONSET_TIME_REC_vCSI(Praw,sf,f_ref,type)

% EXAMPLE.m shows an example of how a well arranged matlab file is
% designed. In this "head" a brief description of what the function does
% and how it works should be placed.
%
%
% Input:
%       - varin1       Praw=signal;
%       - varin2       f_p=samp freq (Hz);
%       - varin3      type=['D','K'];
%       - varin4      f_ref= frame of t_rENG from Xsens or KTR;
% Output:
%       - varout1      put the characteristics of output variable #1 and
%                      the characteristics it must have
%       - varout2      put the characteristics of output variable #1 and
%                      the characteristics it must have
%
%
% Recalled functions or nested functions ( >> put all the recalled functions)
%       - ADF_BUTTERWORTH.m
%
%       - FINDMIN.m
%           + RECALLED_21.m
%               + RECALLED_211.m
%               + RECALLED_212.m
%           + RECALLED_22.m
%       - RECALLED_3.m
%
%
% Useful References:
% 1) Soda P, Mazzoleni S, Cavallo G, Guglielmelli E, Iannello G. (2010) - Human movement onset detection from isometric force and torque measurements: a supervised pattern recognition approach.
% Artif Intell Med. 2010 Sep;50(1):55-61. doi: 10.1016/j.artmed.2010.04.008. Epub 2010 May 26.
% 2)Botzer L, Karniel A.(2009 Jun 29;3:61. doi: 10.3389/neuro.20.002.2009.)Front Neurosci.SourceDepartment of Biomedical Engineering, Ben-Gurion University of the Negev Israel.
%  - A Simple and Accurate onset Detection Method for a Measured Bell-shaped Speed Profile.
%
%
%
% Written by:
% -Dario Cazzola-   University of Bath, Jan 2013
%
% Modified by:
% -Ezio Preatoni-   University of Bath, Jan 2013

%% input check
%check arguments number and content
error(nargchk(4,4,nargin,'struct'))
if (isempty(Praw) || isempty(sf) || isempty(type))
    errordlg('Check input: some input parameters cannot be an empty matrix!');
    return
end

%check data dimensions
if (~isvector(Praw)) || (~isnumeric(sf)) || (~ischar(type))
    errordlg('Check input dimensions and/or type!');
    whos('Praw','sf','type')
    return
end

%make sure input are column arrays
Praw= Praw(:);

%% identification of area of search
%definition of time vector
t= ((1:length(Praw))'-1)/sf;

if isempty(f_ref) || isnan(f_ref) || ~isnumeric(f_ref)
    f_ref= 1;
end

%definition of search area
switch type
    case 'D1'
        df_start= min(150,f_ref);
        df_end= min(50,length(Praw)-f_ref);
    case 'D2'
        df_start= 1;
        df_end= min(200,length(Praw)-f_ref);
    case 'D3'
        df_start= 1;
        df_end= min(750,length(Praw)-f_ref);
    case 'K'
        df_start= 1;
        df_end= round(length(Praw)/2)-f_ref;
    otherwise
        e= errordlg('BAD INPUT! Please check inputs: available options are D1, D2, D3 or K',...
            'Input Error');
        uiwait(e); return
end
Pcut=Praw((f_ref-df_start+1):(f_ref+df_end));
tcut= t((f_ref-df_start+1):(f_ref+df_end));

%% identification of the time of onset
%chooses different methods depending on type of input
if strcmpi(type(1),'D'),
    %raw signal filtering
    [Pcutf,~]= ADF_BUTTERWORTH(Pcut,sf,[10 20],1,'N','N');
    %Second Derivative Method
    [~,maxIndex_Pcutf]= max(Pcutf(df_start:end));
    D1= gradient(Pcutf(1:maxIndex_Pcutf+df_start-1));
    bline= mean(D1(1:round(length(D1)*0.1)));                   %baseline (initial 10% of the points in D1)
    [maxVal,~]= max(D1(:));
    k3=(mean(D1(round(length(D1)*0.5):length(D1)))/max(D1)); % % % % %
    % % % % %     k3=(mean(D1(20:length(D1)))/max(D1))*(mean(gradient(D1(20:length(D1))))/max(gradient(D1))); % % % % %
    % % % % %     k3= 0.15; % % % % %
    Th= (maxVal-bline)*k3;
    D1_Th= nan(size(D1));
    D1_Th(D1<=Th)= D1(D1<=Th);
    D2= gradient(D1_Th);
    [~,maxIndexSDM_D2]= max(D2(:));
    f_ENG_2DM= maxIndexSDM_D2+(f_ref-df_start);    
    %Minimum Density Point (Kernel density function)
    ks_d= nan(length(Pcut),1);
    for i= 1:length(Pcut),
        ks_d(i,1)= (1/sqrt(2*pi*var(Pcut)^2))*exp(-(Pcutf(i)^2/(2*var(Pcut)^2)));
    end
    [~,MDP_t]= FINDMIN(ks_d(1:length(ks_d)));             %Change array chunck selection
    f_ENG_MDP= MDP_t+(f_ref-df_start);
    %MACC-based Onset Detection Algorithm
    U=gradient(gradient(gradient(Pcut)));                 %Jerk calculation
    [~,Umin_f]= FINDMIN(U);                               %minimum Jerk
    f_ENG_Umin= Umin_f+(f_ref-df_start);
elseif type=='K',
    %raw signal filtering
    [Pcutf,~]= ADF_BUTTERWORTH(Pcut,sf,[10 20],1,'N','N');
    [Praw3,~]= ADF_BUTTERWORTH(Pcut,sf,[3 10],1,'N','N');
    %Second Derivative Method (3 Hz)
    D1= gradient(Praw3);
    bline= 0;                   %baseline for kinematic data
    [maxVal,~]= max(D1(:));
    k3=0.10;
    Th=(maxVal-bline)*k3;
    D1_ThK=nan(length(D1),1);
    for i=1:length(D1),
        if D1(i)<Th,
            D1_ThK(i)=D1(i);
        else
            break;
        end
    end
    D2=gradient(D1_ThK);
    [~,maxIndexSDM_D2]=max(D2(:));
    f_ENG_2DM=maxIndexSDM_D2+(f_ref-df_start);
    
    %Minimum Density Point (Kernel density function)
    ks_d= nan(round((length(Pcut)*2/3)),1);
    for i= 1:round((length(Pcut)*2/3)),
        ks_d(i,1)= (1/sqrt(2*pi*var(Pcut(1:round(length(Pcut)*2/3),:))^2))*exp(-(Pcutf(i)^2/(2*var(Pcut(1:round(length(Pcut)*2/3),:))^2)));
    end
    [~,MDP_t]= FINDMIN(ks_d(1:length(ks_d))); % Change array chunck selection
    f_ENG_MDP= MDP_t+(f_ref-df_start);
    %MACC-based Onset Detection Algorithm
    U= gradient(gradient(gradient(Pcut(1:length(Pcut)*2/3,:))));        %Jerk calculation
    [~,Umin_t]= FINDMIN(U);                                             %minimum Jerk
    f_ENG_Umin= Umin_t+(f_ref-df_start);
else
    e= errordlg('Insert D for dynamic data and K for kynematic data',...
        'Input error');
    uiwait(e);
end

%% OUTPUT
f_ONP.SDM= f_ENG_2DM;
f_ONP.MDP= f_ENG_MDP;
f_ONP.Umin= f_ENG_Umin;

%% GRAPHICS
% % % % % if type=='K',   
% % % % % else
% % % % %     figure('units','normalized','position',[0.1 0.1 0.8 0.8]); hold on;
% % % % %     tcutD1= tcut(1:length(D1));
% % % % %     plot(tcutD1,D1,'c-.','linewidth',2);
% % % % %     plot(tcutD1,D1_Th,'b.');
% % % % %     plot(tcutD1,D2,'m-.');
% % % % %     hline(bline,'y');
% % % % %     hline(Th,'b');
% % % % %     plot(t,Praw,'k.-');
% % % % %     plot(tcut,Pcutf,'b--','linewidth',2);
% % % % %     plot(tcut,Pcut,'ko','markerfacecolor','w');
% % % % %     vline(t(f_ref),'k--');
% % % % %     plot(t(f_ENG_2DM),Praw(f_ENG_2DM),'ro','linewidth',2,'markerfacecolor','r');
% % % % %     plot(t(f_ENG_MDP),Praw(f_ENG_MDP),'co','linewidth',2,'markerfacecolor','c');
% % % % %     plot(t(f_ENG_Umin),Praw(f_ENG_Umin),'go','linewidth',2,'markerfacecolor','g');
% % % % %     hline(0,'k');vline(0,'k');box on;
% % % % %     xlabel('t [s]'); ylabel('input signal [au]')    
% % % % % end

