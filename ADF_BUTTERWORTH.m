function [Yf,f_BW]=ADF_BUTTERWORTH(Y,fs,cutoff,niter,pos,graph)

% [varargout]=ADF_BUTTERWORTH(varargin)
%
% ADF_BUTTERWORTH.m ... (to be completed).
%
% Input:
%       - varin1       put the characteristics of input variable #1 and
%                      the characteristics it must have
%       - varin2       put the characteristics of input variable #2 and
%                      the characteristics it must have
%
% Output:
%       - varout1      put the characteristics of output variable #1 and
%                      the characteristics it must have
%       - varout2      put the characteristics of output variable #1 and
%                      the characteristics it must have
%
%
% Recalled functions or nested functions ( >> put all the recalled functions)
%       - hline.m
%       - vline.m
%       - fCONTENT.m
%           + hline.m
%           + vline.m
%       - ADAPTIVE_BUTTER.m
%
%
% Useful References:
% 1) Erer (2007) Adaptive usage of the Butterworth digital filter. Journal
%           of Biomechanics. 40,2934-2943.
% 2) Robertson & Dowling (2003) Design and responses of Butterworth and critically 
%           damped digital filters. Journal of Electromyography and
%           Kinesiology. 13, 569-573
%
%
% Notes/Observation:
% -->   make the initial cutoff frequency interactive! (e.g. with ginput or
%       rect)
% -->   put npasses in the inputs
% -->   insert possibility to chooose between butterworth and critically
%       damped
% -->   improve drawing a rectangle in every subplot
%
%
% Written by:
% -Ezio Preatoni-   University of Bath, May 2011
%
% Modified by:
% -xxx-             xxx - changed xxx


%% input & input check
npasses= 2;                                 %number of passes of the filter @each iteration. By default =2
c_BW= 1/((2^(1/npasses)-1)^(1/4));          %correcting factor cutoff frequency (for multiple passes, "cascading")

%check arguments number and content
error(nargchk(2,6,nargin,'struct'))
if (isempty(Y) || isempty(fs))
    Yf= [];
    return
end

if nargin==2 || isempty(cutoff)
    %ask for UI filtering parameters
    w1= 8;
    w2= floor(fs/(2*c_BW))-w1;
    prompt={'MIN cutoff freq value [Hz]:','MAX cutoff freq value [Hz]:'};
    name='Limits for cutoff freq range';
    defaultanswer={num2str(w1),num2str(w2+w1)};
    answer=inputdlg(prompt,name,1,defaultanswer);
    [val status] = str2num(char(answer));  %#ok<ST2NM>
    if ~status || ~(val(1)>1 && (val(2)*c_BW)<fs/2 && val(1)<val(2)),
        e=errordlg(['BAD INPUT!! please check that w1+w2 < fs*' num2str(c_BW)],'Input Error');
        uiwait(e); return
    end
    cutoff(1)= val(1);
    cutoff(2)= val(2);
else
    switch nargin
        case 3
            niter= 1;
            pos= 'N';
            graph= 'Y';
        case 4
            pos= 'N';
            graph= 'Y';
        case 5
            graph= 'Y';
    end
end
w1= cutoff(1);
w2= cutoff(2)-w1;

if ~(w1>1 && (w1+w2)*c_BW<fs/2),
    e=errordlg(['BAD INPUT!! please check that w1+w2 < fs*' num2str(c_BW)],'Input Error');
    uiwait(e); return
end


%check data dimensions
[r,c]= size(Y);

if (r>1) && (c>1)
    Yf= NaN(r,c);
    f_BW= NaN(r,c);
    for i= 1:c          %loop over columns
        [Yf(:,i), f_BW(:,i)]= ADF_BUTTERWORTH(Y(:,i),fs,cutoff,niter,pos,graph);
    end
    return
end

if r==1
    Y= Y(:);            %convert row to column
end

len= size(Y,1);        %length of input
time= (1:len)/fs;      %time array 

%% -->filtering (filter: y(n)=a0*(x(n)+2*x(n-1)+x(n-2))+b1*y(n-1)+b2*y(n-2))<
%pre-filtering
[Yf_pre,f_BWpre]= ADAPTIVE_BUTTER(Y,fs,w1,w2,'N',pos);          %pre-filtering with a standard Butterworth

%iterative-adaptive filtering
Yf_temp= Yf_pre;
Yf_iter= NaN(len,niter);
f_BWiter= NaN(len,niter);
for i=1:niter;
    [Yf_temp,f_BWtemp]= ADAPTIVE_BUTTER(Yf_temp,fs,w1,w2,'Y',pos);      %iterative filtering with an adaptive Butterworth
    Yf_iter(:,i)= Yf_temp;
    f_BWiter(:,i)= f_BWtemp;
end
    
%% OUTPUT    
Yf= Yf_iter(:,niter);
f_BW= f_BWiter(:,niter);

% % % % % %% GRAPHICS
% % % % % %graphical representation of pre-filtering with initial cutoff freq (w1+w2)
% % % % % if strcmpi(graph,'Y')
% % % % %     figure('units','normalized','position',[0.1 0.1 0.8 0.8]); hold on;
% % % % %     
% % % % %     subplot(3,1,1); hold on;
% % % % %     plot(time,ones(len,1)*(fs/2),'k');
% % % % %     plot(time,f_BWpre,'b');
% % % % %     plot(time,f_BWiter,'r');
% % % % %     hline(0,'k'); vline(0,'k'); box on; legend('1/2 samp freq','pre cutoff f','adapt cutoff f');
% % % % %     xlabel('t [s]'); ylabel('f [Hz]'); title('SIGNAL AFTER FILTERING');
% % % % %     ylim([-2 (fs/2)+5]);
% % % % %     
% % % % %     subplot(3,1,2);hold on;
% % % % %     plot(time,Y,'ko-');
% % % % %     plot(time,Yf_pre,'bx-');
% % % % %     plot(time,Yf_iter,'r.-');
% % % % %     hline(0,'k'); vline(0,'k'); box on; legend('orig','pre cutoff f','adapt cutoff f');
% % % % %     xlabel('t [s]'); ylabel('y(t)');
% % % % %     
% % % % %     subplot(3,1,3);hold on;
% % % % % % %     plot(time,((Yf_pre-Y)*100)./Y,'bx-');
% % % % %     plot(time,(Yf_pre-Y),'bx-');
% % % % %     err= bsxfun(@minus,Yf_iter,Y);
% % % % % % %    err100= bsxfun(@times,err,100./Y);
% % % % %     avg_err=median(abs(err));
% % % % % % %     plot(time,err100,'r.-');
% % % % %     plot(time,err,'r.-');
% % % % %     hline(avg_err,'m--',num2str(avg_err));
% % % % %     hline(0,'k'); vline(0,'k'); box on; legend('change pre-filt','change adapt-filt');
% % % % %     xlabel('t [s]'); ylabel('diff [N]');
% % % % %     ylim([-200 200]);
% % % % %     
% % % % %     %show curves, spectrograms and cumulative PSD areas of the raw signal
% % % % %     [h1]= fCONTENT(Y,fs);
% % % % %     hline([w1 w2],'k',[],[],h1(3));                 %to do-->improve drawing a rectangle in every subplot
% % % % %  
% % % % %     %show curves, spectrograms and cumulative PSD areas after pre-filtering + initial cutoff bands
% % % % %     [h2]= fCONTENT(Yf_pre,fs);
% % % % %     hline([w1 w2],'k',[],[],h2(3));                 %-->improve drawing a rectangle in every subplot
% % % % %     
% % % % %     %show curves, spectrograms and cumulative PSD areas after pre-filtering + initial cutoff bands
% % % % %     [h3]= fCONTENT(Yf,fs);
% % % % %     hline([w1 w2],'k',[],[],h3(3));                 %-->improve drawing a rectangle in every subplot
end