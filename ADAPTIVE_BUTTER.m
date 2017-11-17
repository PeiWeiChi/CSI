function [y,f_BWs]=ADAPTIVE_BUTTER(x,fs,w1,w2,adapt,pos)

% [varargout]=ADAPTIVE_BUTTER(varargin)
%
% ADAPTIVE_BUTTER.m ... (to be completed).
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
%       - /
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
% /
%
%
% Written by:
% -Ezio Preatoni-   University of Bath, May 2011
%
% Modified by:
% -xxx-             xxx - changed xxx


%% input check
%check arguments number and content
error(nargchk(2,6,nargin,'struct'))
if (isempty(x) || isempty(fs))
    y= [];
    f_BWs= [];
    return
end

switch nargin
    case 2
        w1= 2;
        w2= fs/4;
        adapt= 'Y';
        pos= 'N';
    case 3
        w2= fs/4;
        adapt= 'Y';
        pos= 'N';
    case 4
        adapt= 'Y';
        pos= 'N';
    case 5
        pos= 'N';
end

%check data dimensions
[r,c]= size(x);

if (r>1) && (c>1)
    errordlg('Input data "x" must be a vector, not a matrix!');
    return
end

if r==1
    x= x(:);            %convert row to column
end

len= length(x);

%% filter definition: double pass 2nd order Butterworth,  y(n)=a0*(x(n)+2*x(n-1)+x(n-2))+b1*y(n-1)+b2*y(n-2)
npasses= 2;                                     %number of passes of the filter @each iteration. By default =2
c_BW= 1/((2^(1/npasses)-1)^(1/4));              %correcting factor cutoff frequency (for multiple passes, "cascading")
padl= 2;                                        %pad length= 2= filter order
x= padarray(x,padl,'symmetric','both');         %padding before and after the signal

%% -->FORWARD pass
xff= x;
if strcmpi(adapt,'Y')                           %adaptive filter
    xff1= gradient(xff);                        %1st derifative of pre-filtered data
    xff2= gradient(xff1);                       %2nd derivative of pre-filtered data
    %xffn= abs(xff)/max(abs(xff));
    xff1n= abs(xff1)/max(abs(xff1));
    xff2n= abs(xff2)/max(abs(xff2));
    d= xff1n+xff2n;
    %d= xffn+xff1n;
    %d= xff/max(xff);
    k_BW= d/max(d);                             %adaptive coefficient
    f_BW= (w1+w2*k_BW)*c_BW;                    %adjusted+adaptive cutoff freq
    %smooth the adjusted+adaptive cutoff freq
    [b,a]= butter(2, floor((w1+w2)/2)*c_BW/(fs/2));
    f_BWs= filtfilt(b,a,f_BW);
else
    f_BWs= (w1+w2*ones(size(x)))*c_BW;
end
f_BWs(f_BWs>fs/2)= fs/2;

%calculate coefficients after smoothing
OMEGAc= tan(pi*f_BWs/fs);                       %corrected angular cutoff freq
K1= (2^0.5)*OMEGAc;
K2= OMEGAc.^2;
a0= (K2)./(1+K1+K2);
a1= 2*a0;
a2= a0;
b1= (2*a0).*((1./K2)-1);
b2= 1-(a0+a1+a2+b1);

%filter data (xff is unfiltered data, yff filtered data... ff= filter fwd)
yff= xff;
for n= (1:len)+padl,
    yff(n)= a0(n)*(xff(n)+2*xff(n-1)+xff(n-2))+b1(n)*yff(n-1)+b2(n)*yff(n-2);  %yff stands for "filtered farward"
    if strcmpi(pos,'Y') && yff(n)<0             %if values cannot be <0...
        yff(n)= 0;
    end
end

%% <--BACKWARD pass
xfb= yff(length(xff):-1:1);                     %reverse output from forward
if strcmpi(adapt,'Y')                           %adaptive filter
    xfb1= gradient(xfb);                        %1st derifative of pre-filtered data
    xfb2= gradient(xfb1);                       %2nd derivative of pre-filtered data
    %xfbn= abs(xfb)/max(abs(xfb));
    xfb1n= abs(xfb1)/max(abs(xfb1));
    xfb2n= abs(xfb2)/max(abs(xfb2));
    d= xfb1n+xfb2n;
    %d= xfbn+xfb1n;
    %d= xfb/max(xfb);
    k_BW= d/max(d);                             %adaptive coefficient
    f_BW= (w1+w2*k_BW)*c_BW;                    %adjusted+adaptive cutoff freq
    %smooth the adjusted+adaptive cutoff freq
    [b,a]= butter(2, floor((w1+w2)/2)*c_BW/(fs/2));
    f_BWs= filtfilt(b,a,f_BW);
else
    f_BWs= (w1+w2*ones(size(x)))*c_BW;
end
f_BWs(f_BWs>fs/2)= fs/2;
    
%calculate coefficients after smoothing
OMEGAc= tan(pi*f_BWs/fs);                       %corrected angular cutoff freq
K1= (2^0.5)*OMEGAc;
K2= OMEGAc.^2;
a0= (K2)./(1+K1+K2);
a1= 2*a0;
a2= a0;
b1= (2*a0).*((1./K2)-1);
b2= 1-(a0+a1+a2+b1);

yfb= xfb;
%filter data
for n= (1:len)+padl,
    yfb(n)= a0(n)*(xfb(n)+2*xfb(n-1)+xfb(n-2))+b1(n)*yfb(n-1)+b2(n)*yfb(n-2); %xfb stands for "filtered backward"
    if strcmpi(pos,'Y') && yfb(n)<0             %if values cannot be <0...
        yfb(n)= 0;
    end
end
yfb= yfb(length(yfb):-1:1);
f_BWs= f_BWs(length(yfb):-1:1);

%% OUTPUT
y= yfb((1:len)+padl);
f_BWs= f_BWs((1:len)+padl);

