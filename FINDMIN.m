function[m,r]=FINDMIN(data,ini,fin)


% [m,r]=FINDMIN(data), finds minima and the corresponding rows for every column
% of the input data matrix.
% 
% Input:
%       - data              is the [r x c] input data matrix
%       - ini               is the [1 x c] vector containing the initial
%                           row of every column which the search for minimum value must be started from 
%       - fin               is the [1 x c] vector containing the final
%                           row of every column which the search for minimum value must be stopped 
% 
% Output:
%       - m                 is the [1 x c] vector containing minima
%       - r                 is the [1 x c] vector containing rows where
%                           minima where found
%
% Called functions:
%       /
%
% -Ezio Preatoni- Politecnico di Milano, October 2006
% -Ezio Preatoni- Politecnico di Milano, April 2008 - search for min can be
%                       done in different sections for each column of input
%                       matrix


%% controls over input arguments

if nargin<1,
    errordlg('Not enough input arguments','FINDMIN.m - input error');
end
[rows cols]=size(data);
if nargin==1,
    ini=ones(1,cols);
    fin=ones(1,cols)*rows;
elseif nargin==2,
    fin=ones(1,cols)*rows;
end


%% finds minima

m=NaN(1,cols);
r=NaN(1,cols);

for i=1:cols,
    m(i)=min(data(ini(i):fin(i),i),[],1);
    if ~isnan(m(i)),
        r(i)= find(data(ini(i):fin(i),i)==m(i),1);
    else
        r(i)= NaN;
    end
end