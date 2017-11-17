function segment_mass=read_rra_mass(file)

%[file, pname] = uigetfile('*.log', 'Select C3D file');

var=importdata(file, '\t');
index = 1;

% matrix=[];
for i =1:length(var)
    if strfind(var{i},'new mass')
        ind1=find((var{i}==','), 1 );
        x=length(var{i});
        ind3=index;        
        num(ind3,:)=str2num(var{i}(ind1+12:end));     
        [segment_mass]=num;
        index = index+1;
    end
end

end