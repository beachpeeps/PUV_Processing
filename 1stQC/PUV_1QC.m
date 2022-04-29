%% Import raw PUV data and do first QC
%
%
%
% CAUTION: work still TBD for input beginning & inspection times -
%           currently values hard coded
%
%
%% Author: 
% Athina Lange, SIO July 2021
clear all
tic
% input directory where data is - one deloyment per folder
cd /Volumes/LANGE_Passport/PUV/Ruby2D/
files = dir(pwd)
for ii = 1:length(files);names{ii}=files(ii).name;end
deploy_notes = readtable(string(names(contains(names, 'Notes'))))

[~,ii]=sort(deploy_notes.MOP);
deploy_notes=deploy_notes(ii,:);
aa=tabulate(deploy_notes.MOP); 

if length(find(aa(:,2)>1))>0
    mm = find(aa(:,2)>1); % if multiple sensors on single mop line
    for jj = 1:length(mm)
        ll = mm(jj);
        nn = find(aa(ll,1)==deploy_notes.MOP);
        [~,ii]=sort(deploy_notes.Depth(nn));
        deploy_notes(nn,:)=deploy_notes(nn(ii),:);
    end
end

for nn = length(files):-1:1 
    if files(nn).name(1) == '.' | files(nn).isdir == 0
        files(nn)=[];
    end
end
%%
for ii = 1:length(files)
    directory = ['MOP' num2str(deploy_notes.MOP(ii)) '-' num2str(deploy_notes.Depth(ii)) 'm'];
    filename = ['Torrey_Ruby2D_' num2str(deploy_notes.MOP(ii)) '_' num2str(deploy_notes.Depth(ii)) 'm'];
    % input lat/lon of sensor
    LATLON = [deploy_notes.Latitude(ii) deploy_notes.Longitude(ii)];
    % input compass heading of beam 1 of ADV
    rot_angle = [deploy_notes.Heading(ii)];
    % input clock drift between beginning and end of deployment
    clockdrift = [deploy_notes.ClockDrift(ii)];
    [PUV] = PUV_raw_process(directory, filename, LATLON, rot_angle, clockdrift)

end
