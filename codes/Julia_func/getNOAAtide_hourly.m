function [t,verified,tideInfo] = getNOAAtide_hourly(begin_date,end_date,station)
% get tides from NOAA buoy
% [t,verified,tideInfo] = getNOAAtide_hourly('20120204','20140206','9410230')

datum = 'NAVD';
% station = '9410230'; % La Jolla 
time_zone = 'GMT';
units = 'metric';
product = 'hourly_height'; % for more products, see:
% https://tidesandcurrents.noaa.gov/api/

tideInfo = struct('datum',datum,'station',station,...
    'time_zone',time_zone,'units',units,...
    'begin_date',begin_date,'end_date',end_date);

% check to see length of timeseries
ta = datenum(begin_date,'yyyymmdd');
tb = datenum(end_date,'yyyymmdd');

% preassign variables
t = [];
verified = [];

if tb-ta < 30
    S = getSdata(product,begin_date,end_date,datum,station,time_zone,units);
    
    t = S(:,1);
     %t = datenum(t,'yyyy-mm-dd HH:MM');
    t = cellfun(@datenum,t);
   
    verified = cell2mat(S(:,2));  
end

begin_datetemp = begin_date;
end_datetemp = ta+29;
while end_datetemp < tb +29 % if this is true we will have to download in chunks
    
    
    end_datetemp = datestr(end_datetemp,'yyyymmdd');
    
    %%
    S = getSdata(product,begin_datetemp,end_datetemp,datum,station,time_zone,units);
    
    ti = S(:,1);
%     ti = datenum(ti,'yyyy-mm-dd HH:MM');
     %ti = datenum(t,'yyyy-mm-dd HH:MM');
    ti = cellfun(@datenum,ti);
    t = [t; ti];
    ver = cell2mat(S(:,2));
    verified = [verified; ver];
    
    begin_datetemp = datenum(begin_datetemp,'yyyymmdd') + 30;
    begin_datetemp = datestr(begin_datetemp,'yyyymmdd');
    end_datetemp = datenum(begin_datetemp,'yyyymmdd')+29;
    
end

    function S = getSdata(product,begin_date,end_date,datum,station,time_zone,units)
        api = 'https://tidesandcurrents.noaa.gov/api/datagetter?';
        ss = strcat('product=',product,'&application=NOS.COOPS.TAC.WL&begin_date=',...
            begin_date,'&end_date=',end_date,'&datum=',datum,'&station=',station,...
            '&time_zone=',time_zone,'&units=',units,'&format=csv');
        url = [api ss];
        options = weboptions('timeout',20);
        S = webread(url,options);
        S = table2cell(S);
    end
end

