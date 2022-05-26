function colorlegend(legendhandle,plothandle)
%% color legend strings
n = 1;
for i=1:numel(plothandle)
    legendhandle.String{n} = ['\color[rgb]{' num2str(plothandle(i).Color) '} ' legendhandle.String{n}];
    n = n+1;
end