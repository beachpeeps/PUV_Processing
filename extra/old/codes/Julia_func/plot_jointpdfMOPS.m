function [zdata, hdata, ZBIN] = plot_jointpdfMOPS(X,Y,Z,HC,fileName,printYN)
% joinpdfMOPS(xdatalabel,xdata,ydatalabel,ydata,Z,hcolor,HC)
% plots the joint pdf of Z as binned into X and Y, with histograms of values
% used to form the bins on the x and y axes. Histogram bars are colored
% by HC.
%
% INPUTS: 
% fileName: location where output figure is to be saved
% printYN: 0 if not saved, 1 if saved
%
% ***Input for the main plots (X,Y,Z)***
% Z, structure with the following components
%   .label = string of variable name, example: '$\sqrt (HoLo)$'
%   .labelInterpreter = 'tex' or 'latex'
%   .array = variable to put on the z axis, size(1, number of obs)
%   .binMethod = 'max', 'min', or 'mean'
%   .scale = 'log' or 'linear'
% 
% X/Y, structure with the following components:
%   .label = string of variable name, example: 'Lowest f peak (Hz)'
%   .labelInterpreter = 'tex' or 'latex'
%   .array = variable to put on the x or y axis, size(1, number of obs)
%   .bin = bins for the x/y axis, example: [245:3:295]
%   .scale = 'log' or 'linear'
% 
% *** Input for the histogram colors (HC) ****
% HC, structure with the following components
%   .xtick = array to set bounds of colorbar, example: [0 1]
%   .xticklabel = label for the colorbar ticks, ex: {'summer','winter'}
%                 **note these must match the # of ticks!!** 
% 	.label = string of variable name, example: 'Avg. Season'
% 	.array = variable to color the histograms, size(1, number of obs)
%   .binMethod = 'mean' or 'fraction'  
%   .yLimit = maximum count number, for showing extremes
%
% OUTPUTS:
% zdata: binned Z.array
% hdata: binned HC.array
 
% Julia Fiedler, 2019 jfiedler@ucsd.edu


%% BIN DATA
yind = discretize(Y.array,Y.bin); % bin ydata
xind = discretize(X.array,X.bin); % bin xdata
nx = length(X.bin);
ny = length(Y.bin);

% preallocate variables for binning
zdata = nan(nx,ny);
counts = nan(nx,ny);
hdata = nan(nx,ny);
indZ = cell(nx,ny);
%
for i=1:length(X.bin)
    for j = 1:length(Y.bin)
        zind = intersect(find(xind==i),find(yind==j));
        counts(i,j) = numel(zind);
        indZ(i,j) = {zind};
        if numel(zind)>0
            switch Z.binMethod
                case 'mean'
                    zdata(i,j) = nanmean(Z.array(zind));
                case 'max'
                    zdata(i,j) = nanmax(Z.array(zind));
                case 'min'
                    zdata(i,j) = nanmin(Z.array(zind));
                case 'fraction'
                    zdata(i,j) = sum(Z.array(zind))./numel(zind);
            end

            switch HC.binMethod
                case 'mean'
                    hdata(i,j) = nanmean(HC.array(zind));
                case 'fraction'
                    hdata(i,j) = sum(HC.array(zind))./numel(zind);
            end
        end
    end
end
%%


disp(['Making plot of ' Z.label ' as binned by ' X.label ' and ' Y.label]);
figwidth = 7;
figheight = 7;
nrow = 2;
ncol = 2;
units = 'inches';
cbar = 1;
[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,units,cbar);

ax(1).Position = [0.1 0.71 0.6 0.2];
ax(3).Position = [0.1 0.1 0.6 0.6];
ax(4).Position = [0.71 0.1 0.2 0.6];
ax(2).Visible = 'off';

ax(5) = axes('Position',ax(3).Position);
ax(5).Visible = 'off';


cmap = cmocean('-gray');
axes(ax(3))
pcolor(X.bin,Y.bin,zdata'); shading flat

ax(3).XLabel.Interpreter = X.labelInterpreter;
ax(3).XLabel.String = X.label;
ax(3).YLabel.Interpreter = Y.labelInterpreter;
ax(3).YLabel.String = Y.label;


colormap(ax(3),Z.cmap)

axes(ax(5))
contour(X.bin,Y.bin,counts','linewidth',1)
ax(5).Color = 'none';
colormap(ax(5),gray);

axes(ax(1))
hx = histogram(X.array,'BinEdges',X.bin);
bx = bar(X.bin(1:end-1),hx.BinCounts,1);
bx.FaceColor = 'flat';

axes(ax(4))
hy = histogram(Y.array,'BinEdges',Y.bin);
by = bar(Y.bin(1:end-1),hy.BinCounts,1);
by.FaceColor = 'flat';

switch HC.binMethod
    case 'mean'
        colorindX = nanmean(hdata,2)/max(hdata(:));
        colorindY = nanmean(hdata)/max(hdata(:));
    case 'fraction'
        colorindX = nanmean(hdata,2);
        colorindY = nanmean(hdata);
end

colorindX(isnan(colorindX)) = 0;
colorindY(isnan(colorindY)) = 0;

% colorindX(colorindX == 0) = 1;
% 
% colorindY(isnan(colorindY)) = 1;
% colorindY(colorindY == 0) = 1;

colorindXY = [colorindX(:) ; colorindY(:)];
colorindXYind = round(colorindXY./max(colorindXY)*256);
colorindXYind(colorindXYind==0) = 1;

colorX = colorindXYind(1:length(colorindX)-1);
colorY = colorindXYind(length(colorindX)+1:end-1);

bx.CData = cmap(colorX,:);
by.CData = cmap(colorY,:);

ax(4).View = [90 -90];
if isfield(HC,'yLimit')
    ax(4).YLim(2) = HC.yLimit;
end

ax(6) = axes('Position',[0.71 0.85 0.2 0.03]);
ax(6).Visible = 'off'; 
plot(hdata)
ax(6).CLim = [0 max(colorindXY)];
colormap(ax(6),cmap)
c = colorbar(ax(6),'location','northoutside');
c.Position = ax(6).Position;
ax(6).Position([3 4]) = [0 0];
ax(6).XTickLabel = []; ax(6).YTickLabel = [];

c.Ticks = ax(6).CLim;
c.TickLabels = {num2str(c.Limits(1)), num2str(round(c.Limits(2),1),'%2.1f')};
c.Label.String = HC.label;

c2 = colorbar(ax(3),'location','northoutside');
c2.Position = [0.4729    0.5879    0.2000    0.0300];
c2.Label.Interpreter = Z.labelInterpreter;
c2.Label.String = Z.label;
c2.TickDirection = 'out';

for i = [1 4]
ax(i).Box = 'off'; ax(i).Color = 'none';
end

ax(1).XTickLabel = [];
for i= [1 3 5]
    ax(i).XScale = X.scale;
end

linkaxes([ax(1),ax(3),ax(5)],'x')


% in lieu of linkaxes([ax(1),ax(3),ax(5)],'y'):
ax(3).YScale = Y.scale;
ax(5).YScale = Y.scale;

ax(4).XLim = ax(3).YLim;
ax(4).XScale = ax(3).YScale;
ax(4).XTickLabel = [];
ax(4).YLabel.String = 'counts';
ax(1).YLabel.String = 'counts';

ax(3).ColorScale = Z.scale;

if printYN
    print(hFig,'-djpeg',fileName,'-r300')
    print(hFig,'-dpdf',fileName,'-r300')

end

ZBIN.xind = X.bin;
ZBIN.zind = indZ;
ZBIN.yind = Y.bin;
