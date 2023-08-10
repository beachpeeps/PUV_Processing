function [hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,varargin)
% [hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,varargin)
% make a figure with a specified width, height, and number of rows and
% columns. Output is the figure handle and the axes, in a structure.
%
% To add the options to the function, call (...'units','inches') etc.
%
% Example:
%
% [hFig, ax] = makeFig(6,5,2,2,'units','inches','cbar','bottom')
% will make a figure with 2 rows and 2 columns that is 6 inches wide and 5
% inches tall, with a colorbar axis on the bottom.

% copyright 2018, Julia Fiedler jfiedler@ucsd.edu


% These are the default options
options.units = 'inches';
options.cbar   = 'none';
options.margin = 0.05;
options.bmargin = 0.1;
options.space = 0.04;
options.hspace = 0.03;
options.visible = 'on';
options.axesfontsize = 14;

options = parseOptions( options , varargin );


hFig = figure('visible',options.visible);
set(0,'defaultaxesfontsize',options.axesfontsize)
set(0,'defaultaxeslinewidth',1)
hFig.PaperUnits = options.units;
hFig.PaperSize = [figwidth figheight];
hFig.PaperPosition = [0 0 figwidth figheight];

hFig.Position = hFig.PaperPosition*100;
hFig.Position(1:2) = [-1200 250];
hFig.PaperPositionMode = 'manual';
hFig.InvertHardcopy = 'off';


% shortcuts
margin = options.margin;
bmargin = options.bmargin;
space = options.space;
hspace = options.hspace;


totwidth = 1-3*margin;
width = (totwidth-(ncol-1)*space)./ncol;
totheight = 1-bmargin-2*margin-(nrow-1)*hspace;
height = totheight./nrow;

left = zeros(1,ncol);
bottom = zeros(1,nrow);

for nncol=1:ncol
    left(nncol) = 2*margin+(nncol-1)*space+(nncol-1)*width;
end

for nnrow=1:nrow
    bottom(nnrow) = bmargin+(nnrow-1)*hspace+(nnrow-1)*height;
end
bottom = fliplr(bottom);


% Add in the colorbar
% TODO: Allow for side bar (on either side)

cbarheight = 0;
cbaradjust = 0;

switch options.cbar
    case 'bottom'
        cbarheight = 0.1;
        ax(nrow*ncol+1) = axes('position',[left(1) .05 totwidth cbarheight]);
        totheight = totheight-cbarheight-2*hspace;
        cbaradjust = cbarheight+2*hspace;
    case 'top'
        cbarheight = 0.1;
        ax(nrow*ncol+1) = axes('position',[left(1) 0.9 totwidth cbarheight]);
        totheight = totheight-cbarheight-2*hspace;
        cbaradjust = 0;
    case 'none'
        cbaradjust = 0;
        totheight = totheight-cbarheight-2*hspace;
end


height = totheight./nrow;
bottom = bottom+cbaradjust;

% Adjust the rest of the axes
%[left bottom width height]
for nnrow=1:nrow
    for i=1:ncol
        ax(i+(nnrow-1)*ncol) = axes('position',[left(i) bottom(nnrow) width height]);
        box on

    end
end

switch 'cbar'
    case 'top'
ax(nrow*ncol+1).YTickLabel = [];
ax(nrow*ncol+1).XTickLabel = [];
end


    function options = parseOptions( defaults , cellString )
        
        %INPUT PARSING
        p = inputParser;
        p.KeepUnmatched = true;
        
        names = fieldnames( defaults );
        for ii = 1 : length(names)
            %
            addOptional( p , names{ii}, defaults.(names{ii}));
            %
        end
        %
        parse( p , cellString{:} );
        options = p.Results;
        %
    end

end
