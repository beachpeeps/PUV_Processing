function adjustVspace(ax)
nrow = length(ax);

for nn = 1:nrow-1
    ax(nn).XTickLabel = [];

end


% set all yaxes equal
ymax = nan(1,nrow);
ymin = nan(1,nrow);

for nn = 1:nrow
    ymax(nn) = ax(nn).YLim(2);
    ymin(nn) = ax(nn).YLim(1);
end

ymax = max(ymax);
ymin = min(ymin);

for nn=1:nrow
    ax(nn).YLim = [ymin ymax];
end


% make boxes closer
% for nn=1:nrow
%     ybot(nn) = ax(nn).Position(2);
%     height(nn) = ax(nn).Position(4);
% end
% 
% space = ybot(1)-ybot(2)-height(1);
% nspace = nrow-1;
totfill = ax(1).Position(2)+ax(1).Position(4)-ax(nrow).Position(2);
newheight = totfill/nrow;

for nn=1:nrow
    ax(nn).Position(4) = newheight;
end

    