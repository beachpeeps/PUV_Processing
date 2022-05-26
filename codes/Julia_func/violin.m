% violin plot
function [hv, hm] = violin(x,fac,center,alp,fc,ax)
% calculate kernal density


% pd = fitdist(Sig(:),'Normal');
% % ci = paramci(pd);
% cilo = pd.mu-2*pd.sigma;
% ciup = pd.mu+2*pd.sigma;
% y = icdf(pd,[0.25 0.75]);

[f,xi,bw] = ksdensity(x);



f = fac*f/max(f);
MED = nanmedian(x);
MX = nanmean(x);

yy = [xi flipud(xi)];
xx = [f -flipud(f)]+center;
hv = fill(xx,yy,fc,'EdgeColor','none','FaceAlpha',alp,'parent',ax);
hold on

hm  = scatter(center,MED,100,'o','linewidth',2,'MarkerFaceColor',fc,'MarkerEdgeColor','none');
% xx = [-0.2 -0.2 0.2 0.2]+center;
% yy = [ci(:,1); flipud(ci(:,1))];
% yy = pd.mu-2*pd.sigma
%     
% fill(xx,yy,'r')

% hold on
% hv(2) = fill(-[f flipud(f)],[xi flipud(xi)],'g','EdgeColor','none','FaceAlpha',alp);

% for i=1:2
%     hv(i).EdgeColor = hv(i).FaceColor;
%     hv(i).EdgeAlpha = hv(i).FaceAlpha;
% end
