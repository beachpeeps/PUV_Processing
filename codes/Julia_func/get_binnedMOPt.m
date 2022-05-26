%% Pick out spectra for testing
% this is probably obsolete??

clear
load('~/Documents/NeoStockdon/data/binnedMOPS.mat','ZBIN','MOP')
% HoLoN = 38;
for HoLoN = 31 %run for HoLo bins #31 m
    for fspN=5:7 %run for frequency spread bins 5-7 Hz
        fspbin = ZBIN.xind(fspN);
        HoLobin = HoLoN;
        % close all
        figure
        tol = eps(0.5);
        xind = find(abs(ZBIN.xind-fspbin) < tol); % find fspread in bin = fspbin-fspbin+spacing Hz
        yind = find(abs(ZBIN.yind-HoLobin) < tol); % find sqrt(HoLo)in bin = HoLobin-HoLobin+spacing
        zind = ZBIN.zind{xind,yind}; % get indices for spectra in this bin
        if numel(zind)>0
            %%
            cmap = parula(numel(zind));
            
            for i=1:numel(zind)
                hp(i) = plot(MOP.frequency,MOP.spec1D(zind(i),:),'linewidth',2,'color',cmap(i,:)) ;
                hold on
            end
            hl = legend(datestr(MOP.time(zind)));
            hl.Box = 'off';
            hl.FontSize = 14;
            colorlegend(hl,hp)
            
            ss1 = [num2str(ZBIN.yind(yind)) ' - ' num2str(ZBIN.yind(yind+1))];
            ss2 = [num2str(ZBIN.xind(xind)) ' - ' num2str(ZBIN.xind(xind+1))];
            
            title(['$\sqrt(H_0L_0)$ = ' ss1 ' m, fspread = ' ss2 ' Hz'], 'interpreter','latex')
            xlabel('Energy (m^2/Hz)')
            ylabel('Frequency (Hz)')
            
            t = MOP.time(zind);
            % save('~/Documents/NeoStockdon/data/binnedMOPS.mat','ZBIN','xind','yind','zind','bimodalmat','hdatabimodal','Ho','Lo','LoC','sqrtHoLo','sqrtHoLoC')
            run_makeInputFilesMOPS(t,xind,yind)
            figureDir = '~/GoogleDriveUCSD/MOPS_Climatology/';
            
            figureName = [figureDir 'specCASES_' num2str(yind,'%02.0f') '_'  num2str(xind,'%02.0f') '.jpeg'];
            print(gcf, '-djpeg', figureName,'-r300');
        end
        clearvars -except ZBIN MOP HoLoN
    end
    close all
end