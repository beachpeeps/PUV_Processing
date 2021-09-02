% %% 
% 
% %% Initial QC
% % PUV_1QC
% load Torrey_1920_582_10_processed.mat
% %% Rotate to shorenormal coordinates
% [PUV.ShorenormalCoord.U, PUV.ShorenormalCoord.V, shorenormal] = rotate_shorenormal(PUV.BuoyCoord.U, PUV.BuoyCoord.V);
% 
% %%
% tic
% % 66cm sand to top of pressure port on deployment. 78cm on recovery. 
% doffp = 0.72; % splitting the difference between the two
%     doffu = doffp+0.371; %adding 37.1 cm to u,v sample point
% 
%     PUV_work = PUV;
% 
% %% Interp NaNs if less than 2sec gap
% PUV_work.P = fillmissing(PUV.P,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
% PUV_work.BuoyCoord.U = fillmissing(PUV.BuoyCoord.U,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
% PUV_work.BuoyCoord.V = fillmissing(PUV.BuoyCoord.V,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
% PUV_work.BuoyCoord.W = fillmissing(PUV.BuoyCoord.W,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
% PUV_work.ShorenormalCoord.U = fillmissing(PUV.ShorenormalCoord.U,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
% PUV_work.ShorenormalCoord.V = fillmissing(PUV.ShorenormalCoord.V,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
% PUV_work.T = fillmissing(PUV.T,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
% toc
% %% Split data into hour segments
% hourlen = (60*60*PUV_work.fs)
% segtotal = floor(length(PUV_work.time)/hourlen)
% PUV_work.time(segtotal*hourlen+1:end)=[];
% PUV_work.P(segtotal*hourlen+1:end)=[];
% PUV_work.BuoyCoord.U(segtotal*hourlen+1:end)=[];
% PUV_work.BuoyCoord.V(segtotal*hourlen+1:end)=[];
% PUV_work.BuoyCoord.W(segtotal*hourlen+1:end)=[];
% PUV_work.ShorenormalCoord.U(segtotal*hourlen+1:end)=[];
% PUV_work.ShorenormalCoord.V(segtotal*hourlen+1:end)=[];
% 
% PUV_work.time = reshape(PUV_work.time, hourlen, segtotal);
% PUV_work.P = reshape(PUV_work.P, hourlen, segtotal);
% PUV_work.BuoyCoord.U = reshape(PUV_work.BuoyCoord.U, hourlen, segtotal);
% PUV_work.BuoyCoord.V = reshape(PUV_work.BuoyCoord.V, hourlen, segtotal);
% PUV_work.BuoyCoord.W = reshape(PUV_work.BuoyCoord.W, hourlen, segtotal);
% 
% PUV_work.ShorenormalCoord.U = reshape(PUV_work.ShorenormalCoord.U, hourlen, segtotal);
% PUV_work.ShorenormalCoord.V = reshape(PUV_work.ShorenormalCoord.V, hourlen, segtotal);
% 
% P = PUV_work.P;
% U = PUV_work.BuoyCoord.U;
% V = PUV_work.BuoyCoord.V;
% W = PUV_work.BuoyCoord.W;
% fs = PUV_work.fs;
% time = PUV_work.time;
% toc
% %% Get tides to check
% disp('Getting tides for a sanity check')
% datum = 'NAVD';
% station =  '9410230'; % La Jolla
% [tide_t_verified,tide_t_predicted,tide_z_verified,tide_z_predicted,tideInfo] = getNOAAtide(time(1,1),time(end,end),station,datum);
% % A1=split(importdata('Nov.txt', ' '));
% % A2=split(importdata('Dec.txt', ' '));
% % A3=split(importdata('Jan.txt', ' '));
% % A4=split(importdata('Feb.txt', ' '));
% % A5=split(importdata('Mar.txt', ' '));
% % 
% % tide_time = [datetime(A1(2:end,1), 'InputFormat', 'yyyy/MM/dd') + timeofday(datetime(A1(2:end,3), 'InputFormat', 'HH:mm'));...
% % 	datetime(A2(2:end,1), 'InputFormat', 'yyyy/MM/dd') + timeofday(datetime(A2(2:end,3), 'InputFormat', 'HH:mm'));...
% % 	datetime(A3(2:end,1), 'InputFormat', 'yyyy/MM/dd') + timeofday(datetime(A3(2:end,3), 'InputFormat', 'HH:mm'));...
% % 	datetime(A4(2:end,1), 'InputFormat', 'yyyy/MM/dd') + timeofday(datetime(A4(2:end,3), 'InputFormat', 'HH:mm'));...
% % 	datetime(A5(2:end,1), 'InputFormat', 'yyyy/MM/dd') + timeofday(datetime(A5(2:end,3), 'InputFormat', 'HH:mm'))];
% % 
% % tide = [str2num(char(A1(2:end,4)));...
% %     str2num(char(A2(2:end,4)));...
% %     str2num(char(A3(2:end,4)));...
% %     str2num(char(A4(2:end,4)));...
% %     str2num(char(A5(2:end,4)))];
% 
% %% Get useable values
% disp('Getting wave stats, writing to PUV_process_ENU')
% tic
% Uprime = PUV_work.BuoyCoord.U;
% Vprime = PUV_work.BuoyCoord.V;
% for ii = 1:segtotal
%     if isempty(find(isnan(P(:,ii))))
%         [PUV_process_ENU(ii)] = vector_wave_stats_mspec(Uprime(:,ii), Vprime(:,ii), P(:,ii), doffp, doffu);
%     end
% end
% toc 
% %%
% disp('Getting wave stats, writing to PUV_process')
% tic
% Uprime = PUV_work.ShorenormalCoord.U;
% Vprime = PUV_work.ShorenormalCoord.V;
% for ii = 1866:segtotal
%     if isempty(find(isnan(P(:,ii))))
%         [PUV_process(ii)] = vector_wave_stats_mspec(Uprime(:,ii), Vprime(:,ii), P(:,ii), doffp, doffu);
%     end
% end
% toc

%% load variables up to this point
load temp

%% Extract values
disp('Extracting values from PUV_process structure')
tic
for ii = 1:segtotal
    if ~isempty(PUV_process(ii).ztest)
        
        Spec.SSE(:,ii) = PUV_process(ii).Spec.SSE;
        Spec.Spp(:,ii) = PUV_process(ii).Spec.Spp;
        Spec.Suu(:,ii) = PUV_process(ii).Spec.Suu;
        Spec.Svv(:,ii) = PUV_process(ii).Spec.Svv;
        Spec.Suv(:,ii) = PUV_process(ii).Spec.Suv;
        
        Ztest.ztest_all(:,ii) = PUV_process(ii).ztest.ztest_all;
        Ztest.ztest_ss(:,ii) = PUV_process(ii).ztest.ztest_ss;
        Ztest.ztest_ig(:,ii) = PUV_process(ii).ztest.ztest_ig;
        Ztest.ztest_ss_sum(ii) = PUV_process(ii).ztest.ztest_ss_sum;
        Ztest.ztest_ig_sum(ii) = PUV_process(ii).ztest.ztest_ig_sum;
        
        Hsig.Hs(ii) = PUV_process(ii).Hsig.Hs;
        Hsig.Hsss(ii) = PUV_process(ii).Hsig.Hs_ss;
        Hsig.Hsig(ii) = PUV_process(ii).Hsig.Hs_ig;
        
        Dir.dir1(:,ii) = PUV_process(ii).dir.dir1;
        Dir.spread1(:,ii) = PUV_process(ii).dir.spread1;
        Dir.dir2(:,ii) = PUV_process(ii).dir.dir2;
        Dir.spread2(:,ii) = PUV_process(ii).dir.spread2;
        Dir.dir1_ss_sum(ii) = PUV_process(ii).dir.dir1_ss_sum;
        Dir.spread1_ss_sum(ii) = PUV_process(ii).dir.spread1_ss_sum;
        Dir.dir2_ss_sum(ii) = PUV_process(ii).dir.dir2_ss_sum;
        Dir.spread2_ss_sum(ii) = PUV_process(ii).dir.spread2_ss_sum;
        
        FC.a1(:,ii) = PUV_process(ii).FC.a1;
        FC.b1(:,ii) = PUV_process(ii).FC.b1;
        FC.a2(:,ii) = PUV_process(ii).FC.a2;
        FC.b2(:,ii) = PUV_process(ii).FC.b2;
        
        RS.Sxx(:,ii) = PUV_process(ii).RS.Sxx;
        RS.Syy(:,ii) = PUV_process(ii).RS.Syy;
        RS.Sxy(:,ii) = PUV_process(ii).RS.Sxy;
        RS.Sxx_ss(ii) = PUV_process(ii).RS.Sxx_ss;
        RS.Syy_ss(ii) = PUV_process(ii).RS.Syy_ss;
        RS.Sxy_ss(ii) = PUV_process(ii).RS.Sxy_ss;
        
        Tp(ii) = 1./PUV_process(ii).ids.fpeak;
        
        Eflux.Fpos(:,ii) = PUV_process(ii).Eflux.posX;
        Eflux.Fneg(:,ii) = PUV_process(ii).Eflux.negX;
        Eflux.Epos(:,ii) = PUV_process(ii).Eflux.posX./PUV_process(ii).Eflux.Cg;
        Eflux.Eneg(:,ii) = PUV_process(ii).Eflux.negX./PUV_process(ii).Eflux.Cg;
        Eflux.depth(ii) = PUV_process(ii).Eflux.depth;
    else
        
        Spec.SSE(:,ii) = NaN(size(Spec.SSE(:,ii-1),1),1);
        Spec.Spp(:,ii) = NaN(size(Spec.Spp(:,ii-1),1),1);
        Spec.Suu(:,ii) = NaN(size(Spec.Suu(:,ii-1),1),1);
        Spec.Svv(:,ii) = NaN(size(Spec.Svv(:,ii-1),1),1);
        Spec.Suv(:,ii) = NaN(size(Spec.Suv(:,ii-1),1),1);
        
        Ztest.ztest_all(:,ii) = NaN(size(Ztest.ztest_all(:,ii-1),1),1);
        Ztest.ztest_ss(:,ii) = NaN(size(Ztest.ztest_ss(:,ii-1),1),1);
        Ztest.ztest_ig(:,ii) = NaN(size(Ztest.ztest_ig(:,ii-1),1),1);
        Ztest.ztest_ss_sum(ii) = NaN;
        Ztest.ztest_ig_sum(ii) = NaN;
        
        Hsig.Hs(ii) = NaN;
        Hsig.Hsss(ii) = NaN;
        Hsig.Hsig(ii) = NaN;
        
        Dir.dir1(:,ii) = NaN(size(Dir.dir1(:,ii-1),1),1);
        Dir.spread1(:,ii) = NaN(size(Dir.spread1(:,ii-1),1),1);
        Dir.dir2(:,ii) = NaN(size(Dir.dir2(:,ii-1),1),1);
        Dir.spread2(:,ii) = NaN(size(Dir.spread2(:,ii-1),1),1);
        Dir.dir1_ss_sum(ii) = NaN;
        Dir.spread1_ss_sum(ii) = NaN;
        Dir.dir2_ss_sum(ii) = NaN;
        Dir.spread2_ss_sum(ii) = NaN;
                
        FC.a1(:,ii) = NaN(size(FC.a1(:,ii-1),1),1);
        FC.b1(:,ii) = NaN(size(FC.b1(:,ii-1),1),1);
        FC.a2(:,ii) = NaN(size(FC.a2(:,ii-1),1),1);
        FC.b2(:,ii) = NaN(size(FC.b2(:,ii-1),1),1);
        
        RS.Sxx(:,ii) = NaN(size(RS.Sxx(:,ii-1),1),1);
        RS.Syy(:,ii) = NaN(size(RS.Syy(:,ii-1),1),1);
        RS.Sxy(:,ii) = NaN(size(RS.Sxy(:,ii-1),1),1);
        RS.Sxx_ss(ii) = NaN;
        RS.Syy_ss(ii) = NaN;
        RS.Sxy_ss(ii) = NaN;
        
        Tp(ii) = NaN;
        
        Eflux.Epos(:,ii) = NaN(size(Eflux.Epos(:,ii-1),1),1);
        Eflux.Eneg(:,ii) = NaN(size(Eflux.Eneg(:,ii-1),1),1);
        Eflux.Fpos(:,ii) = NaN(size(Eflux.Fpos(:,ii-1),1),1);
        Eflux.Fneg(:,ii) = NaN(size(Eflux.Fneg(:,ii-1),1),1);
        Eflux.depth(ii) = NaN;
    end    
end
i_ig = PUV_process(1).ids.i_ig;
i_swell = PUV_process(1).ids.i_swell;
fm = PUV_process(1).Spec.fm;
df = fm(2)-fm(1);
toc
disp('Computing reflection coefficient')
R2_ig = NaN(1,segtotal);
R2_ss = NaN(1,segtotal);
for ii = 1:segtotal
    if ~isempty(PUV_process(ii).Eflux)
        
        R2_ig(ii) = sum(Eflux.Fneg(i_ig,ii).*df)/sum(Eflux.Fpos(i_ig,ii).*df);
        R2_ss(ii) = sum(Eflux.Fneg(i_swell,ii).*df)/sum(Eflux.Fpos(i_swell,ii).*df);
    
    end
    if R2_ig(ii)==0; R2_ig(ii)=NaN;end
    if R2_ss(ii)==0; R2_ss(ii)=NaN;end
end
toc

%% Extract values ENU FOR
disp('Extracting values from PUV_process_ENU structure')
tic
for ii = 1:segtotal
    if ~isempty(PUV_process_ENU(ii).ztest)
        
        Spec_ENU.SSE(:,ii) = PUV_process_ENU(ii).Spec.SSE;
        Spec_ENU.Spp(:,ii) = PUV_process_ENU(ii).Spec.Spp;
        Spec_ENU.Suu(:,ii) = PUV_process_ENU(ii).Spec.Suu;
        Spec_ENU.Svv(:,ii) = PUV_process_ENU(ii).Spec.Svv;
        Spec_ENU.Suv(:,ii) = PUV_process_ENU(ii).Spec.Suv;
        
        Ztest_ENU.ztest_all(:,ii) = PUV_process_ENU(ii).ztest.ztest_all;
        Ztest_ENU.ztest_ss(:,ii) = PUV_process_ENU(ii).ztest.ztest_ss;
        Ztest_ENU.ztest_ig(:,ii) = PUV_process_ENU(ii).ztest.ztest_ig;
        Ztest_ENU.ztest_ss_sum(ii) = PUV_process_ENU(ii).ztest.ztest_ss_sum;
        Ztest_ENU.ztest_ig_sum(ii) = PUV_process_ENU(ii).ztest.ztest_ig_sum;
        
        Hsig_ENU.Hs(ii) = PUV_process_ENU(ii).Hsig.Hs;
        Hsig_ENU.Hsss(ii) = PUV_process_ENU(ii).Hsig.Hs_ss;
        Hsig_ENU.Hsig(ii) = PUV_process_ENU(ii).Hsig.Hs_ig;
        
        Dir_ENU.dir1(:,ii) = PUV_process_ENU(ii).dir.dir1;
        Dir_ENU.spread1(:,ii) = PUV_process_ENU(ii).dir.spread1;
        Dir_ENU.dir2(:,ii) = PUV_process_ENU(ii).dir.dir2;
        Dir_ENU.spread2(:,ii) = PUV_process_ENU(ii).dir.spread2;
        Dir_ENU.dir1_ss_sum(ii) = PUV_process_ENU(ii).dir.dir1_ss_sum;
        Dir_ENU.spread1_ss_sum(ii) = PUV_process_ENU(ii).dir.spread1_ss_sum;
        Dir_ENU.dir2_ss_sum(ii) = PUV_process_ENU(ii).dir.dir2_ss_sum;
        Dir_ENU.spread2_ss_sum(ii) = PUV_process_ENU(ii).dir.spread2_ss_sum;
        
        FC_ENU.a1(:,ii) = PUV_process_ENU(ii).FC.a1;
        FC_ENU.b1(:,ii) = PUV_process_ENU(ii).FC.b1;
        FC_ENU.a2(:,ii) = PUV_process_ENU(ii).FC.a2;
        FC_ENU.b2(:,ii) = PUV_process_ENU(ii).FC.b2;
        
        RS_ENU.Sxx(:,ii) = PUV_process_ENU(ii).RS.Sxx;
        RS_ENU.Syy(:,ii) = PUV_process_ENU(ii).RS.Syy;
        RS_ENU.Sxy(:,ii) = PUV_process_ENU(ii).RS.Sxy;
        RS_ENU.Sxx_ss(ii) = PUV_process_ENU(ii).RS.Sxx_ss;
        RS_ENU.Syy_ss(ii) = PUV_process_ENU(ii).RS.Syy_ss;
        RS_ENU.Sxy_ss(ii) = PUV_process_ENU(ii).RS.Sxy_ss;
        
        Tp_ENU(ii) = 1./PUV_process_ENU(ii).ids.fpeak;
        
        Eflux_ENU.Fpos(:,ii) = PUV_process_ENU(ii).Eflux.posX;
        Eflux_ENU.Fneg(:,ii) = PUV_process_ENU(ii).Eflux.negX;
        Eflux_ENU.Epos(:,ii) = PUV_process_ENU(ii).Eflux.posX./PUV_process_ENU(ii).Eflux.Cg;
        Eflux_ENU.Eneg(:,ii) = PUV_process_ENU(ii).Eflux.negX./PUV_process_ENU(ii).Eflux.Cg;
        Eflux_ENU.depth(ii) = PUV_process_ENU(ii).Eflux.depth;
    else
        
        Spec_ENU.SSE(:,ii) = NaN(size(Spec_ENU.SSE(:,ii-1),1),1);
        Spec_ENU.Spp(:,ii) = NaN(size(Spec_ENU.Spp(:,ii-1),1),1);
        Spec_ENU.Suu(:,ii) = NaN(size(Spec_ENU.Suu(:,ii-1),1),1);
        Spec_ENU.Svv(:,ii) = NaN(size(Spec_ENU.Svv(:,ii-1),1),1);
        Spec_ENU.Suv(:,ii) = NaN(size(Spec_ENU.Suv(:,ii-1),1),1);
        
        Ztest_ENU.ztest_all(:,ii) = NaN(size(Ztest_ENU.ztest_all(:,ii-1),1),1);
        Ztest_ENU.ztest_ss(:,ii) = NaN(size(Ztest_ENU.ztest_ss(:,ii-1),1),1);
        Ztest_ENU.ztest_ig(:,ii) = NaN(size(Ztest_ENU.ztest_ig(:,ii-1),1),1);
        Ztest_ENU.ztest_ss_sum(ii) = NaN;
        Ztest_ENU.ztest_ig_sum(ii) = NaN;
        
        Hsig_ENU.Hs(ii) = NaN;
        Hsig_ENU.Hsss(ii) = NaN;
        Hsig_ENU.Hsig(ii) = NaN;
        
        Dir_ENU.dir1(:,ii) = NaN(size(Dir_ENU.dir1(:,ii-1),1),1);
        Dir_ENU.spread1(:,ii) = NaN(size(Dir_ENU.spread1(:,ii-1),1),1);
        Dir_ENU.dir2(:,ii) = NaN(size(Dir_ENU.dir2(:,ii-1),1),1);
        Dir_ENU.spread2(:,ii) = NaN(size(Dir_ENU.spread2(:,ii-1),1),1);
        Dir_ENU.dir1_ss_sum(ii) = NaN;
        Dir_ENU.spread1_ss_sum(ii) = NaN;
        Dir_ENU.dir2_ss_sum(ii) = NaN;
        Dir_ENU.spread2_ss_sum(ii) = NaN;
                
        FC_ENU.a1(:,ii) = NaN(size(FC_ENU.a1(:,ii-1),1),1);
        FC_ENU.b1(:,ii) = NaN(size(FC_ENU.b1(:,ii-1),1),1);
        FC_ENU.a2(:,ii) = NaN(size(FC_ENU.a2(:,ii-1),1),1);
        FC_ENU.b2(:,ii) = NaN(size(FC_ENU.b2(:,ii-1),1),1);
        
        RS_ENU.Sxx(:,ii) = NaN(size(RS_ENU.Sxx(:,ii-1),1),1);
        RS_ENU.Syy(:,ii) = NaN(size(RS_ENU.Syy(:,ii-1),1),1);
        RS_ENU.Sxy(:,ii) = NaN(size(RS_ENU.Sxy(:,ii-1),1),1);
        RS_ENU.Sxx_ss(ii) = NaN;
        RS_ENU.Syy_ss(ii) = NaN;
        RS_ENU.Sxy_ss(ii) = NaN;
        
        Tp_ENU(ii) = NaN;
        
        Eflux_ENU.Epos(:,ii) = NaN(size(Eflux_ENU.Epos(:,ii-1),1),1);
        Eflux_ENU.Eneg(:,ii) = NaN(size(Eflux_ENU.Eneg(:,ii-1),1),1);
        Eflux_ENU.Fpos(:,ii) = NaN(size(Eflux_ENU.Fpos(:,ii-1),1),1);
        Eflux_ENU.Fneg(:,ii) = NaN(size(Eflux_ENU.Fneg(:,ii-1),1),1);
        Eflux_ENU.depth(ii) = NaN;
    end    
end
i_ig = PUV_process_ENU(1).ids.i_ig;
i_swell = PUV_process_ENU(1).ids.i_swell;
fm_ENU = PUV_process_ENU(1).Spec.fm;
df_ENU = fm_ENU(2)-fm_ENU(1);
toc
R2_ig_ENU = NaN(1,segtotal);
R2_ss_ENU = NaN(1,segtotal);
for ii = 1:segtotal
    if ~isempty(PUV_process_ENU(ii).Eflux)
        
        R2_ig_ENU(ii) = sum(Eflux_ENU.Fneg(i_ig,ii).*df)/sum(Eflux_ENU.Fpos(i_ig,ii).*df);
        R2_ss_ENU(ii) = sum(Eflux_ENU.Fneg(i_swell,ii).*df)/sum(Eflux_ENU.Fpos(i_swell,ii).*df);
    
    end
    if R2_ig_ENU(ii)==0; R2_ig_ENU(ii)=NaN;end
    if R2_ss_ENU(ii)==0; R2_ss_ENU(ii)=NaN;end
end
toc
%% Extract values MOP FOR
disp('Extracting values from PUV_process_MOP structure')
tic
for ii = 1:segtotal
    if ~isempty(PUV_process_MOP(ii).ztest)
        
        Spec_MOP.SSE(:,ii) = PUV_process_MOP(ii).Spec.SSE;
        Spec_MOP.Spp(:,ii) = PUV_process_MOP(ii).Spec.Spp;
        Spec_MOP.Suu(:,ii) = PUV_process_MOP(ii).Spec.Suu;
        Spec_MOP.Svv(:,ii) = PUV_process_MOP(ii).Spec.Svv;
        Spec_MOP.Suv(:,ii) = PUV_process_MOP(ii).Spec.Suv;
        
        Ztest_MOP.ztest_all(:,ii) = PUV_process_MOP(ii).ztest.ztest_all;
        Ztest_MOP.ztest_ss(:,ii) = PUV_process_MOP(ii).ztest.ztest_ss;
        Ztest_MOP.ztest_ig(:,ii) = PUV_process_MOP(ii).ztest.ztest_ig;
        Ztest_MOP.ztest_ss_sum(ii) = PUV_process_MOP(ii).ztest.ztest_ss_sum;
        Ztest_MOP.ztest_ig_sum(ii) = PUV_process_MOP(ii).ztest.ztest_ig_sum;
        
        Hsig_MOP.Hs(ii) = PUV_process_MOP(ii).Hsig.Hs;
        Hsig_MOP.Hsss(ii) = PUV_process_MOP(ii).Hsig.Hs_ss;
        Hsig_MOP.Hsig(ii) = PUV_process_MOP(ii).Hsig.Hs_ig;
        
        Dir_MOP.dir1(:,ii) = PUV_process_MOP(ii).dir.dir1;
        Dir_MOP.spread1(:,ii) = PUV_process_MOP(ii).dir.spread1;
        Dir_MOP.dir2(:,ii) = PUV_process_MOP(ii).dir.dir2;
        Dir_MOP.spread2(:,ii) = PUV_process_MOP(ii).dir.spread2;
        Dir_MOP.dir1_ss_sum(ii) = PUV_process_MOP(ii).dir.dir1_ss_sum;
        Dir_MOP.spread1_ss_sum(ii) = PUV_process_MOP(ii).dir.spread1_ss_sum;
        Dir_MOP.dir2_ss_sum(ii) = PUV_process_MOP(ii).dir.dir2_ss_sum;
        Dir_MOP.spread2_ss_sum(ii) = PUV_process_MOP(ii).dir.spread2_ss_sum;
        
        FC_MOP.a1(:,ii) = PUV_process_MOP(ii).FC.a1;
        FC_MOP.b1(:,ii) = PUV_process_MOP(ii).FC.b1;
        FC_MOP.a2(:,ii) = PUV_process_MOP(ii).FC.a2;
        FC_MOP.b2(:,ii) = PUV_process_MOP(ii).FC.b2;
        
        RS_MOP.Sxx(:,ii) = PUV_process_MOP(ii).RS.Sxx;
        RS_MOP.Syy(:,ii) = PUV_process_MOP(ii).RS.Syy;
        RS_MOP.Sxy(:,ii) = PUV_process_MOP(ii).RS.Sxy;
        RS_MOP.Sxx_ss(ii) = PUV_process_MOP(ii).RS.Sxx_ss;
        RS_MOP.Syy_ss(ii) = PUV_process_MOP(ii).RS.Syy_ss;
        RS_MOP.Sxy_ss(ii) = PUV_process_MOP(ii).RS.Sxy_ss;
        
        Tp_MOP(ii) = 1./PUV_process_MOP(ii).ids.fpeak;
        
        Eflux_MOP.Fpos(:,ii) = PUV_process_MOP(ii).Eflux.posX;
        Eflux_MOP.Fneg(:,ii) = PUV_process_MOP(ii).Eflux.negX;
        Eflux_MOP.Epos(:,ii) = PUV_process_MOP(ii).Eflux.posX./PUV_process(ii).Eflux.Cg;
        Eflux_MOP.Eneg(:,ii) = PUV_process_MOP(ii).Eflux.negX./PUV_process(ii).Eflux.Cg;
        Eflux_MOP.depth(ii) = PUV_process_MOP(ii).Eflux.depth;
    else
        
        Spec_MOP.SSE(:,ii) = NaN(size(Spec_MOP.SSE(:,ii-1),1),1);
        Spec_MOP.Spp(:,ii) = NaN(size(Spec_MOP.Spp(:,ii-1),1),1);
        Spec_MOP.Suu(:,ii) = NaN(size(Spec_MOP.Suu(:,ii-1),1),1);
        Spec_MOP.Svv(:,ii) = NaN(size(Spec_MOP.Svv(:,ii-1),1),1);
        Spec_MOP.Suv(:,ii) = NaN(size(Spec_MOP.Suv(:,ii-1),1),1);
        
        Ztest_MOP.ztest_all(:,ii) = NaN(size(Ztest_MOP.ztest_all(:,ii-1),1),1);
        Ztest_MOP.ztest_ss(:,ii) = NaN(size(Ztest_MOP.ztest_ss(:,ii-1),1),1);
        Ztest_MOP.ztest_ig(:,ii) = NaN(size(Ztest_MOP.ztest_ig(:,ii-1),1),1);
        Ztest_MOP.ztest_ss_sum(ii) = NaN;
        Ztest_MOP.ztest_ig_sum(ii) = NaN;
        
        Hsig_MOP.Hs(ii) = NaN;
        Hsig_MOP.Hsss(ii) = NaN;
        Hsig_MOP.Hsig(ii) = NaN;
        
        Dir_MOP.dir1(:,ii) = NaN(size(Dir_MOP.dir1(:,ii-1),1),1);
        Dir_MOP.spread1(:,ii) = NaN(size(Dir_MOP.spread1(:,ii-1),1),1);
        Dir_MOP.dir2(:,ii) = NaN(size(Dir_MOP.dir2(:,ii-1),1),1);
        Dir_MOP.spread2(:,ii) = NaN(size(Dir_MOP.spread2(:,ii-1),1),1);
        Dir_MOP.dir1_ss_sum(ii) = NaN;
        Dir_MOP.spread1_ss_sum(ii) = NaN;
        Dir_MOP.dir2_ss_sum(ii) = NaN;
        Dir_MOP.spread2_ss_sum(ii) = NaN;
                
        FC_MOP.a1(:,ii) = NaN(size(FC_MOP.a1(:,ii-1),1),1);
        FC_MOP.b1(:,ii) = NaN(size(FC_MOP.b1(:,ii-1),1),1);
        FC_MOP.a2(:,ii) = NaN(size(FC_MOP.a2(:,ii-1),1),1);
        FC_MOP.b2(:,ii) = NaN(size(FC_MOP.b2(:,ii-1),1),1);
        
        RS_MOP.Sxx(:,ii) = NaN(size(RS_MOP.Sxx(:,ii-1),1),1);
        RS_MOP.Syy(:,ii) = NaN(size(RS_MOP.Syy(:,ii-1),1),1);
        RS_MOP.Sxy(:,ii) = NaN(size(RS_MOP.Sxy(:,ii-1),1),1);
        RS_MOP.Sxx_ss(ii) = NaN;
        RS_MOP.Syy_ss(ii) = NaN;
        RS_MOP.Sxy_ss(ii) = NaN;
        
        Tp_MOP(ii) = NaN;
        
        Eflux_MOP.Epos(:,ii) = NaN(size(Eflux_MOP.Epos(:,ii-1),1),1);
        Eflux_MOP.Eneg(:,ii) = NaN(size(Eflux_MOP.Eneg(:,ii-1),1),1);
        Eflux_MOP.Fpos(:,ii) = NaN(size(Eflux_MOP.Fpos(:,ii-1),1),1);
        Eflux_MOP.Fneg(:,ii) = NaN(size(Eflux_MOP.Fneg(:,ii-1),1),1);
        Eflux_MOP.depth(ii) = NaN;
    end    
end
i_ig = PUV_process_MOP(1).ids.i_ig;
i_swell = PUV_process_MOP(1).ids.i_swell;
fm_MOP = PUV_process_MOP(1).Spec.fm;
df_MOP = fm_MOP(2)-fm_MOP(1);
toc


toc
disp('Making QA/QC plots')
%% Ztest check
figure(1);clf
subplot(211)
plot(time(1,:), Ztest.ztest_ss_sum)
hold on
title('Ztest SS')
set(gca, 'FontSize', 20)
grid on
subplot(212)
plot(time(1,:), Ztest.ztest_ig_sum)
title('Ztest IG')
set(gca, 'FontSize', 20)
grid on
hline(1.2)
hline(0.8)
set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'Ztest.png')
%% Mean direction (f)
id = 112; disp(['Selected comparison time is ' datestr(PUV_work.time(1,id)) ' , Hs = ' num2str(Hsig.Hs(id),'%2.2f') 'm'])

figure(2);clf
subplot(121)
polarplot(deg2rad(Dir.dir1(i_ig,id)), fm(i_ig), '.')
hold on
polarplot(deg2rad(Dir.dir2(i_ig,id)), fm(i_ig), '.')
ax = gca;
ax.ThetaDir = 'counterclockwise';
title('IG')
subplot(122)
polarplot(deg2rad(Dir.dir1(i_swell,id)), fm(i_swell), '.')
hold on
ax=gca;
polarplot(deg2rad(Dir.dir2(i_swell,id)), fm(i_swell), '.')
ax.ThetaDir = 'counterclockwise';
ax.RLim=[0 0.25];
ax.RTick = [0 0.05 0.1 0.15 0.2 0.25];
title('SS')
%legend('dir1', 'dir2')
sgtitle([string(time(1,id)) '- SN Coordinates'])

set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,['Dir_polar_id' num2str(id) '_SN_CCW.png'])
%% Reflection Coeff
figure(3);clf
subplot(311)
plot(time(1,:), Hsig.Hsss)
hold on
plot(time(1,:), Hsig.Hsig)
subplot(312)
plot(time(1,:), R2_ss)
subplot(313)
plot(time(1),1)
hold on
plot(time(1,:), R2_ig)
ylabel('R^2'); set(gca, 'FontSize', 20);title('R^2_{IG}')
subplot(311); ylabel('Hsig m'); set(gca, 'FontSize', 20);title('Hsig')
subplot(312); ylabel('R^2'); set(gca, 'FontSize', 20); title('R^2_{SS}')
subplot(313);grid on
subplot(311);grid on
subplot(312);grid on
set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'Reflection_coeff_vHs.png')
%% Mean Direction
figure(4);clf

subplot(211)
plot(time(1,:), Dir.dir1_ss_sum)
hold on
plot(time(1,:), Dir.dir2_ss_sum)
title('Mean direction SS')
legend('Mean direction (b1/a1)', 'Peak direction (b2/a2)')
set(gca, 'FontSize', 20)
grid on
subplot(212)
plot(time(1,:), Dir.spread1_ss_sum)
hold on
plot(time(1,:), Dir.spread2_ss_sum)
title('Mean spread SS')
legend('Mean spread 1st', 'Mean spread 2nd')

set(gca, 'FontSize', 20)
grid on

set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'Mean_dir_spread.png')
%% Phase and Coherence
for ii = 1:segtotal
if ~isempty(PUV_process(ii).ztest)
ph(:,ii)=PUV_process(ii).coh.phPUpres;
end
end
ph(ph==0)=NaN;
for ii = 1:segtotal
if ~isempty(PUV_process(ii).ztest)
coh(:,ii)=PUV_process(ii).coh.cohPUpres;
end
end
coh(coh==0)=NaN;
figure(5);clf
subplot(211)
h=pcolor(time(1,:)',fm, ph(1:901,:))
c =colorbar('location','WestOutside');
set(h, 'EdgeColor', 'none')
c.Label.String = 'Phase'
set(gca, 'FontSize', 20)
caxis([-180 180])
subplot(212)
h=pcolor(time(1,:)',fm, coh(1:901,:))
c =colorbar('location','WestOutside');
set(h, 'EdgeColor', 'none')
c.Label.String = 'Phase'
set(gca, 'FontSize', 20)
c.Label.String = 'Coherence'
sgtitle('PU phase and coherence')

set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'PU_phase_coherence.png')
%%
MOP = read_MOPline('D0582',time(1,1)-hours(1),time(end,end));

% find(min(abs(time(1)-timeseries_mop)) == abs(time(1)-timeseries_mop))


%% Compare bulk parameters to MOP
figure(6);clf
subplot(411)
% plot(timeseries_mop, Hs_mop)
plot(MOP.time,MOP.Hs)
hold on
plot(time(1,:), Hsig.Hsss)
%plot(time(1,:), Hsig_MOP.Hsss)
legend('MOP', 'PUV', 'PUV MOP FOR')
xlim([time(1) time(end)])
title('Hs')
ylabel('m')
set(gca, 'FontSize', 20)

subplot(412)
% plot(timeseries_mop, Tp_mop)
plot(MOP.time,1./MOP.fp)
hold on
plot(time(1,:), Tp)
%plot(time(1,:), Tp_MOP)
legend('MOP', 'PUV', 'PUV MOP FOR')
xlim([time(1) time(end)])
title('Tp')
ylabel('s')
set(gca, 'FontSize', 20)

subplot(413)
% plot(timeseries_mop, Sxx_mop)
plot(MOP.time,MOP.Sxx)
hold on
plot(time(1,:), RS.Sxx_ss)
%plot(time(1,:), RS_MOP.Sxx_ss)
legend('MOP', 'PUV', 'PUV MOP FOR')
xlim([time(1) time(end)])
title('Sxx')
ylabel('m^2')
set(gca, 'FontSize', 20)

subplot(414)
% plot(timeseries_mop, Sxy_mop)
plot(MOP.time,MOP.Sxy)
hold on
plot(time(1,:), RS.Sxy_ss)
%plot(time(1,:), RS_MOP.Sxy_ss)
legend('MOP', 'PUV', 'PUV MOP FOR')
xlim([time(1) time(end)])
title('Sxy')
ylabel('m^2')
set(gca, 'FontSize', 20)

set(gcf, 'Position', [100,100,1500,1300])
saveas(gcf,'MOPvPUV_bulk.png')
%%
figure(7);clf
subplot(211)
scatter(Eflux.depth, R2_ss, 50, Dir.spread1_ss_sum,'filled')
ylabel('R^2_{SS}')
set(gca, 'FontSize', 20)
grid on
cb = colorbar;
cb.Label.String = 'Spread1 (deg)';
xlabel('<-- Low tide --- Water Depth (m) --- High tide -->')
title('SS')
subplot(212)
scatter(Eflux.depth, R2_ig, 'filled')
ylabel('R^2_{IG}')
set(gca, 'FontSize', 20)
grid on
title('IG')
xlabel('<-- Low tide --- Water Depth (m) --- High tide -->')

set(gcf, 'Position', [100,100,1000,800])
saveas(gcf,'R2_vtide_cspread1.png')
%%
figure(8);clf
subplot(211)
plot(time(1,:),sum(Eflux.Epos(i_ig,:),1)./sum(Eflux.Epos(i_swell,:),1))
hold on
ylim([-0.015 0.03])
ylabel('E^+_{IG}/ E^+_{SS}');set(gca, 'FontSize', 20)
yyaxis right
plot(time(1,:),Eflux.depth)
ylim([10 13])
xlim([time(1,10) time(1,400)])
ylabel('Water depth (m)') 
subplot(212)
scatter(sum(Eflux.Epos(i_swell,:),1),sum(Eflux.Epos(i_ig,:),1), 50, Eflux.depth, 'filled')
cb = colorbar;
cb.Label.String = 'Water depth (m)';
xlabel('E^+_{SS}')
ylabel('E^+{IG}')
set(gca, 'FontSize', 20)

set(gcf, 'Position', [100,100,1000,800])
saveas(gcf,'EigEss_v_depth.png')
%% Directional Spreading Function
clear d ds
%compute directional spreding function using MEM estimator
for ii = id
    
    d(:,:,ii)= mem_est(FC.a1(:,ii)', FC.a2(:,ii)', FC.b1(:,ii)', FC.b2(:,ii)');
    if isempty(PUV_process(ii).Spec)
        PUV_process(ii).Spec.SSE = NaN(length(fm),1);
    end
    for i=1:length(fm) % loop through freq bands
        ds(i,:,ii)=d(i,:,ii)*PUV_process(ii).Spec.SSE(i); % mutiply by the freq band total energy
    end

end
ds(:,361,:)=ds(:,1,:);

%% Directional Spreading Function
clear d_ENU ds_ENU
%compute directional spreding function using MEM estimator
for ii = 1:segtotal
    
    d_ENU(:,:,ii)= getmem(FC_ENU.a1(:,ii)', FC_ENU.a2(:,ii)', FC_ENU.b1(:,ii)', FC_ENU.b2(:,ii)');
    if isempty(PUV_process_ENU(ii).Spec)
        PUV_process_ENU(ii).Spec.SSE = NaN(length(fm),1);
    end
    for i=1:length(fm) % loop through freq bands
        ds_ENU(i,:,ii)=d_ENU(i,:,ii)*PUV_process_ENU(ii).Spec.SSE(i); % mutiply by the freq band total energy
    end

end
%% Directional Spreading Function MOP FOR
clear d_MOP ds_MOP
%compute directional spreding function using MEM estimator
for ii = 1:segtotal
    
    d_MOP(:,:,ii)= mem_est(FC_MOP.a1(:,ii), FC_MOP.a2(:,ii), FC_MOP.b1(:,ii), FC_MOP.b2(:,ii));
    if isempty(PUV_process_MOP(ii).Spec)
        PUV_process_MOP(ii).Spec.SSE = NaN(length(fm_MOP),1);
    end
    for i=1:length(fm_MOP) % loop through freq bands
        ds_MOP(i,:,ii)=d_MOP(i,:,ii)*PUV_process_MOP(ii).Spec.SSE(i); % mutiply by the freq band total energy
    end

end
%% Checks on directional spectra
aa = ds(:,1:360,id);
ab = sum(aa,2);

figure(9);clf
plot(fm,ab, 'LineWidth', 5)
hold on
plot(fm,PUV_process(id).Spec.SSE, 'LineWidth', 3)
title('E(f) = \int_0^{2\pi}E(f,\theta)d\theta')
legend('\int_0^{2\pi}E(f,\theta)d\theta', 'E(f)')
set(gca, 'FontSize',  20)

set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'Ef_vEftheta.png')

figure(10);clf
plot(fm, sum(d(:,:,id),2))
title('\int_0^{2\pi}G(\theta,f)d\theta = 1')
set(gca, 'FontSize', 20)

set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'intG_check.png')

% id = find(min(abs(time(1)-timeseries_mop)) == abs(time(1)-timeseries_mop))
% idmop = 1;
a1 = FC.a1(:,id);
b1 = FC.b1(:,id);
a2 = FC.a2(:,id);
b2 = FC.b2(:,id);
%%
figure(11);clf
theta = 0:359;
subplot(221)
ab = sum(aa.*cosd(theta),2);
plot(fm,ab)
hold on
plot(fm, a1)
plot(MOP.frequency, MOP.a1(id,:))
legend('\int_0^{2\pi}G(f,\theta)cos{\theta}d\theta', 'a1', 'mop')
title('a1')
set(gca, 'FontSize', 20)

subplot(222)
ab = sum(aa.*cosd(2*theta),2);
plot(fm,ab)
hold on
plot(fm, a2)
plot(MOP.frequency, MOP.a2(id,:))
legend('\int_0^{2\pi}G(f,\theta)cos{2\theta}d\theta', 'a2', 'mop')
title('a2')
set(gca, 'FontSize', 20)

subplot(223)
ab = sum(aa.*sind(theta),2);
plot(fm,ab)
hold on
plot(fm, b1)
plot(MOP.frequency, MOP.b1(id,:))
legend('\int_0^{2\pi}G(f,\theta)sin{\theta}d\theta', 'b1', 'mop')
title('b1')
set(gca, 'FontSize', 20)

subplot(224)
ab = sum(aa.*sind(2*theta),2);
plot(fm,ab)
hold on
plot(fm, b2)
plot(MOP.frequency, MOP.b2(id,:))
legend('\int_0^{2\pi}G(f,\theta)sin{2\theta}d\theta', 'b2', 'mop')
title('b2')
set(gca, 'FontSize', 20)

% set(gcf, 'Position', [100,100,2000,1000])
saveas(gcf,'mop_puv_dirspec_FC.png')
%% Plot Dirspectra
ds(:,361,:)=ds(:,1,:);
figure(12);clf
polarPcolor(fm', [0:360], ds(:,:,id), 'Nspokes',13,'typeRose','default')

%h1 = pcolor(fm', [0:360], log(ds(:,:,83))')
sgtitle('Dirspec from PUV in Shorenormal Coords')
%set(h1, 'EdgeColor', 'none')

set(gcf, 'Position', [100,100,600,600])
saveas(gcf,['Dirspec_id' num2str(id) '_SN.png'])
%%
%ds_MOP(:,361,:)=ds_MOP(:,1,:);
% figure(2);clf
% polarPcolor(fm', [0:360], (log(ds_MOP(:,:,83))'), 'XAngle', 270-shorenormal+180,'Nspokes',13)
% %h2=pcolor(fm', [0:360], log(ds_MOP(:,:,83))')
% sgtitle('Dirspec from PUV in MOP Coords')
% caxis([-15 -5])
%set(h2, 'EdgeColor', 'none')

ds_ENU(:,361,:)=ds_ENU(:,1,:);
figure(13);clf
polarPcolor(fm',  [0:360], log(ds_ENU(:,:,id))', 'XAngle',0,'Nspokes',13)
%h3 = pcolor(fm', [0:360], log(ds_ENU(:,:,83))')
%set(h3, 'EdgeColor', 'none')
sgtitle('Dirspec from PUV in ENU Coords')
% caxis([-15 -5])
set(gcf, 'Position', [100,100,600,600])
saveas(gcf,['Dirspec_id' num2str(id) '_ENU.png'])
%%
% id = find(min(abs(time(83)-timeseries_mop)) == abs(time(83)-timeseries_mop))
clear ds_mopmop d_mopmop
    d_mopmop=getmem(MOP.a1(id,:), MOP.a2(id,:), MOP.b1(id,:), MOP.b2(id,:));
    
    for i=1:length(MOP.fbw) % loop through freq bands
        ds_mopmop(i,:)=d_mopmop(i,:)*MOP.spec1D(id,i)'; % mutiply by the freq band total energy
    end

%%
ds_mopmop(:,361)=ds_mopmop(:,1);
figure(14);clf
% polarPcolor(double(Fq_mop'), [0:360], log(ds_mopmop)', 'XAngle', 270+180,'Nspokes',13)
polarPcolor(double(MOP.frequency)', [0:360], ds_mopmop', 'Nspokes',13,'typeRose','meteo')

sgtitle('Dirspec from MOP')
% caxis([-15 -5])
