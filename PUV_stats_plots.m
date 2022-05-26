% %% 
% 
% %% Initial QC
% % PUV_1QC

%% load variables up to this point
load temp

%% Extract values - Shorenormal Coordinates
disp('Extracting values from PUV_process structure')
tic
segtotal = length(PUV_process);

for ii = 1:segtotal
    if ~isempty(PUV_process(ii).ztest)
        ii
        time(:,ii) = PUV_process(ii).time;
        
        Spec_sn.SSE(:,ii) = PUV_process(ii).Spec.SSE;
        Spec_sn.Spp(:,ii) = PUV_process(ii).Spec.Spp;
        Spec_sn.Suu(:,ii) = PUV_process(ii).Spec.Suu;
        Spec_sn.Svv(:,ii) = PUV_process(ii).Spec.Svv;
        Spec_sn.Suv(:,ii) = PUV_process(ii).Spec.Suv;
        
        Ztest_sn.ztest_all(:,ii) = PUV_process(ii).ztest.ztest_all;
        Ztest_sn.ztest_ss(:,ii) = PUV_process(ii).ztest.ztest_ss;
        Ztest_sn.ztest_ig(:,ii) = PUV_process(ii).ztest.ztest_ig;
        Ztest_sn.ztest_ss_sum(ii) = PUV_process(ii).ztest.ztest_ss_sum;
        Ztest_sn.ztest_ig_sum(ii) = PUV_process(ii).ztest.ztest_ig_sum;
        
        Hsig_sn.Hs(ii) = PUV_process(ii).Hsig.Hs;
        Hsig_sn.Hsss(ii) = PUV_process(ii).Hsig.Hs_ss;
        Hsig_sn.Hsig(ii) = PUV_process(ii).Hsig.Hs_ig;
        
        Dir_sn.dir1(:,ii) = PUV_process(ii).dir.dir1;
        Dir_sn.spread1(:,ii) = PUV_process(ii).dir.spread1;
        Dir_sn.dir2(:,ii) = PUV_process(ii).dir.dir2;
        Dir_sn.spread2(:,ii) = PUV_process(ii).dir.spread2;
        Dir_sn.dir1_ss_sum(ii) = PUV_process(ii).dir.dir1_ss_sum;
        Dir_sn.spread1_ss_sum(ii) = PUV_process(ii).dir.spread1_ss_sum;
        Dir_sn.dir2_ss_sum(ii) = PUV_process(ii).dir.dir2_ss_sum;
        Dir_sn.spread2_ss_sum(ii) = PUV_process(ii).dir.spread2_ss_sum;
        
        FC_sn.a1(:,ii) = PUV_process(ii).FC.a1;
        FC_sn.b1(:,ii) = PUV_process(ii).FC.b1;
        FC_sn.a2(:,ii) = PUV_process(ii).FC.a2;
        FC_sn.b2(:,ii) = PUV_process(ii).FC.b2;
        
        RS_sn.Sxx(:,ii) = PUV_process(ii).RS.Sxx;
        RS_sn.Syy(:,ii) = PUV_process(ii).RS.Syy;
        RS_sn.Sxy(:,ii) = PUV_process(ii).RS.Sxy;
        RS_sn.Sxx_ss(ii) = PUV_process(ii).RS.Sxx_ss;
        RS_sn.Syy_ss(ii) = PUV_process(ii).RS.Syy_ss;
        RS_sn.Sxy_ss(ii) = PUV_process(ii).RS.Sxy_ss;
        
        Tp_sn(ii) = 1./PUV_process(ii).ids.fpeak;
        
        Eflux_sn.Fpos(:,ii) = PUV_process(ii).Eflux.posX;
        Eflux_sn.Fneg(:,ii) = PUV_process(ii).Eflux.negX;
        Eflux_sn.Epos(:,ii) = PUV_process(ii).Eflux.posX./PUV_process(ii).Eflux.Cg;
        Eflux_sn.Eneg(:,ii) = PUV_process(ii).Eflux.negX./PUV_process(ii).Eflux.Cg;
        Eflux_sn.depth(ii) = PUV_process(ii).Eflux.depth;
    else

        time(:,ii) = NaT(7200,1);
        
        Spec_sn.SSE(:,ii) = NaN(size(Spec_sn.SSE(:,ii-1),1),1);
        Spec_sn.Spp(:,ii) = NaN(size(Spec_sn.Spp(:,ii-1),1),1);
        Spec_sn.Suu(:,ii) = NaN(size(Spec_sn.Suu(:,ii-1),1),1);
        Spec_sn.Svv(:,ii) = NaN(size(Spec_sn.Svv(:,ii-1),1),1);
        Spec_sn.Suv(:,ii) = NaN(size(Spec_sn.Suv(:,ii-1),1),1);
        
        Ztest_sn.ztest_all(:,ii) = NaN(size(Ztest_sn.ztest_all(:,ii-1),1),1);
        Ztest_sn.ztest_ss(:,ii) = NaN(size(Ztest_sn.ztest_ss(:,ii-1),1),1);
        Ztest_sn.ztest_ig(:,ii) = NaN(size(Ztest_sn.ztest_ig(:,ii-1),1),1);
        Ztest_sn.ztest_ss_sum(ii) = NaN;
        Ztest_sn.ztest_ig_sum(ii) = NaN;
        
        Hsig_sn.Hs(ii) = NaN;
        Hsig_sn.Hsss(ii) = NaN;
        Hsig_sn.Hsig(ii) = NaN;
        
        Dir_sn.dir1(:,ii) = NaN(size(Dir_sn.dir1(:,ii-1),1),1);
        Dir_sn.spread1(:,ii) = NaN(size(Dir_sn.spread1(:,ii-1),1),1);
        Dir_sn.dir2(:,ii) = NaN(size(Dir_sn.dir2(:,ii-1),1),1);
        Dir_sn.spread2(:,ii) = NaN(size(Dir_sn.spread2(:,ii-1),1),1);
        Dir_sn.dir1_ss_sum(ii) = NaN;
        Dir_sn.spread1_ss_sum(ii) = NaN;
        Dir_sn.dir2_ss_sum(ii) = NaN;
        Dir_sn.spread2_ss_sum(ii) = NaN;
                
        FC_sn.a1(:,ii) = NaN(size(FC_sn.a1(:,ii-1),1),1);
        FC_sn.b1(:,ii) = NaN(size(FC_sn.b1(:,ii-1),1),1);
        FC_sn.a2(:,ii) = NaN(size(FC_sn.a2(:,ii-1),1),1);
        FC_sn.b2(:,ii) = NaN(size(FC_sn.b2(:,ii-1),1),1);
        
        RS_sn.Sxx(:,ii) = NaN(size(RS_sn.Sxx(:,ii-1),1),1);
        RS_sn.Syy(:,ii) = NaN(size(RS_sn.Syy(:,ii-1),1),1);
        RS_sn.Sxy(:,ii) = NaN(size(RS_sn.Sxy(:,ii-1),1),1);
        RS_sn.Sxx_ss(ii) = NaN;
        RS_sn.Syy_ss(ii) = NaN;
        RS_sn.Sxy_ss(ii) = NaN;
        
        Tp_sn(ii) = NaN;
        
        Eflux_sn.Epos(:,ii) = NaN(size(Eflux_sn.Epos(:,ii-1),1),1);
        Eflux_sn.Eneg(:,ii) = NaN(size(Eflux_sn.Eneg(:,ii-1),1),1);
        Eflux_sn.Fpos(:,ii) = NaN(size(Eflux_sn.Fpos(:,ii-1),1),1);
        Eflux_sn.Fneg(:,ii) = NaN(size(Eflux_sn.Fneg(:,ii-1),1),1);
        Eflux_sn.depth(ii) = NaN;
    end    
end

i_ig_sn = PUV_process(1).ids.i_ig;
i_swell_sn = PUV_process(1).ids.i_swell;
fm_sn = PUV_process(1).Spec.fm;
df_sn = fm_sn(2)-fm_sn(1);
toc
disp('Computing reflection coefficient')
R2_ig_sn = NaN(1,segtotal);
R2_ss_sn = NaN(1,segtotal);
for ii = 1:segtotal
    if ~isempty(PUV_process(ii).Eflux)
        
        R2_ig_sn(ii) = sum(Eflux_sn.Fneg(i_ig_sn,ii).*df_sn)/sum(Eflux_sn.Fpos(i_ig_sn,ii).*df_sn);
        R2_ss_sn(ii) = sum(Eflux_sn.Fneg(i_swell_sn,ii).*df_sn)/sum(Eflux_sn.Fpos(i_swell_sn,ii).*df_sn);
    
    end
    if R2_ig_sn(ii)==0; R2_ig_sn(ii)=NaN;end
    if R2_ss_sn(ii)==0; R2_ss_sn(ii)=NaN;end
end
toc

%% Extract values ENU FOR
% disp('Extracting values from PUV_process_ENU structure')
% tic
% for ii = 1:segtotal
%     if ~isempty(PUV_process_ENU(ii).ztest)
%         
%         Spec_ENU.SSE(:,ii) = PUV_process_ENU(ii).Spec.SSE;
%         Spec_ENU.Spp(:,ii) = PUV_process_ENU(ii).Spec.Spp;
%         Spec_ENU.Suu(:,ii) = PUV_process_ENU(ii).Spec.Suu;
%         Spec_ENU.Svv(:,ii) = PUV_process_ENU(ii).Spec.Svv;
%         Spec_ENU.Suv(:,ii) = PUV_process_ENU(ii).Spec.Suv;
%         
%         Ztest_ENU.ztest_all(:,ii) = PUV_process_ENU(ii).ztest.ztest_all;
%         Ztest_ENU.ztest_ss(:,ii) = PUV_process_ENU(ii).ztest.ztest_ss;
%         Ztest_ENU.ztest_ig(:,ii) = PUV_process_ENU(ii).ztest.ztest_ig;
%         Ztest_ENU.ztest_ss_sum(ii) = PUV_process_ENU(ii).ztest.ztest_ss_sum;
%         Ztest_ENU.ztest_ig_sum(ii) = PUV_process_ENU(ii).ztest.ztest_ig_sum;
%         
%         Hsig_ENU.Hs(ii) = PUV_process_ENU(ii).Hsig.Hs;
%         Hsig_ENU.Hsss(ii) = PUV_process_ENU(ii).Hsig.Hs_ss;
%         Hsig_ENU.Hsig(ii) = PUV_process_ENU(ii).Hsig.Hs_ig;
%         
%         Dir_ENU.dir1(:,ii) = PUV_process_ENU(ii).dir.dir1;
%         Dir_ENU.spread1(:,ii) = PUV_process_ENU(ii).dir.spread1;
%         Dir_ENU.dir2(:,ii) = PUV_process_ENU(ii).dir.dir2;
%         Dir_ENU.spread2(:,ii) = PUV_process_ENU(ii).dir.spread2;
%         Dir_ENU.dir1_ss_sum(ii) = PUV_process_ENU(ii).dir.dir1_ss_sum;
%         Dir_ENU.spread1_ss_sum(ii) = PUV_process_ENU(ii).dir.spread1_ss_sum;
%         Dir_ENU.dir2_ss_sum(ii) = PUV_process_ENU(ii).dir.dir2_ss_sum;
%         Dir_ENU.spread2_ss_sum(ii) = PUV_process_ENU(ii).dir.spread2_ss_sum;
%         
%         FC_ENU.a1(:,ii) = PUV_process_ENU(ii).FC.a1;
%         FC_ENU.b1(:,ii) = PUV_process_ENU(ii).FC.b1;
%         FC_ENU.a2(:,ii) = PUV_process_ENU(ii).FC.a2;
%         FC_ENU.b2(:,ii) = PUV_process_ENU(ii).FC.b2;
%         
%         RS_ENU.Sxx(:,ii) = PUV_process_ENU(ii).RS.Sxx;
%         RS_ENU.Syy(:,ii) = PUV_process_ENU(ii).RS.Syy;
%         RS_ENU.Sxy(:,ii) = PUV_process_ENU(ii).RS.Sxy;
%         RS_ENU.Sxx_ss(ii) = PUV_process_ENU(ii).RS.Sxx_ss;
%         RS_ENU.Syy_ss(ii) = PUV_process_ENU(ii).RS.Syy_ss;
%         RS_ENU.Sxy_ss(ii) = PUV_process_ENU(ii).RS.Sxy_ss;
%         
%         Tp_ENU(ii) = 1./PUV_process_ENU(ii).ids.fpeak;
%         
%         Eflux_ENU.Fpos(:,ii) = PUV_process_ENU(ii).Eflux.posX;
%         Eflux_ENU.Fneg(:,ii) = PUV_process_ENU(ii).Eflux.negX;
%         Eflux_ENU.Epos(:,ii) = PUV_process_ENU(ii).Eflux.posX./PUV_process_ENU(ii).Eflux.Cg;
%         Eflux_ENU.Eneg(:,ii) = PUV_process_ENU(ii).Eflux.negX./PUV_process_ENU(ii).Eflux.Cg;
%         Eflux_ENU.depth(ii) = PUV_process_ENU(ii).Eflux.depth;
%     else
%         
%         Spec_ENU.SSE(:,ii) = NaN(size(Spec_ENU.SSE(:,ii-1),1),1);
%         Spec_ENU.Spp(:,ii) = NaN(size(Spec_ENU.Spp(:,ii-1),1),1);
%         Spec_ENU.Suu(:,ii) = NaN(size(Spec_ENU.Suu(:,ii-1),1),1);
%         Spec_ENU.Svv(:,ii) = NaN(size(Spec_ENU.Svv(:,ii-1),1),1);
%         Spec_ENU.Suv(:,ii) = NaN(size(Spec_ENU.Suv(:,ii-1),1),1);
%         
%         Ztest_ENU.ztest_all(:,ii) = NaN(size(Ztest_ENU.ztest_all(:,ii-1),1),1);
%         Ztest_ENU.ztest_ss(:,ii) = NaN(size(Ztest_ENU.ztest_ss(:,ii-1),1),1);
%         Ztest_ENU.ztest_ig(:,ii) = NaN(size(Ztest_ENU.ztest_ig(:,ii-1),1),1);
%         Ztest_ENU.ztest_ss_sum(ii) = NaN;
%         Ztest_ENU.ztest_ig_sum(ii) = NaN;
%         
%         Hsig_ENU.Hs(ii) = NaN;
%         Hsig_ENU.Hsss(ii) = NaN;
%         Hsig_ENU.Hsig(ii) = NaN;
%         
%         Dir_ENU.dir1(:,ii) = NaN(size(Dir_ENU.dir1(:,ii-1),1),1);
%         Dir_ENU.spread1(:,ii) = NaN(size(Dir_ENU.spread1(:,ii-1),1),1);
%         Dir_ENU.dir2(:,ii) = NaN(size(Dir_ENU.dir2(:,ii-1),1),1);
%         Dir_ENU.spread2(:,ii) = NaN(size(Dir_ENU.spread2(:,ii-1),1),1);
%         Dir_ENU.dir1_ss_sum(ii) = NaN;
%         Dir_ENU.spread1_ss_sum(ii) = NaN;
%         Dir_ENU.dir2_ss_sum(ii) = NaN;
%         Dir_ENU.spread2_ss_sum(ii) = NaN;
%                 
%         FC_ENU.a1(:,ii) = NaN(size(FC_ENU.a1(:,ii-1),1),1);
%         FC_ENU.b1(:,ii) = NaN(size(FC_ENU.b1(:,ii-1),1),1);
%         FC_ENU.a2(:,ii) = NaN(size(FC_ENU.a2(:,ii-1),1),1);
%         FC_ENU.b2(:,ii) = NaN(size(FC_ENU.b2(:,ii-1),1),1);
%         
%         RS_ENU.Sxx(:,ii) = NaN(size(RS_ENU.Sxx(:,ii-1),1),1);
%         RS_ENU.Syy(:,ii) = NaN(size(RS_ENU.Syy(:,ii-1),1),1);
%         RS_ENU.Sxy(:,ii) = NaN(size(RS_ENU.Sxy(:,ii-1),1),1);
%         RS_ENU.Sxx_ss(ii) = NaN;
%         RS_ENU.Syy_ss(ii) = NaN;
%         RS_ENU.Sxy_ss(ii) = NaN;
%         
%         Tp_ENU(ii) = NaN;
%         
%         Eflux_ENU.Epos(:,ii) = NaN(size(Eflux_ENU.Epos(:,ii-1),1),1);
%         Eflux_ENU.Eneg(:,ii) = NaN(size(Eflux_ENU.Eneg(:,ii-1),1),1);
%         Eflux_ENU.Fpos(:,ii) = NaN(size(Eflux_ENU.Fpos(:,ii-1),1),1);
%         Eflux_ENU.Fneg(:,ii) = NaN(size(Eflux_ENU.Fneg(:,ii-1),1),1);
%         Eflux_ENU.depth(ii) = NaN;
%     end    
% end
% i_ig = PUV_process_ENU(1).ids.i_ig;
% i_swell = PUV_process_ENU(1).ids.i_swell;
% fm_ENU = PUV_process_ENU(1).Spec.fm;
% df_ENU = fm_ENU(2)-fm_ENU(1);
% toc
% R2_ig_ENU = NaN(1,segtotal);
% R2_ss_ENU = NaN(1,segtotal);
% for ii = 1:segtotal
%     if ~isempty(PUV_process_ENU(ii).Eflux)
%         
%         R2_ig_ENU(ii) = sum(Eflux_ENU.Fneg(i_ig,ii).*df)/sum(Eflux_ENU.Fpos(i_ig,ii).*df);
%         R2_ss_ENU(ii) = sum(Eflux_ENU.Fneg(i_swell,ii).*df)/sum(Eflux_ENU.Fpos(i_swell,ii).*df);
%     
%     end
%     if R2_ig_ENU(ii)==0; R2_ig_ENU(ii)=NaN;end
%     if R2_ss_ENU(ii)==0; R2_ss_ENU(ii)=NaN;end
% end
% toc
%% Extract values MOP FOR
disp('Extracting values from PUV_process_buoy structure')
tic
for ii = 1:segtotal
    if ~isempty(PUV_process_buoy(ii).ztest)
        
        Spec_buoy.SSE(:,ii) = PUV_process_buoy(ii).Spec.SSE;
        Spec_buoy.Spp(:,ii) = PUV_process_buoy(ii).Spec.Spp;
        Spec_buoy.Suu(:,ii) = PUV_process_buoy(ii).Spec.Suu;
        Spec_buoy.Svv(:,ii) = PUV_process_buoy(ii).Spec.Svv;
        Spec_buoy.Suv(:,ii) = PUV_process_buoy(ii).Spec.Suv;
        
        Ztest_buoy.ztest_all(:,ii) = PUV_process_buoy(ii).ztest.ztest_all;
        Ztest_buoy.ztest_ss(:,ii) = PUV_process_buoy(ii).ztest.ztest_ss;
        Ztest_buoy.ztest_ig(:,ii) = PUV_process_buoy(ii).ztest.ztest_ig;
        Ztest_buoy.ztest_ss_sum(ii) = PUV_process_buoy(ii).ztest.ztest_ss_sum;
        Ztest_buoy.ztest_ig_sum(ii) = PUV_process_buoy(ii).ztest.ztest_ig_sum;
        
        Hsig_buoy.Hs(ii) = PUV_process_buoy(ii).Hsig.Hs;
        Hsig_buoy.Hsss(ii) = PUV_process_buoy(ii).Hsig.Hs_ss;
        Hsig_buoy.Hsig(ii) = PUV_process_buoy(ii).Hsig.Hs_ig;
        
        Dir_buoy.dir1(:,ii) = PUV_process_buoy(ii).dir.dir1;
        Dir_buoy.spread1(:,ii) = PUV_process_buoy(ii).dir.spread1;
        Dir_buoy.dir2(:,ii) = PUV_process_buoy(ii).dir.dir2;
        Dir_buoy.spread2(:,ii) = PUV_process_buoy(ii).dir.spread2;
        Dir_buoy.dir1_ss_sum(ii) = PUV_process_buoy(ii).dir.dir1_ss_sum;
        Dir_buoy.spread1_ss_sum(ii) = PUV_process_buoy(ii).dir.spread1_ss_sum;
        Dir_buoy.dir2_ss_sum(ii) = PUV_process_buoy(ii).dir.dir2_ss_sum;
        Dir_buoy.spread2_ss_sum(ii) = PUV_process_buoy(ii).dir.spread2_ss_sum;
        
        FC_buoy.a1(:,ii) = PUV_process_buoy(ii).FC.a1;
        FC_buoy.b1(:,ii) = PUV_process_buoy(ii).FC.b1;
        FC_buoy.a2(:,ii) = PUV_process_buoy(ii).FC.a2;
        FC_buoy.b2(:,ii) = PUV_process_buoy(ii).FC.b2;
        
        RS_buoy.Sxx(:,ii) = PUV_process_buoy(ii).RS.Sxx;
        RS_buoy.Syy(:,ii) = PUV_process_buoy(ii).RS.Syy;
        RS_buoy.Sxy(:,ii) = PUV_process_buoy(ii).RS.Sxy;
        RS_buoy.Sxx_ss(ii) = PUV_process_buoy(ii).RS.Sxx_ss;
        RS_buoy.Syy_ss(ii) = PUV_process_buoy(ii).RS.Syy_ss;
        RS_buoy.Sxy_ss(ii) = PUV_process_buoy(ii).RS.Sxy_ss;
        
        Tp_buoy(ii) = 1./PUV_process_buoy(ii).ids.fpeak;
        
        Eflux_buoy.Fpos(:,ii) = PUV_process_buoy(ii).Eflux.posX;
        Eflux_buoy.Fneg(:,ii) = PUV_process_buoy(ii).Eflux.negX;
        Eflux_buoy.Epos(:,ii) = PUV_process_buoy(ii).Eflux.posX./PUV_process_buoy(ii).Eflux.Cg;
        Eflux_buoy.Eneg(:,ii) = PUV_process_buoy(ii).Eflux.negX./PUV_process_buoy(ii).Eflux.Cg;
        Eflux_buoy.depth(ii) = PUV_process_buoy(ii).Eflux.depth;
    else
        
        Spec_buoy.SSE(:,ii) = NaN(size(Spec_buoy.SSE(:,ii-1),1),1);
        Spec_buoy.Spp(:,ii) = NaN(size(Spec_buoy.Spp(:,ii-1),1),1);
        Spec_buoy.Suu(:,ii) = NaN(size(Spec_buoy.Suu(:,ii-1),1),1);
        Spec_buoy.Svv(:,ii) = NaN(size(Spec_buoy.Svv(:,ii-1),1),1);
        Spec_buoy.Suv(:,ii) = NaN(size(Spec_buoy.Suv(:,ii-1),1),1);
        
        Ztest_buoy.ztest_all(:,ii) = NaN(size(Ztest_buoy.ztest_all(:,ii-1),1),1);
        Ztest_buoy.ztest_ss(:,ii) = NaN(size(Ztest_buoy.ztest_ss(:,ii-1),1),1);
        Ztest_buoy.ztest_ig(:,ii) = NaN(size(Ztest_buoy.ztest_ig(:,ii-1),1),1);
        Ztest_buoy.ztest_ss_sum(ii) = NaN;
        Ztest_buoy.ztest_ig_sum(ii) = NaN;
        
        Hsig_buoy.Hs(ii) = NaN;
        Hsig_buoy.Hsss(ii) = NaN;
        Hsig_buoy.Hsig(ii) = NaN;
        
        Dir_buoy.dir1(:,ii) = NaN(size(Dir_buoy.dir1(:,ii-1),1),1);
        Dir_buoy.spread1(:,ii) = NaN(size(Dir_buoy.spread1(:,ii-1),1),1);
        Dir_buoy.dir2(:,ii) = NaN(size(Dir_buoy.dir2(:,ii-1),1),1);
        Dir_buoy.spread2(:,ii) = NaN(size(Dir_buoy.spread2(:,ii-1),1),1);
        Dir_buoy.dir1_ss_sum(ii) = NaN;
        Dir_buoy.spread1_ss_sum(ii) = NaN;
        Dir_buoy.dir2_ss_sum(ii) = NaN;
        Dir_buoy.spread2_ss_sum(ii) = NaN;
                
        FC_buoy.a1(:,ii) = NaN(size(FC_buoy.a1(:,ii-1),1),1);
        FC_buoy.b1(:,ii) = NaN(size(FC_buoy.b1(:,ii-1),1),1);
        FC_buoy.a2(:,ii) = NaN(size(FC_buoy.a2(:,ii-1),1),1);
        FC_buoy.b2(:,ii) = NaN(size(FC_buoy.b2(:,ii-1),1),1);
        
        RS_buoy.Sxx(:,ii) = NaN(size(RS_buoy.Sxx(:,ii-1),1),1);
        RS_buoy.Syy(:,ii) = NaN(size(RS_buoy.Syy(:,ii-1),1),1);
        RS_buoy.Sxy(:,ii) = NaN(size(RS_buoy.Sxy(:,ii-1),1),1);
        RS_buoy.Sxx_ss(ii) = NaN;
        RS_buoy.Syy_ss(ii) = NaN;
        RS_buoy.Sxy_ss(ii) = NaN;
        
        Tp_buoy(ii) = NaN;
        
        Eflux_buoy.Epos(:,ii) = NaN(size(Eflux_buoy.Epos(:,ii-1),1),1);
        Eflux_buoy.Eneg(:,ii) = NaN(size(Eflux_buoy.Eneg(:,ii-1),1),1);
        Eflux_buoy.Fpos(:,ii) = NaN(size(Eflux_buoy.Fpos(:,ii-1),1),1);
        Eflux_buoy.Fneg(:,ii) = NaN(size(Eflux_buoy.Fneg(:,ii-1),1),1);
        Eflux_buoy.depth(ii) = NaN;
    end    
end
i_ig_buoy = PUV_process_buoy(1).ids.i_ig;
i_swell_buoy = PUV_process_buoy(1).ids.i_swell;
fm_buoy = PUV_process_buoy(1).Spec.fm;
df_buoy = fm_buoy(2)-fm_buoy(1);
toc


toc
disp('Making QA/QC plots')
%% Ztest check
figure(1);clf
subplot(211)
plot(time(1,:), Ztest_sn.ztest_ss_sum)
hold on
title('Ztest SS')
set(gca, 'FontSize', 20)
grid on
subplot(212)
plot(time(1,:), Ztest_sn.ztest_ig_sum)
title('Ztest IG')
set(gca, 'FontSize', 20)
grid on
hline(1.2)
hline(0.8)
set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'Ztest_sn.png')
%% Mean direction (f)
id=1072
disp(['Selected comparison time is ' datestr(time(1,id)) ' , Hs = ' num2str(Hsig_sn.Hs(id),'%2.2f') 'm'])

figure(2);clf
subplot(121)
polarplot(deg2rad(Dir_sn.dir1(i_ig_sn,id)), fm_sn(i_ig_sn), '.')
hold on
polarplot(deg2rad(Dir_sn.dir2(i_ig_sn,id)), fm_sn(i_ig_sn), '.')
ax = gca;
ax.ThetaDir = 'counterclockwise';
title('IG')
subplot(122)
polarplot(deg2rad(Dir_sn.dir1(i_swell_sn,id)), fm_sn(i_swell_sn), '.')
hold on
ax=gca;
polarplot(deg2rad(Dir_sn.dir2(i_swell_sn,id)), fm_sn(i_swell_sn), '.')
ax.ThetaDir = 'counterclockwise';
ax.RLim=[0 0.25];
ax.RTick = [0 0.05 0.1 0.15 0.2 0.25];
title('SS')
%legend('dir1', 'dir2')
sgtitle([string(time(1,id)) '- SN Coordinates'])

set(gcf, 'Position', [100,100,1000,600])
%saveas(gcf,['Dir_polar_id' num2str(id) '_SN_CCW.png'])
%% Reflection Coeff
figure(3);clf
subplot(311)
plot(time(1,:), Hsig_sn.Hsss)
hold on
plot(time(1,:), Hsig_sn.Hsig)
subplot(312)
plot(time(1,:), R2_ss_sn)
subplot(313)
plot(time(1),1)
hold on
plot(time(1,:), R2_ig_sn)
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
plot(time(1,:), Dir_sn.dir1_ss_sum)
hold on
plot(time(1,:), Dir_sn.dir2_ss_sum)
title('Mean direction SS')
legend('Mean direction (b1/a1)', 'Peak direction (b2/a2)')
set(gca, 'FontSize', 20)
grid on
subplot(212)
plot(time(1,:), Dir_sn.spread1_ss_sum)
hold on
plot(time(1,:), Dir_sn.spread2_ss_sum)
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
h = pcolor(time(1,:)',fm_sn, ph(1:901,:));
c = colorbar('location','WestOutside');
set(h, 'EdgeColor', 'none')
c.Label.String = 'Phase'
set(gca, 'FontSize', 20)
caxis([-180 180])
subplot(212)
h = pcolor(time(1,:)',fm_sn, coh(1:901,:));
c = colorbar('location','WestOutside');
set(h, 'EdgeColor', 'none')
c.Label.String = 'Phase';
set(gca, 'FontSize', 20)
c.Label.String = 'Coherence';
sgtitle('PU phase and coherence')

set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'PU_phase_coherence.png')
%% Get MOP and Buoy data
%mop = filename(15:17);
mop='582'
MOP = read_MOPline(['D0' char(mop)],time(1,1)-hours(1),time(end,end));

dsurl = 'https://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/100p1/100p1_historic.nc'
%ds=ncdataset(dsurl)
%GET TIMES OF DATA
%timevar = ncread(dsurl,'waveTime');
timeconvert = ncread(dsurl,'waveTime'); % Convert UNIX timestamps to Matlab serial units
Buoy.time = datetime(timeconvert, 'ConvertFrom', 'datenum');
Buoy.Hs = ncread(dsurl,'waveHs');
Buoy.spec1D = ncread(dsurl,'waveEnergyDensity');
Buoy.Fq = ncread(dsurl,'waveFrequency');
Buoy.Bw = ncread(dsurl,'waveBandwidth')';
Buoy.Dp = ncread(dsurl,'waveDp');
Buoy.Tp = ncread(dsurl,'waveTp');
Buoy.a1 = ncread(dsurl,'waveA1Value');
Buoy.a2 = ncread(dsurl,'waveA2Value');
Buoy.b1 = ncread(dsurl,'waveB1Value');
Buoy.b2 = ncread(dsurl,'waveB2Value');


%% Compare bulk parameters to MOP
figure(6);clf
subplot(411)
plot(MOP.time,MOP.Hs)
hold on
plot(Buoy.time, Buoy.Hs)
plot(time(1,:), Hsig_sn.Hsss)
plot(time(1,:), Hsig_buoy.Hsss)
legend('MOP','Buoy', 'PUV SN', 'PUV Buoy')
xlim([time(1) time(end)])
title('Hs')
ylabel('m')
set(gca, 'FontSize', 20)

subplot(412)
plot(MOP.time,1./MOP.fp)
hold on
plot(Buoy.time, Buoy.Tp)
plot(time(1,:), Tp_sn)
plot(time(1,:), Tp_buoy)
legend('MOP','Buoy', 'PUV SN', 'PUV Buoy')
xlim([time(1) time(end)])
ylim([0 25])
title('Tp')
ylabel('s')
set(gca, 'FontSize', 20)

subplot(413)
plot(MOP.time,MOP.Sxx)
hold on
plot(time(1,:), RS_sn.Sxx_ss)
plot(time(1,:), RS_buoy.Sxx_ss)
legend('MOP', 'PUV SN', 'PUV Buoy')
xlim([time(1) time(end)])
title('Sxx')
ylabel('m^2')
set(gca, 'FontSize', 20)

subplot(414)
plot(MOP.time,MOP.Sxy)
hold on
plot(time(1,:), RS_sn.Sxy_ss)
plot(time(1,:), RS_buoy.Sxy_ss)
legend('MOP', 'PUV SN', 'PUV Buoy')
xlim([time(1) time(end)])
title('Sxy')
ylabel('m^2')
set(gca, 'FontSize', 20)

set(gcf, 'Position', [100,100,1500,1300])
saveas(gcf,'MOPvPUV_bulk.png')
%% Reflection vs Water Depth
figure(7);clf
subplot(211)
scatter(Eflux_sn.depth, R2_ss_sn, 50, Dir_sn.spread1_ss_sum,'filled')
ylabel('R^2_{SS}')
set(gca, 'FontSize', 20)
grid on
cb = colorbar;
cb.Label.String = 'Spread1 (deg)';
xlabel('<-- Low tide --- Water Depth (m) --- High tide -->')
title('SS')
subplot(212)
scatter(Eflux_sn.depth, R2_ig_sn, 'filled')
ylabel('R^2_{IG}')
set(gca, 'FontSize', 20)
grid on
title('IG')
xlabel('<-- Low tide --- Water Depth (m) --- High tide -->')

set(gcf, 'Position', [100,100,1000,800])
saveas(gcf,'R2_vtide_cspread1.png')
%% Eflux vs Tide
figure(8);clf
subplot(211)
plot(time(1,:),sum(Eflux_sn.Epos(i_ig_sn,:),1)./sum(Eflux_sn.Epos(i_swell_sn,:),1))
hold on
ylim([-0.015 0.03])
ylabel('E^+_{IG}/ E^+_{SS}');set(gca, 'FontSize', 20)
yyaxis right
plot(time(1,:),Eflux_sn.depth)
ylim([10 13])
xlim([time(1,10) time(1,400)])
ylabel('Water depth (m)') 
subplot(212)
scatter(sum(Eflux_sn.Epos(i_swell_sn,:),1),sum(Eflux_sn.Epos(i_ig_sn,:),1), 50, Eflux_sn.depth, 'filled')
cb = colorbar;
cb.Label.String = 'Water depth (m)';
xlabel('E^+_{SS}')
ylabel('E^+{IG}')
set(gca, 'FontSize', 20)

set(gcf, 'Position', [100,100,1000,800])
saveas(gcf,'EigEss_v_depth.png')
%% Directional Spreading Function
clear d_sn ds_sn
%compute directional spreding function using MEM estimator
for ii = 1:segtotal
    ii
    d_sn(:,:,ii)= mem_est(FC_sn.a1(:,ii), FC_sn.a2(:,ii), FC_sn.b1(:,ii), FC_sn.b2(:,ii));
    if isempty(PUV_process(ii).Spec)
        PUV_process(ii).Spec.SSE = NaN(length(fm_sn),1);
    end
    for i=1:length(fm_sn) % loop through freq bands
        ds_sn(i,:,ii)=d_sn(i,:,ii)*PUV_process(ii).Spec.SSE(i); % mutiply by the freq band total energy
    end

end

%% Directional Spreading Function ENU FOR
% clear d_ENU ds_ENU
% %compute directional spreding function using MEM estimator
% for ii = 1:segtotal
%     
%     d_ENU(:,:,ii)= getmem(FC_ENU.a1(:,ii)', FC_ENU.a2(:,ii)', FC_ENU.b1(:,ii)', FC_ENU.b2(:,ii)');
%     if isempty(PUV_process_ENU(ii).Spec)
%         PUV_process_ENU(ii).Spec.SSE = NaN(length(fm),1);
%     end
%     for i=1:length(fm) % loop through freq bands
%         ds_ENU(i,:,ii)=d_ENU(i,:,ii)*PUV_process_ENU(ii).Spec.SSE(i); % mutiply by the freq band total energy
%     end
% 
% end
%% Directional Spreading Function MOP FOR
clear d_buoy ds_buoy
%compute directional spreding function using MEM estimator
for ii = 1:segtotal
    ii
    d_buoy(:,:,ii)= mem_est(FC_buoy.a1(:,ii), FC_buoy.a2(:,ii), FC_buoy.b1(:,ii), FC_buoy.b2(:,ii));
    if isempty(PUV_process_buoy(ii).Spec)
        PUV_process_buoy(ii).Spec.SSE = NaN(length(fm_buoy),1);
    end
    for i=1:length(fm_buoy) % loop through freq bands
        ds_buoy(i,:,ii)=d_buoy(i,:,ii)*PUV_process_buoy(ii).Spec.SSE(i); % mutiply by the freq band total energy
    end

end

%% Checks on directional spectra
aa_sn = ds_sn(:,1:360,id);
ab_sn = sum(aa_sn,2);

figure(9);clf
plot(fm_sn,ab_sn, 'LineWidth', 5)
hold on
plot(fm_sn,PUV_process(id).Spec.SSE, 'LineWidth', 3)
title('E(f) = \int_0^{2\pi}E(f,\theta)d\theta')
legend('\int_0^{2\pi}E(f,\theta)d\theta', 'E(f)')
set(gca, 'FontSize',  20)

set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'Ef_vEftheta.png')

figure(10);clf
plot(fm_sn, sum(d_sn(:,:,id),2))
title('\int_0^{2\pi}G(\theta,f)d\theta = 1')
set(gca, 'FontSize', 20)

set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'intG_check.png')



aa_buoy = ds_buoy(:,1:360,id);
ab_buoy = sum(aa_buoy,2);

figure(10);clf
plot(fm_sn,ab_buoy, 'LineWidth', 5)
hold on
plot(fm_sn,PUV_process_buoy(id).Spec.SSE, 'LineWidth', 3)
title('E(f) = \int_0^{2\pi}E(f,\theta)d\theta')
legend('\int_0^{2\pi}E(f,\theta)d\theta', 'E(f)')
set(gca, 'FontSize',  20)

set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'Ef_vEftheta.png')

figure(10);clf
plot(fm_sn, sum(d_buoy(:,:,id),2))
title('\int_0^{2\pi}G(\theta,f)d\theta = 1')
set(gca, 'FontSize', 20)

set(gcf, 'Position', [100,100,1000,600])
saveas(gcf,'intG_check.png')


%%
id=1072
id_mop = find(min(abs(time(1, id) - MOP.time)) == abs(time(1, id) - MOP.time))
id_buoy = find(min(abs(time(1, id) - Buoy.time)) == abs(time(1, id) - Buoy.time))
figure(11);clf
theta = 0:359;
subplot(221)
ab = sum(ds_sn(:,1:360,id).*cosd(theta),2);
%plot(fm_sn,ab)
hold on
plot(fm_sn, FC_sn.a1(:,id))
plot(fm_buoy, FC_buoy.a1(:,id))
plot(MOP.frequency, MOP.a1(id_mop,:))
plot(Buoy.Fq, Buoy.a1(id_buoy,:))
%legend('\int_0^{2\pi}G(f,\theta)cos{\theta}d\theta', 'a1 (SN)', 'a1 (Buoycoords)', 'mop', 'buoy')
legend('a1 (SN)', 'a1 (Buoycoords)', 'mop', 'buoy')
title('a1')
set(gca, 'FontSize', 20)

subplot(222)
ab = sum(ds_sn(:,1:360,id).*cosd(2*theta),2);
%plot(fm_sn,ab)
hold on
plot(fm_sn, FC_sn.a2(:,id))
plot(fm_buoy, FC_buoy.a2(:,id))
plot(MOP.frequency, MOP.a2(id_mop,:))
plot(Buoy.Fq, Buoy.a2(id_buoy,:))
%legend('\int_0^{2\pi}G(f,\theta)cos{2\theta}d\theta', 'a2 (SN)', 'a2 (Buoycoords)', 'mop', 'buoy')
legend('a2 (SN)', 'a2 (Buoycoords)', 'mop', 'buoy')
title('a2')
set(gca, 'FontSize', 20)

subplot(223)
ab = sum(ds_sn(:,1:360,id).*sind(theta),2);
%plot(fm_sn,ab)
hold on
plot(fm_sn, FC_sn.b1(:,id))
plot(fm_buoy, FC_buoy.b1(:,id))
plot(MOP.frequency, MOP.b1(id_mop,:))
plot(Buoy.Fq, Buoy.b1(id_buoy,:))
%legend('\int_0^{2\pi}G(f,\theta)sin{\theta}d\theta', 'b1 (SN)', 'b1 (Buoycoords)', 'mop', 'buoy')
legend( 'b1 (SN)', 'b1 (Buoycoords)', 'mop', 'buoy')
title('b1')
set(gca, 'FontSize', 20)

subplot(224)
ab = sum(ds_sn(:,1:360,id).*sind(2*theta),2);
%plot(fm_sn,ab)
hold on
plot(fm_sn, FC_sn.b2(:,id))
plot(fm_buoy, FC_buoy.b2(:,id))
plot(MOP.frequency, MOP.b2(id_mop,:))
plot(Buoy.Fq, Buoy.b2(id_buoy,:))
%legend('\int_0^{2\pi}G(f,\theta)sin{2\theta}d\theta', 'b2 (SN)', 'b2 (Buoycoords)', 'mop', 'buoy')
legend( 'b2 (SN)', 'b2 (Buoycoords)', 'mop', 'buoy')
title('b2')
set(gca, 'FontSize', 20)

% set(gcf, 'Position', [100,100,2000,1000])
%saveas(gcf,'mop_puv_dirspec_FC.png')

%% Plot Dirspectra
ds_sn(:,361,:)=ds_sn(:,1,:);
figure(13);clf
polarPcolor(fm_sn', [0:360], ds_sn(:,:,id), 'Nspokes',13)%,'typeRose','default')

%h1 = pcolor(fm', [0:360], log(ds(:,:,83))')
sgtitle('Dirspec from PUV in Shorenormal Coords')
%set(h1, 'EdgeColor', 'none')

set(gcf, 'Position', [100,100,600,600])
saveas(gcf,['Dirspec_id' num2str(id) '_SN.png'])
%%
ds_buoy(:,361,:)=ds_buoy(:,1,:);
figure(2);clf
polarPcolor(fm_sn', [0:360], (ds_buoy(:,:,id)'),'Nspokes',13)%, 'XAngle', 270-shorenormal+180
%h2=pcolor(fm', [0:360], log(ds_MOP(:,:,83))')%
sgtitle('Dirspec from PUV in Buoy Coords')
%caxis([-15 -5])
%set(h, 'EdgeColor', 'none')
%%
clear ds_mopmop d_mopmop
    d_mopmop=getmem(MOP.a1(id_mop,:), MOP.a2(id_mop,:), MOP.b1(id_mop,:), MOP.b2(id_mop,:));
    
    for i=1:length(MOP.fbw) % loop through freq bands
        ds_mopmop(i,:)=d_mopmop(i,:)*MOP.spec1D(id_mop,i)'; % mutiply by the freq band total energy
    end
ds_mopmop(:,361)=ds_mopmop(:,1);
figure(14);clf
% polarPcolor(double(Fq_mop'), [0:360], log(ds_mopmop)', 'XAngle', 270+180,'Nspokes',13)
polarPcolor(double(MOP.frequency)', [0:360], ds_mopmop', 'Nspokes',13)%,'typeRose','meteo')
sgtitle('Dirspec from MOP')
%%
clear ds_buoybuoy d_buoybuoy
    d_buoybuoy=getmem(Buoy.a1(id_buoy,:), Buoy.a2(id_buoy,:), Buoy.b1(id_buoy,:), Buoy.b2(id_buoy,:));
    
    for i=1:length(Buoy.Bw) % loop through freq bands
        ds_buoybuoy(i,:)=d_buoybuoy(i,:)*Buoy.spec1D(id_buoy,i)'; % mutiply by the freq band total energy
    end
ds_buoybuoy(:,361)=ds_buoybuoy(:,1);
figure(15);clf
% polarPcolor(double(Fq_mop'), [0:360], log(ds_mopmop)', 'XAngle', 270+180,'Nspokes',13)
polarPcolor(double(Buoy.Fq)', [0:360], ds_buoybuoy', 'Nspokes',13)%,'typeRose','meteo')
sgtitle('Dirspec from Buoy')
% caxis([-15 -5])
