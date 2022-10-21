%% Load (Tumble. Scroll down for cross tumble)
load('../dataformatting/x20180706_Tumble_CR12p5_T1_C33_DVA_Motored_Processed_all_masked.mat')
load('T1rng/ccm_T1_mot.mat');

%%
% Run first few sections of tumble_dmd.m, then
dmd1u = dmdplotu;
dmd1v = dmdplotv;

ccmu = ccmdata.m285.u;
ccmv = ccmdata.m285.w;

emu = mean(MaskedData.U(:,:,CrankAngleNo_plot,:),4,'omitnan');
emv = mean(MaskedData.V(:,:,CrankAngleNo_plot,:),4,'omitnan');
PIVem_SpeedMap = abs( complex( emu, emv ) );

PIV_mask = ~isnan( PIVem_SpeedMap );
PIV_mask = double( PIV_mask );
PIV_mask( PIV_mask==0 ) = NaN;
ccmu = ccmu .* PIV_mask;
ccmv = ccmv .* PIV_mask;

ri1 = Find_Relevance_Index( complex( ccmu, ccmv ), complex( emu, emv ), 'normal_num' )
ri2 = Find_Relevance_Index( complex( ccmu, ccmv ), complex( dmd1u, dmd1v ), 'normal_num' )
ri3 = Find_Relevance_Index( complex( dmd1u, dmd1v ), complex( emu, emv ), 'normal_num' )

ri2 = Find_Relevance_Index( complex( ccmu, ccmv ), complex( dmdplotu, dmdplotv ), 'normal_num' )

%% Load cross-tumble
load('rrT1rng/ccm_T1_mot_CTP.mat')
load('../dataformatting/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored_Processed_all.mat')

%%
dmd1u = dmdplotu;
dmd1v = dmdplotv;

ccmu = ccmdata.m285.v;
ccmv = ccmdata.m285.w;

emu = mean(MaskedData.U(:,:,CrankAngleNo_plot,:),4,'omitnan');
emv = mean(MaskedData.V(:,:,CrankAngleNo_plot,:),4,'omitnan');
PIVem_SpeedMap = abs( complex( emu, emv ) );

PIV_mask = ~isnan( PIVem_SpeedMap );
PIV_mask = double( PIV_mask );
PIV_mask( PIV_mask==0 ) = NaN;
ccmu = ccmu .* PIV_mask;
ccmv = ccmv .* PIV_mask;

ri1 = Find_Relevance_Index( complex( ccmu, ccmv ), complex( emu, emv ), 'normal_num' )
ri2 = Find_Relevance_Index( complex( ccmu, ccmv ), complex( dmd1u, dmd1v ), 'normal_num' )
ri3 = Find_Relevance_Index( complex( dmd1u, dmd1v ), complex( emu, emv ), 'normal_num' )

%% Cycles
ri_dmd = zeros(300, 2);

for k = 1:300
    
    temp_piv_u = MaskedData.U( :, :, CrankAngleNo_plot, k );
    temp_piv_v = MaskedData.V( :, :, CrankAngleNo_plot, k );
%     temp_PIVem_SpeedMap = abs( complex( temp_piv_u, temp_piv_v ) );

    temp_RI_piv_dmd = Find_Relevance_Index( complex( temp_piv_u, temp_piv_v ), complex( dmdplotu, dmdplotv ), 'normal_num' );
    ri_dmd(k,1) = k;
    ri_dmd(k,2) = temp_RI_piv_dmd;

end

%% Plot cycles
figure
axes('fontname','times','fontsize',24,'linewidth',1)
hold on
scatter(ri_dmd(:,1),ri_dmd(:,2),'filled')
% yline(0.7)
box on
% l = legend('ccm vs ind cycle','em vs ind cycle');
% set(l,'location','southwest')
% title([str_cad,' CAD'])
xlabel('Cycle Number')
ylabel('Relevance Index')
ylim([0 1])
mean(ri_dmd(:,2))
