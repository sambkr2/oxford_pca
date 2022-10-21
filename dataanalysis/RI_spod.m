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

ccmu = ccmdata.m270.v;
ccmv = ccmdata.m270.w;

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
