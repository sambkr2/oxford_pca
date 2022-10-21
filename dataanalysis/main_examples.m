% Example codes for PIV and CFD vector field comparison
% Author: Li (Sam) Shen
% Email: sam-li.shen@eng.ox.ac.uk
% Last update date: Dec.21st, 2021

clear
close all
clc

%% load data
% TumbleData = load( 'x20180823_Tumble_CR12p5_T7_C33_DVA_Motored.mat' );
TumbleData = load( "../dataformatting/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored.mat");

%% Plot settings
Plot_FontSize = 14;
Plot_FontName = 'Times New Roman';
Plot_xlabel = 'x';
Plot_ylabel = 'z';
Plot_clabel = 'Flow Speed (m/s)';
Plot_Prop.XLim = [ -30 30 ];
Plot_Prop.YLim = [ -30 10 ];
Plot_Prop.Clim_max = 8;
PIV_SpacingReduce = 1; % reduce PIV vector spacing by a factor for a better visualisation
CFD_SpacingReduce = 1; % reduce CFD vector spacing by a factor for a better visualisation

NorScale = 1; % normalise factor of flow speed, no need for internal use but needed for publication (JLR requires to normalise the flow speed in publications)
BWvector_threshold = 0.7; % thresholds

xPcolorShift = TumbleData.PIVData.Info.y_mmPerPixel * TumbleData.PIVData.Info.y_NoPixelsBetweenVectors / 2;
yPcolorShift = TumbleData.PIVData.Info.z_mmPerPixel * TumbleData.PIVData.Info.z_NoPixelsBetweenVectors / 2;

plot_crankangle = -285;
t_cad = num2str(plot_crankangle);
ccm_cad = strrep(t_cad,'-','m');

%% 
figureprop.axes_lim = [ -25 25 -20 2 ];
figureprop.xlabel = '{\it x} (mm)';
figureprop.ylabel = '{\it z} (mm)';
figureprop.velocity_normalisation = 1;
figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 50 ];

%% Plot example PIV ensemble mean flow field
[ ~, PIV_CAindex ] = ismember( plot_crankangle, TumbleData.PIVData.Data.CrankAngle );

temp_x = TumbleData.PIVData.Data.y;
temp_y = TumbleData.PIVData.Data.z;
temp_PIVem_u = mean( TumbleData.PIVData.Data.v( :,:,PIV_CAindex,: ), 4,'omitnan' );
temp_PIVem_v = mean( TumbleData.PIVData.Data.w( :,:,PIV_CAindex,: ), 4,'omitnan' );
temp_PIVem_SpeedMap = abs( complex( temp_PIVem_u, temp_PIVem_v ) );

figure_output = ColourQuiver( temp_x, temp_y, temp_PIVem_u, temp_PIVem_v, figureprop );


%% Load CFD
load('rrT1rng/ccm_T1_mot_CTP.mat')

%% Plot CFD 
temp_x = TumbleData.CFDData.Data.y_PIVGrid;
temp_y = TumbleData.CFDData.Data.z_PIVGrid;
temp_CCM_u = ccmdata.(ccm_cad).v;
temp_CCM_v = ccmdata.(ccm_cad).w;
temp_CCM_SpeedMap = abs( complex( temp_CCM_u, temp_CCM_v ) );

% Add PIV mask to CFD 
[ ~, PIV_CAindex ] = ismember( plot_crankangle, TumbleData.PIVData.Data.CrankAngle );

temp_PIVem_u = mean( TumbleData.PIVData.Data.v( :,:,PIV_CAindex,: ), 4,'omitnan' );
temp_PIVem_v = mean( TumbleData.PIVData.Data.w( :,:,PIV_CAindex,: ), 4,'omitnan' );
temp_PIVem_SpeedMap = abs( complex( temp_PIVem_u, temp_PIVem_v ) );

temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
temp_PIV_mask = double( temp_PIV_mask );
temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
temp_CCM_u = temp_CCM_u .* temp_PIV_mask;
temp_CCM_v = temp_CCM_v .* temp_PIV_mask;
temp_CCM_SpeedMap = temp_CCM_SpeedMap .* temp_PIV_mask;

ColourQuiver(temp_x, temp_y, temp_CCM_u ,temp_CCM_v, figureprop)
ylim([-20 2])
% title(['CFD ',num2str(CurrentCrankAngle),' CAD'])
% set(gcf, 'position', [440 377 560 290])
% name = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/CTP_CFD_',ccm_cad,'.png'];
% exportgraphics(gcf,name,'resolution',600)


%% Calculate vector comparison metrics
RI.PIVem_StarCD = Find_Relevance_Index( complex( temp_PIVem_u, temp_PIVem_v ), complex( temp_CCM_u, temp_CCM_v ), 'normal_num' ); % relevance index
WRI.PIVem_StarCD = calculate_WRI( temp_PIVem_u, temp_PIVem_v, temp_CCM_u, temp_CCM_v, ones( size(temp_PIVem_u) ) ); % weighted relevance index
WMI.PIVem_StarCD = calculate_WMI( temp_PIVem_u, temp_PIVem_v, temp_CCM_u, temp_CCM_v, ones( size(temp_PIVem_u) ) ); % weighted magnitude index

%% Calculate vector comparison metrics
RI.DMD_StarCD = Find_Relevance_Index( complex( dmdplotu, dmdplotv ), complex( temp_CCM_u, temp_CCM_v ), 'normal_num' ); % relevance index
WRI.DMD_StarCD = calculate_WRI( dmdplotu, dmdplotv, temp_CCM_u, temp_CCM_v, ones( size(dmdplotu) ) ); % weighted relevance index
WMI.DMD_StarCD = calculate_WMI( dmdplotu, dmdplotv, temp_CCM_u, temp_CCM_v, ones( size(dmdplotu) ) ); % weighted magnitude index

%% Calculate vector comparison metrics
RI.EM_StarCD = Find_Relevance_Index( complex( dmdplotu, dmdplotv ), complex( temp_PIVem_u, temp_PIVem_v ), 'normal_num' ); % relevance index
WRI.EM_StarCD = calculate_WRI( dmdplotu, dmdplotv, temp_PIVem_u, temp_PIVem_v, ones( size(dmdplotu) ) ); % weighted relevance index
WMI.EM_StarCD = calculate_WMI( dmdplotu, dmdplotv, temp_PIVem_u, temp_PIVem_v, ones( size(dmdplotu) ) ); % weighted magnitude index

%% Plot vector comparison metrics
temp_x = TumbleData.PIVData.Data.y;
temp_y = TumbleData.PIVData.Data.z;

% WRI
% temp_data = WRI.PIVem_StarCD;
% temp_data = WRI.PIVem_StarCD;
temp_data = WRI.EM_StarCD;
Plot_ColourMapandQuiver_CTP( 5, temp_x, temp_y, [], [], temp_data,...
    xPcolorShift, yPcolorShift, 1,...
    ['{\it{', Plot_xlabel, '}} (mm)'], ['{\it{', Plot_ylabel, '}} (mm)'], 'Weighted Relevance Index', Plot_FontSize, Plot_FontName,...
    Plot_Prop.XLim, Plot_Prop.YLim, 0.1, [], 0, [] );
% ColourQuiver(temp_x, temp_y, [] ,[], figureprop)

% WMI
% temp_data = WMI.PIVem_StarCD;
% temp_data = WMI.PIVem_StarCD;
temp_data = WMI.EM_StarCD;
Plot_ColourMapandQuiver( 6, temp_x, temp_y, [], [], temp_data,...
    xPcolorShift, yPcolorShift, 1,...
    ['{\it{', Plot_xlabel, '}} (mm)'], ['{\it{', Plot_ylabel, '}} (mm)'], 'Weighted Magnitude Index', Plot_FontSize, Plot_FontName,...
    Plot_Prop.XLim, Plot_Prop.YLim, 2, [], 0, [] );



