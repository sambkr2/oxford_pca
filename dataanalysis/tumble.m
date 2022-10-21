% Example codes for PIV and CFD vector field comparison
% Author: Li (Sam) Shen
% Email: sam-li.shen@eng.ox.ac.uk
% Last update date: Dec.21st, 2021

% clear
% close all
% clc

%% load data
% TumbleData = load( 'x20180823_Tumble_CR12p5_T7_C33_DVA_Motored.mat' );
% load('x20180706_Tumble_CR12p5_T1_C33_DVA_Motored.mat');
load('../../JLR/baseline/matlab/x20180706_Tumble_CR12p5_T1_C33_DVA_Motored.mat')

%% Plot settings
Plot_FontSize = 14;
Plot_FontName = 'serif';
Plot_xlabel = 'x';
Plot_ylabel = 'z';
Plot_clabel = 'Flow Speed (m/s)';
Plot_Prop.XLim = [ -30 30 ];
Plot_Prop.YLim = [ -30 10 ];
Plot_Prop.Clim_max = 15;
PIV_SpacingReduce = 2; % reduce PIV vector spacing by a factor for a better visualisation
CFD_SpacingReduce = 2; % reduce CFD vector spacing by a factor for a better visualisation

NorScale = 1; % normalise factor of flow speed, no need for internal use but needed for publication (JLR requires to normalise the flow speed in publications)
BWvector_threshold = 0.7; % thresholds

xPcolorShift = PIVData.Info.x_mmPerPixel * PIVData.Info.x_NoPixelsBetweenVectors / 2;
yPcolorShift = PIVData.Info.z_mmPerPixel * PIVData.Info.z_NoPixelsBetweenVectors / 2;

plot_crankangle = -90;
t_cad = num2str(plot_crankangle);
ccm_cad = strrep(t_cad,'-','m');

%% Plot example PIV ensemble mean flow field
[ ~, PIV_CAindex ] = ismember( plot_crankangle, PIVData.Data.CrankAngle );

temp_x = PIVData.Data.x;
temp_y = PIVData.Data.z;
temp_PIVem_u = nanmean( PIVData.Data.u( :,:,PIV_CAindex,: ), 4 );
temp_PIVem_v = nanmean( PIVData.Data.w( :,:,PIV_CAindex,: ), 4 );
temp_PIVem_SpeedMap = abs( complex( temp_PIVem_u, temp_PIVem_v ) );

Plot_ColourMapandQuiver( 1, temp_x, temp_y, temp_PIVem_u, temp_PIVem_v, temp_PIVem_SpeedMap,...
                            xPcolorShift, yPcolorShift, NorScale,...
                            ['{\it{', Plot_xlabel, '}} (mm)'], ['{\it{', Plot_ylabel, '}} (mm)'], Plot_clabel, Plot_FontSize, Plot_FontName,...
                            Plot_Prop.XLim, Plot_Prop.YLim, Plot_Prop.Clim_max, PIV_SpacingReduce, 1, BWvector_threshold );
set(gcf, 'position', [-700 199 566 370])
labelPIV = 'PIV';
text(-2, 15, labelPIV, 'fontname', 'serif', 'fontsize', 18);
% colorbar('northoutside')

%% Plot example CFD flow field from Star-CD (from JLR)
[~,CFD_CAindex] = ismember( plot_crankangle, CFDData.Data.CrankAngle );

temp_x = CFDData.Data.x_PIVGrid;
temp_y = CFDData.Data.z_PIVGrid;
temp_CFD_u = CFDData.Data.u_PIVGrid( :,:,CFD_CAindex );
temp_CFD_v = CFDData.Data.w_PIVGrid( :,:,CFD_CAindex );
temp_CFD_SpeedMap = abs( complex( temp_CFD_u, temp_CFD_v ) );

% Add PIV mask to CFD
temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
temp_PIV_mask = double( temp_PIV_mask );
temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
temp_CFD_u = temp_CFD_u .* temp_PIV_mask;
temp_CFD_v = temp_CFD_v .* temp_PIV_mask;
temp_CFD_SpeedMap = temp_CFD_SpeedMap .* temp_PIV_mask;

Plot_ColourMapandQuiver( 2, temp_x, temp_y, temp_CFD_u, temp_CFD_v, temp_CFD_SpeedMap,...
    xPcolorShift, yPcolorShift, NorScale,...
    ['{\it{', Plot_xlabel, '}} (mm)'], ['{\it{', Plot_ylabel, '}} (mm)'], Plot_clabel, Plot_FontSize, Plot_FontName,...
    Plot_Prop.XLim, Plot_Prop.YLim, Plot_Prop.Clim_max, CFD_SpacingReduce, 1, BWvector_threshold );

set(gcf, 'position', [-700 199 566 283])
labelCD = '(b) STAR-CD';
text(-13, 15, labelCD, 'fontname', 'serif', 'fontsize', 18);

%% STAR-CCM+
temp_x = PIVData.Data.x;
temp_y = PIVData.Data.z;
temp_CCM_u = ccmdata.(ccm_cad).u;
temp_CCM_v = ccmdata.(ccm_cad).w;
temp_CCM_SpeedMap = abs( complex( temp_CCM_u, temp_CCM_v ) );

% Add PIV mask to CFD 
temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
temp_PIV_mask = double( temp_PIV_mask );
temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
temp_CCM_u = temp_CCM_u .* temp_PIV_mask;
temp_CCM_v = temp_CCM_v .* temp_PIV_mask;
temp_CCM_SpeedMap = temp_CCM_SpeedMap .* temp_PIV_mask;

Plot_ColourMapandQuiver( 3, temp_x, temp_y, temp_CCM_u, temp_CCM_v, temp_CCM_SpeedMap,...
    xPcolorShift, yPcolorShift, NorScale,...
    ['{\it{', Plot_xlabel, '}} (mm)'], ['{\it{', Plot_ylabel, '}} (mm)'], Plot_clabel, Plot_FontSize, Plot_FontName,...
    Plot_Prop.XLim, Plot_Prop.YLim, Plot_Prop.Clim_max, CFD_SpacingReduce, 1, BWvector_threshold );

set(gcf, 'position', [-700 199 566 370])
labelCCM = 'STAR-CCM+ (RNG)';
text(-10, 15, labelCCM, 'fontname', 'serif', 'fontsize', 18);

%% Calculate vector comparison metrics
% PIV and STAR-CD
RI.piv_cd = Find_Relevance_Index( complex( temp_PIVem_u, temp_PIVem_v ), complex( temp_CFD_u, temp_CFD_v ), 'normal_num' ); % relevance index
WRI.piv_cd = calculate_WRI( temp_PIVem_u, temp_PIVem_v, temp_CFD_u, temp_CFD_v, ones( size(temp_PIVem_u) ) ); % weighted relevance index
WMI.piv_cd = calculate_WMI( temp_PIVem_u, temp_PIVem_v, temp_CFD_u, temp_CFD_v, ones( size(temp_PIVem_u) ) ); % weighted magnitude index
WRI.piv_cd_av = nanmean(WRI.piv_cd,'all');
WMI.piv_cd_av = nanmean(WMI.piv_cd,'all');
MSI_piv_cd = calc_MSI(complex( temp_PIVem_u, temp_PIVem_v ), complex( temp_CFD_u, temp_CFD_v ), 'normal_num');

%%
% PIV and STAR-CCM+ 
RI.piv_ccm = Find_Relevance_Index( complex( temp_PIVem_u, temp_PIVem_v ), complex( temp_CCM_u, temp_CCM_v ), 'normal_num' ); % relevance index
WRI.piv_ccm = calculate_WRI( temp_PIVem_u, temp_PIVem_v, temp_CCM_u, temp_CCM_v, ones( size(temp_PIVem_u) ) ); % weighted relevance index
WMI.piv_ccm = calculate_WMI( temp_PIVem_u, temp_PIVem_v, temp_CCM_u, temp_CCM_v, ones( size(temp_PIVem_u) ) ); % weighted magnitude index
WRI.piv_ccm_av = nanmean(WRI.piv_ccm,'all');
WMI.piv_ccm_av = nanmean(WMI.piv_ccm,'all');
MSI_piv_ccm = calc_MSI(complex( temp_PIVem_u, temp_PIVem_v ), complex( temp_CCM_u, temp_CCM_v ), 'normal_num');

%%
% STAR-CD and STAR-CCM+ 
RI.cd_ccm = Find_Relevance_Index( complex( temp_CFD_u, temp_CFD_v ), complex( temp_CCM_u, temp_CCM_v ), 'normal_num' ); % relevance index
WRI.cd_ccm = calculate_WRI( temp_CFD_u, temp_CFD_v, temp_CCM_u, temp_CCM_v, ones( size(temp_PIVem_u) ) ); % weighted relevance index
WMI.cd_ccm = calculate_WMI( temp_CFD_u, temp_CFD_v, temp_CCM_u, temp_CCM_v, ones( size(temp_PIVem_u) ) ); % weighted magnitude index
WRI.cd_ccm_av = nanmean(WRI.cd_ccm,'all');
WMI.cd_ccm_av = nanmean(WMI.cd_ccm,'all');

%% Plot vector comparison metrics
% temp_x = PIVData.Data.x;
% temp_y = PIVData.Data.z;
% 
% % WRI
% temp_data = WRI.PIVem_StarCD;
% Plot_ColourMapandQuiver( 4, temp_x, temp_y, [], [], temp_data,...
%     xPcolorShift, yPcolorShift, 1,...
%     ['{\it{', Plot_xlabel, '}} (mm)'], ['{\it{', Plot_ylabel, '}} (mm)'], 'Weighted Relevance Index', Plot_FontSize, Plot_FontName,...
%     Plot_Prop.XLim, Plot_Prop.YLim, 0.1, [], 0, [] );
% 
% % WMI
% temp_data = WMI.PIVem_StarCD;
% Plot_ColourMapandQuiver( 5, temp_x, temp_y, [], [], temp_data,...
%     xPcolorShift, yPcolorShift, 1,...
%     ['{\it{', Plot_xlabel, '}} (mm)'], ['{\it{', Plot_ylabel, '}} (mm)'], 'Weighted Magnitude Index', Plot_FontSize, Plot_FontName,...
%     Plot_Prop.XLim, Plot_Prop.YLim, 2, [], 0, [] );



