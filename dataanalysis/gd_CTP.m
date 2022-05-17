%% Load data


%% Parameters setting
AnalysisResult.CrankAngle = [ -270 ];                                   % Change this line to allow more crank angles (avaiable from -295 to -60 CAD aTDCf)
AnalysisResult.CycleNo = 1:300;                                   % Do not change this line

[ ~, AnalysisResult.CrankAngleIndex ] = ismember( AnalysisResult.CrankAngle, PODData.CrankAngle );

t_cad = num2str(AnalysisResult.CrankAngle);
ccm_cad = strrep(t_cad,'-','m');

%% POD Analysis
AnalysisResult.PODResult = cell( length( AnalysisResult.CrankAngleIndex ), 1 );
for ca_No = 1 : length( AnalysisResult.CrankAngleIndex )
    % Perform POD
    CurrentCrankAngle = AnalysisResult.CrankAngle( ca_No );
    fprintf( 'CA = %.0f CAD aTDCf \n', CurrentCrankAngle )

    temp_velo_data = complex( PODData.U{ AnalysisResult.CrankAngleIndex( ca_No ) }, PODData.V{ AnalysisResult.CrankAngleIndex( ca_No ) } );
    temp_PODResult = Perform_POD_SB( temp_velo_data, 'Centered', 'Direct' );

    temp_PODResult.X_PODGrid = PODData.X{ AnalysisResult.CrankAngleIndex( ca_No ) };
    temp_PODResult.Y_PODGrid = PODData.Y{ AnalysisResult.CrankAngleIndex( ca_No ) };
    temp_PODResult.IndexInOriginalGrid = PODData.IndexInOriginal{ AnalysisResult.CrankAngleIndex( ca_No ) };
    temp_PODResult.nRowsInOriginalGrid = size( MaskedData.X, 1 );
    temp_PODResult.nColsInOriginalGrid = size( MaskedData.X, 2 );
    temp_PODResult.X_OriginalGrid = MaskedData.X;
    temp_PODResult.Y_OriginalGrid = MaskedData.Y;

    AnalysisResult.PODResult{ ca_No } = temp_PODResult;
end

%% Gavish Donoho
svs = diag(AnalysisResult.PODResult{1,1}.svdS);
beta = PODResult{1,1}.nColsInOriginalGrid/PODResult{1,1}.nRowsInOriginalGrid;
tau = optimal_SVHT_coef(beta,0) * median(svs); % find cut-off tau
modes = svs(svs>tau);
GDmode = length(modes); % Gavish Donoho threshold mode

%% POD approx paramters
nModes = [ 0 1 2 5 8 20 299 GDmode ];
CycleNo = 95;

%% 
figureprop.axes_lim = [ -25 25 -30 10 ];
figureprop.xlabel = '{\it x} (mm)';
figureprop.ylabel = '{\it z} (mm)';

figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 50 ];

%%
PODResult = AnalysisResult.PODResult;
% InterpolatedData = myData.InterpolatedData;

%% Plot POD approx
for mm = 1 : length( nModes )
    [ PODApprox ] = Calc_PODApprox( PODResult{1,1}, nModes(mm), CycleNo );
%     PODApprox.X = MaskedData.X;
%     PODApprox.Y = MaskedData.Y;
    PODApprox.X = InterpolatedData.X;
    PODApprox.Y = InterpolatedData.Y;
    PODVel.(['POD',num2str(nModes(mm))]).u = PODApprox.U;
    PODVel.(['POD',num2str(nModes(mm))]).v = PODApprox.V;
    figure_output = ColourQuiver_SB( PODApprox.X, PODApprox.Y, PODApprox.U, PODApprox.V, figureprop );
    ylim([-20 2])
    title( [ 'POD Approx., Order = ', num2str( nModes(mm) )] );
%     export_fig( [ 'TP Cycle ', num2str( cycle_No ), ' POD Approx at -270 CAD aTDCf' ], '-pdf', '-nocrop', '-append' )
%     close all
end

%% Load CFD
load('rrT1rng/ccm_T1_mot_CTP.mat')

%% Plot CFD 
temp_x = CFDData.Data.y_PIVGrid;
temp_y = CFDData.Data.z_PIVGrid;
temp_CCM_u = ccmdata.(ccm_cad).v;
temp_CCM_v = ccmdata.(ccm_cad).w;
temp_CCM_SpeedMap = abs( complex( temp_CCM_u, temp_CCM_v ) );

% Add PIV mask to CFD 
[ ~, PIV_CAindex ] = ismember( AnalysisResult.CrankAngle, InterpolatedData.CrankAngle );

temp_PIVem_u = nanmean( InterpolatedData.U( :,:,PIV_CAindex,: ), 4 );
temp_PIVem_v = nanmean( InterpolatedData.V( :,:,PIV_CAindex,: ), 4 );
temp_PIVem_SpeedMap = abs( complex( temp_PIVem_u, temp_PIVem_v ) );

temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
temp_PIV_mask = double( temp_PIV_mask );
temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
temp_CCM_u = temp_CCM_u .* temp_PIV_mask;
temp_CCM_v = temp_CCM_v .* temp_PIV_mask;
temp_CCM_SpeedMap = temp_CCM_SpeedMap .* temp_PIV_mask;

ColourQuiver_SB(temp_x, temp_y, temp_CCM_u ,temp_CCM_v, figureprop)
ylim([-20 2])
title('RNG')