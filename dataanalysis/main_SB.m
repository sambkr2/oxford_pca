%% Load
myData = load( '../dataformatting/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored_Processed.mat' );
DataNamePrefix = 'CTPCR12p5T1C33DVA';
PODData = myData.PODData;
MaskedData = myData.MaskedData;
% Tumble_280 = load( 'TumbleCR11T2C33DVA_CA_280.mat' );

%% Set Parameters
AnalysisResult.CrankAngle = [ -90 ];                                   % Change this line to allow more crank angles (avaiable from -295 to -60 CAD aTDCf)
AnalysisResult.CycleNo = 1:300;                                   % Do not change this line

[ ~, AnalysisResult.CrankAngleIndex ] = ismember( AnalysisResult.CrankAngle, PODData.CrankAngle );

t_cad = num2str(AnalysisResult.CrankAngle);
ccm_cad = strrep(t_cad,'-','m');

%% Perform POD
for ca_No = 1 : length( AnalysisResult.CrankAngleIndex )
    
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

%     AnalysisResult.PODResult{ ca_No } = temp_PODResult;
    PODResult{ ca_No } = temp_PODResult;
    
    temp_save_name_file = [ DataNamePrefix, '_CA_', num2str( abs( CurrentCrankAngle ) ) ];

    fprintf( 'Saving data to current folder... \n' )
    save( matlab.lang.makeValidName( temp_save_name_file ), 'CurrentCrankAngle', 'PODResult' );
    fprintf( 'Data saved... \n' )
end

%% POD approx paramters
nModes = [ 0 1 2 5 8 20 299 ];
CycleNo = 96;

% Create matrix of ensemble mean for subtraction
meanSB = repmat(PODResult{1,1}.EnsembleMean, 1, 300);

%% 
figureprop.axes_lim = [ -25 25 -30 10 ];
figureprop.xlabel = '{\it x} (mm)';
figureprop.ylabel = '{\it z} (mm)';

figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 15 ];

%% Plot POD approx
% structArray1 = struct('field1',[], 'field2',[], 'field3',[]);

for mm = 1 : length( nModes )
    [ PODApprox ] = Calc_PODApprox_SB( PODResult{1,1}, nModes(mm), CycleNo, PODResult{1,1}.nRowsInOriginalGrid, PODResult{1,1}.nColsInOriginalGrid, PODResult{1,1}.IndexInOriginalGrid, meanSB );
%     PODApprox.X = MaskedData.X;
%     PODApprox.Y = MaskedData.Y;
    PODApprox.X = myData.InterpolatedData.X;
    PODApprox.Y = myData.InterpolatedData.Y;
    PODVel.(['POD',num2str(nModes(mm))]).u = PODApprox.U;
    PODVel.(['POD',num2str(nModes(mm))]).v = PODApprox.V;
    figure_output = ColourQuiver_SB( PODApprox.X, PODApprox.Y, PODApprox.U, PODApprox.V, figureprop );
    ylim([-20 2])
    title( [ 'POD Approx., Order = ', num2str( nModes(mm) )] );
%     export_fig( [ 'TP Cycle ', num2str( cycle_No ), ' POD Approx at -270 CAD aTDCf' ], '-pdf', '-nocrop', '-append' )
%     close all
end

%% Load CFD
load('T1rke/ccm_T1_mot.mat')

%% Plot CFD
temp_x = myData.CFDData.Data.x_PIVGrid;
temp_y = myData.CFDData.Data.z_PIVGrid;
temp_CCM_u = ccmdata.(ccm_cad).u;
temp_CCM_v = ccmdata.(ccm_cad).w;
temp_CCM_SpeedMap = abs( complex( temp_CCM_u, temp_CCM_v ) );

% Add PIV mask to CFD 
[ ~, PIV_CAindex ] = ismember( AnalysisResult.CrankAngle, myData.InterpolatedData.CrankAngle );

temp_PIVem_u = nanmean( myData.InterpolatedData.U( :,:,PIV_CAindex,: ), 4 );
temp_PIVem_v = nanmean( myData.InterpolatedData.V( :,:,PIV_CAindex,: ), 4 );
temp_PIVem_SpeedMap = abs( complex( temp_PIVem_u, temp_PIVem_v ) );

temp_PIV_mask = ~isnan( temp_PIVem_SpeedMap );
temp_PIV_mask = double( temp_PIV_mask );
temp_PIV_mask( temp_PIV_mask==0 ) = NaN;
temp_CCM_u = temp_CCM_u .* temp_PIV_mask;
temp_CCM_v = temp_CCM_v .* temp_PIV_mask;
temp_CCM_SpeedMap = temp_CCM_SpeedMap .* temp_PIV_mask;

ColourQuiver( temp_x, temp_y, temp_CCM_u, temp_CCM_v, figureprop )
