%% Clean up
clear all %#ok<CLALL>
close all
clc

%% Load data
myData = matfile( '../dataformatting/x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored_Processed.mat' );
DataNamePrefix = 'CTPCR12p5T1C33DVA_PODonly';

MaskedData = myData.MaskedData;
PODData = myData.PODData;

%% Parameters setting
AnalysisResult.CrankAngle = [ -90 ];                                   % Change this line to allow more crank angles (avaiable from -295 to -60 CAD aTDCf)
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
    temp_PODResult = Perform_POD( temp_velo_data, 'Centered', 'Direct' );

    temp_PODResult.X_PODGrid = PODData.X{ AnalysisResult.CrankAngleIndex( ca_No ) };
    temp_PODResult.Y_PODGrid = PODData.Y{ AnalysisResult.CrankAngleIndex( ca_No ) };
    temp_PODResult.IndexInOriginalGrid = PODData.IndexInOriginal{ AnalysisResult.CrankAngleIndex( ca_No ) };
    temp_PODResult.nRowsInOriginalGrid = size( MaskedData.X, 1 );
    temp_PODResult.nColsInOriginalGrid = size( MaskedData.X, 2 );
    temp_PODResult.X_OriginalGrid = MaskedData.X;
    temp_PODResult.Y_OriginalGrid = MaskedData.Y;

    AnalysisResult.PODResult{ ca_No } = temp_PODResult;
end

%% Calculate POD approximations
% calc.nModes = [ 0 1 2 5 8 20 299 ];
% calc.CrankAngle = -90;
% calc.CycleNo = 1:300;
% [ ~, calc.CrankAngleIndex ] = ismember( calc.CrankAngle, AnalysisResult.CrankAngle );
% 
% test_PODResult = AnalysisResult.PODResult{ calc.CrankAngleIndex };
% [ test_PODApprox ] = Calc_PODApprox( test_PODResult, calc.nModes, calc.CycleNo );

%% POD approx paramters
nModes = [ 0 1 2 5 8 20 299 ];
CycleNo = 96;

% Create matrix of ensemble mean for subtraction
% meanSB = repmat(PODResult{1,1}.EnsembleMean, 1, 300);

%% 
figureprop.axes_lim = [ -25 25 -30 10 ];
figureprop.xlabel = '{\it x} (mm)';
figureprop.ylabel = '{\it z} (mm)';

figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 15 ];

%%
PODResult = AnalysisResult.PODResult;
InterpolatedData = myData.InterpolatedData;

%% Plot POD approx
% structArray1 = struct('field1',[], 'field2',[], 'field3',[]);

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