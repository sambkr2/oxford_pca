% Codes to perform PIV dimenstion reduction techniques
% Author(s): Li (Sam) Shen
% sam-li.shen@eng.ox.ac.uk
% Last updated date: 2021.10.31 (checked on 2022.04.26)

%% Clean up
clear all %#ok<CLALL>
close all
clc

%% Load data
myData = matfile( 'x20180810_Tumble_CR11_T2_C33_DVA_Motored_Processed_all.mat' );
DataNamePrefix = 'Apollo19_TumbleCR11T2C33DVA_PODonly';

MaskedData = myData.MaskedData;
PODData = myData.PODData;

%% Parameters setting
AnalysisResult.CrankAngle = [ -90 ];                                   % Change this line to allow more crank angles (avaiable from -295 to -60 CAD aTDCf)
AnalysisResult.CycleNo = 1:300;                                   % Do not change this line

[ ~, AnalysisResult.CrankAngleIndex ] = ismember( AnalysisResult.CrankAngle, PODData.CrankAngle );

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
calc.nModes = 5;
calc.CrankAngle = -280;
calc.CycleNo = 1:300;
[ ~, calc.CrankAngleIndex ] = ismember( calc.CrankAngle, AnalysisResult.CrankAngle );

test_PODResult = AnalysisResult.PODResult{ calc.CrankAngleIndex };
[ test_PODApprox ] = Calc_PODApprox( test_PODResult, calc.nModes, calc.CycleNo );

