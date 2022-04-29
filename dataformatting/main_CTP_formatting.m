%% Clean up
clear all %#ok<CLALL>
close all
clc

%% Load data
myData = matfile( 'x20190517_CrossTumble_CR12p5_T1_C33_DVA_Motored.mat' );

PIVData = myData.PIVData;
Parameters = myData.Parameters;
CFDData = myData.CFDData;
EngineData = myData.EngineData;

RawData.X = PIVData.Data.y;
RawData.Y = PIVData.Data.z;
RawData.U = PIVData.Data.v;
RawData.V = PIVData.Data.w;
RawData.CrankAngle = PIVData.Data.CrankAngle;
RawData.nCycle = size( PIVData.Data.CycleNo, 2 );  

%% Extra tumble plane intake valve mask (manual)
close all
CrankAngle_plot = -90;
CycleNo_plot = 96;
[ ~, CrankAngleNo_plot ] = ismember( CrankAngle_plot, RawData.CrankAngle );

figureprop.axes_lim = [ -25 25 -20 2 ];
figureprop.xlabel = '{\it x} (mm)';
figureprop.ylabel = '{\it z} (mm)';
figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 15 ];
figureprop.velocity_normalisation = 1;

ColourQuiver( RawData.X, RawData.Y, RawData.U( :, :, CrankAngleNo_plot, CycleNo_plot ), RawData.V( :, :, CrankAngleNo_plot, CycleNo_plot ), figureprop );
title('Before Masking')

IntakeValve_Mask.TP_CA280.Point1 = [ -24 0 ]; % Tumble_CR11_T2_C33_DVA, -280 CAD aTDCf
IntakeValve_Mask.TP_CA280.Point2 = [ -24 0.01 ]; % Tumble_CR11_T2_C33_DVA, -280 CAD aTDCf

IntakeValve_Mask.TP_CA280.Slope = ( IntakeValve_Mask.TP_CA280.Point2(2) - IntakeValve_Mask.TP_CA280.Point1(2) ) / ( IntakeValve_Mask.TP_CA280.Point2(1) - IntakeValve_Mask.TP_CA280.Point1(1) );
IntakeValve_Mask.TP_CA280.MaskedOut = ( RawData.Y > ( IntakeValve_Mask.TP_CA280.Point1(2) + IntakeValve_Mask.TP_CA280.Slope * ( RawData.X - IntakeValve_Mask.TP_CA280.Point1(1) ) ) );

IntakeValve_Mask.TP_CA280.Mask = ones( size( RawData.X ) );
IntakeValve_Mask.TP_CA280.Mask( IntakeValve_Mask.TP_CA280.MaskedOut ) = NaN;

% Only apply to -280 CAD aTDCf
% IntakeValve_MaskedData.U = RawData.U( :, :, CrankAngleNo_plot, : ) .* repmat( IntakeValve_Mask.TP_CA280.Mask, [ 1 1 1 RawData.nCycle ] );
% IntakeValve_MaskedData.V = RawData.V( :, :, CrankAngleNo_plot, : ) .* repmat( IntakeValve_Mask.TP_CA280.Mask, [ 1 1 1 RawData.nCycle ] );
% ColourQuiver( RawData.X, RawData.Y, IntakeValve_MaskedData.U( :, :, 1, CycleNo_plot ), IntakeValve_MaskedData.V( :, :, 1, CycleNo_plot ), figureprop );

% Apply to all cycles
IntakeValve_MaskedData.U = RawData.U .* repmat( IntakeValve_Mask.TP_CA280.Mask, [ 1 1 length(RawData.CrankAngle) RawData.nCycle ] );
IntakeValve_MaskedData.V = RawData.V .* repmat( IntakeValve_Mask.TP_CA280.Mask, [ 1 1 length(RawData.CrankAngle) RawData.nCycle ] );
ColourQuiver( RawData.X, RawData.Y, IntakeValve_MaskedData.U( :, :, CrankAngleNo_plot, CycleNo_plot ), IntakeValve_MaskedData.V( :, :, CrankAngleNo_plot, CycleNo_plot ), figureprop );

title('After Masking')

%% Format data
[ InterpolatedData, MaskedData, PODData ] = PIVDataFormatting( RawData.X, RawData.Y, IntakeValve_MaskedData.U, IntakeValve_MaskedData.V, RawData.CrankAngle );

%% Check data
CheckData.CycleNo = 96;
CheckData.CrankAngle = -280;

[ ~, CheckData.CrankAngleNo_Raw ] = ismember( CheckData.CrankAngle, RawData.CrankAngle );
[ ~, CheckData.CrankAngleNo_Using ] = ismember( CheckData.CrankAngle, MaskedData.CrankAngle );

close all
figure
hold on
box on
quiver( RawData.X, RawData.Y, RawData.U(:,:,CheckData.CrankAngleNo_Raw,CheckData.CycleNo), RawData.V(:,:,CheckData.CrankAngleNo_Raw,CheckData.CycleNo),'k' )
title( 'Raw Data' )

figure
hold on
box on
quiver( InterpolatedData.X,InterpolatedData.Y, InterpolatedData.U(:,:, CheckData.CrankAngleNo_Using,CheckData.CycleNo),InterpolatedData.V(:,:, CheckData.CrankAngleNo_Using,CheckData.CycleNo),'b' )
quiver( MaskedData.X,MaskedData.Y, MaskedData.U(:,:, CheckData.CrankAngleNo_Using,CheckData.CycleNo),MaskedData.V(:,:, CheckData.CrankAngleNo_Using,CheckData.CycleNo),'r' )
title( 'interp data (blue) and cycle-masked data (red)' )


figure
box on
quiver( PODData.X{ CheckData.CrankAngleNo_Using },PODData.Y{ CheckData.CrankAngleNo_Using }, PODData.U{ CheckData.CrankAngleNo_Using }(:,CheckData.CycleNo),PODData.V{ CheckData.CrankAngleNo_Using }(:,CheckData.CycleNo),'m' )
title( 'data used for POD (magenta)' )

%% Save data
temp_answer = questdlg( 'Do you want to save the data (into the current MATLAB working folder)',...
    'Data Saving Check', 'Save', 'Do not save', 'Do not save' );
switch temp_answer
    case 'Save'
        savepath = uigetdir;
        save( fullfile( savepath, [ matlab.lang.makeValidName( Parameters.Info.save_name ), '_Processed' ] ), 'Parameters', 'InterpolatedData', 'MaskedData', 'PODData', 'CFDData', 'EngineData' );
        fprintf( 'Data saved under:\n %s \n Data filename is: %s \n', savepath, [ matlab.lang.makeValidName( Parameters.Info.save_name ), '_Processed' ] )
    case 'Do not save'
        fprintf( 'Data saving is skipped by the user.' )
end
