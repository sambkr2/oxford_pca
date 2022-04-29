function [ InterpolatedData, MaskedData, PODData ] = PIVDataFormatting_nomask( X, Y, U, V, CrankAngle )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% 
% if ( X( 1, 1 ) == X( 2, 1 ) ) && ( Y( 1, 1 ) == Y( 1, 2 ) ) % if X has different values in each column and Y has different values in each row
%     if X( 1, 1 ) > X( 1, 2 ) % if X values are decreasing as column number increases
%         X = flip( X, 2 );
%         Y = flip( Y, 2 );
%         U = flip( U, 2 );
%         V = flip( V, 2 );
%     end
%     if Y( 1, 1 ) > Y( 2, 1 ) % if Y values are decreasing as row number increases
%         X = flip( X, 1 );
%         Y = flip( Y, 1 );
%         U = flip( U, 1 );
%         V = flip( V, 1 );
%     end
% elseif ( X( 1, 1 ) == X( 1, 2 ) ) && ( Y( 1, 1 ) == Y( 2, 1 ) ) % if X has different values in each row and Y has different values in each column
%     if X( 1, 1 ) > X( 2, 1 ) % if X values are decreasing as row number increases
%         X = flip( X, 1 );
%         Y = flip( Y, 1 );
%         U = flip( U, 1 );
%         V = flip( V, 1 );
%     end
%     if Y( 1, 1 ) > Y( 1, 2 ) % if Y values are decreasing as column number increases
%         X = flip( X, 2 );
%         Y = flip( Y, 2 );
%         U = flip( U, 2 );
%         V = flip( V, 2 );
%     end
% else
%     error( 'Wrong PIV grid, please check.' )
% end





nCrankAngle_Raw = size( U, 3 );
nCycle = size( U, 4 );

if nCrankAngle_Raw ~= length( CrankAngle )
    error( 'Please check data package format.' )
end

%% Interpolate PIV vectors to avoid holes in the middle
% ValidPoints_beforeInterp = ~( isnan( U ) | isnan( V ) );

% U_interp = zeros( size( U ) );
% V_interp = zeros( size( V ) );
% for ca_No = 1 : nCrankAngle_Raw
%     for cycle_No = 1 : nCycle
%         temp_U = U(:,:,ca_No,cycle_No);
%         temp_V = V(:,:,ca_No,cycle_No);
%         temp_ValidPoints = ValidPoints_beforeInterp(:,:,ca_No,cycle_No);
%         
%         [ U_interp(:,:,ca_No,cycle_No), V_interp(:,:,ca_No,cycle_No) ] =...
%             Interpolate_PIVVectors( X( temp_ValidPoints ), Y( temp_ValidPoints ), temp_U( temp_ValidPoints ), temp_V( temp_ValidPoints ), X, Y );
%         clear temp_*
%     end
% end
% clear ca_No cycle_No
% 
% ValidPoints_afterInterp = ~( isnan( U_interp ) | isnan( V_interp ) );
% nPointsInterp = sum( ValidPoints_afterInterp, [ 1 2 4 ] ) - sum( ValidPoints_beforeInterp, [ 1 2 4 ] ); % number of points interpolated
% nPointsInterpPercentage = nPointsInterp ./ sum( ValidPoints_beforeInterp, [ 1 2 4 ] ) * 100; % in percentage

% Criterion: number of points interpolated for each crank angle must be smaller than 1%
InterpolatedData.X = X;
InterpolatedData.Y = Y;
InterpolatedData.U = U;
InterpolatedData.V = V;
% InterpolatedData.U = U_interp( :,:, ( nPointsInterpPercentage < 1 ) ,: );
% InterpolatedData.V = V_interp( :,:, ( nPointsInterpPercentage < 1 ) ,: );
% InterpolatedData.CrankAngle = CrankAngle( nPointsInterpPercentage < 1 );

InterpolatedData.CrankAngle = CrankAngle;
nCrankAngle_Using = length( InterpolatedData.CrankAngle );

%% Mask data so that each cycle has the same number of vectors (get rid of NaNs on the boundary)
% ValidPoints_CycleMask = ( sum( ValidPoints_afterInterp, 4 ) == nCycle ); % nX * nY * nCA

ValidPoints_CycleMask = 300;
CycleMask = double( ValidPoints_CycleMask ); % nX * nY * nCA
CycleMask( CycleMask == 0 ) = NaN; 
% CycleMask = CycleMask( :,:, ( nPointsInterpPercentage < 1 ) );
CycleMask = CycleMask( :,:, : );

MaskedData.X = X;
MaskedData.Y = Y;
MaskedData.U = U;
MaskedData.V = V;
% MaskedData.U = InterpolatedData.U .* repmat( CycleMask, [ 1 1 1 nCycle ]);
% MaskedData.V = InterpolatedData.V .* repmat( CycleMask, [ 1 1 1 nCycle ]);
% MaskedData.CrankAngle = InterpolatedData.CrankAngle;


%% Prepare data format for POD use
PODData.X = cell( nCrankAngle_Using, 1 );
PODData.Y = cell( nCrankAngle_Using, 1 );
PODData.U = cell( nCrankAngle_Using, 1 );
PODData.V = cell( nCrankAngle_Using, 1 );
PODData.IndexInOriginal = cell( nCrankAngle_Using, 1 );

for ca_No = 1 : nCrankAngle_Using
    PODData.IndexInOriginal{ ca_No } = find( CycleMask( :, :, ca_No ) == 1 );
    temp_nValidPoints = length( PODData.IndexInOriginal{ ca_No } );
    
    PODData.X{ ca_No } = Convert_PODFormat( X, 'Original2POD', [], [], PODData.IndexInOriginal{ ca_No } );
    PODData.Y{ ca_No } = Convert_PODFormat( Y, 'Original2POD', [], [], PODData.IndexInOriginal{ ca_No } );
    
    PODData.U{ ca_No } = zeros( temp_nValidPoints, nCycle );
    PODData.V{ ca_No } = zeros( temp_nValidPoints, nCycle );
    for cycle_No = 1 : nCycle
        PODData.U{ ca_No }( :,cycle_No ) = Convert_PODFormat( MaskedData.U( :,:,ca_No,cycle_No ), 'Original2POD', [], [], PODData.IndexInOriginal{ ca_No } );
        PODData.V{ ca_No }( :,cycle_No ) = Convert_PODFormat( MaskedData.V( :,:,ca_No,cycle_No ), 'Original2POD', [], [], PODData.IndexInOriginal{ ca_No } );
    end
    
end
clear temp_* ca_No cycle_No

PODData.CrankAngle = MaskedData.CrankAngle;
PODData.nRowsInOriginal = size( X, 1 );
PODData.nColsInOriginal = size( X, 2 );

end

