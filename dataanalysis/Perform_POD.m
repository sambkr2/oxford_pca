function [ Result ] = Perform_POD( velo_data, POD_center_index, svd_method )
%Perform_POD performs proper orthogonal decomposition to centered (around 
% ensemble mean) data or not centered data. The result is found by either 
% directly performing svd to the data (assumed to be nLoc*nSnap),
% or by the method of snapshots (which gives the same results as
% directly perform svd to the data).
%
% Author(s): Li (Sam) Shen
% sam-li.shen@eng.ox.ac.uk
% Last updated date: 2020.04.28 (checked on 2022.04.26)
%
%   Required Inputs:
%       velo_data:              Velocity data, in complex form (i.e., u+i*v);
%                               Rows: Locations, Columns: Snapshots;
%                               assuming nLocations >= nSnapshots ('skinny' matrix).
%       POD_center_index:       'Centered'    = Center the data around the 
%                                               ensemble mean before conducting POD
%                               'NotCentered' = The data is not centered around
%                                               the ensemble mean (not recommended)
%       svd_method:             'Direct'      = Direct apply svd to the skinny matrix
%                               'Snapshot'    = Method of Snapshots, a correlation matrix
%                                               is used for svd calculation
%
%   Note: Using either svd_method options should give exactly the same results.
%
%   Available Outputs:
%       Result.nLoc:            No. of locations
%       Result.nSnap:           No. of snapshots
%       Result.nMode:           No. of modes
%       Result.svdInput:        svd input Velocity, Rows: Locations, Columns: Snapshots
%       Result.EnsembleMean:    Ensemble mean velocity, nLoc * 1
%       Result.C:               Corrlation matrix that is sent into svd
%       Result.svdU:            svd result U
%       Result.svdS:            svd result S
%       Result.svdV:            svd result V
%       Result.svdMode:         Modes based on svd result, before rotating
%       Result.svdCoeff:        Mode coefficients based on svd result, before rotating
%       Result.Mode:            Modes, Rows: Locations, Columns: Modes (in complex form)
%       Result.Coeff:           Mode coefficients, Rows: Modes, Columns: Snapshots
%       Result.KE.mode:     	Kinetic Energy for each mode, in J/kg (nMode*1)
%       Result.CenterIndex:     Record the inputted 'POD_center_index' 


%% Data formatting
% Find No. of locations and snapshots
[ nLoc, nSnap ] = size( velo_data );

% Check whether the data need to be centered
ensemble_mean = mean( velo_data, 2 );                                       % Ensemble mean, nLoc * 1
switch POD_center_index
    case 'Centered'
        velo_data_POD = velo_data - repmat( ensemble_mean, 1, nSnap );      % Ensemble mean is subtracted from each snapshot before conducting POD
    case 'NotCentered'
        velo_data_POD = velo_data;                                          % Use original data while conducting POD
        if ensemble_mean ~= 0
            warning( 'The input data for POD is not centered around ensemble mean, which is not recommended.' )
        end
    otherwise
        error( 'Wrong ''POD_center_index'' index. Please type ''help Perform_POD'' for detail.' )
end


% Format velocity matrix
svdInput = [ real( velo_data_POD ); imag( velo_data_POD ) ];                % Rows: Locations * 2, Columns: Snapshots


%% POD analysis
switch svd_method                                                           % Detect svd method
    case 'Direct'                                                           % Direct apply svd to the skinny matrix
        % svd
        C = svdInput;                                                       % Correlation matrix = Velocity (nLoc * nSnap)
        [ svdU, svdS, svdV ] = svd( C, 0 );                                 % C = U * S * V'
        nMode = rank( svdS );                                               % No. of modes
        % Eliminate extra zero rows/columns according to No. of modes
        svdU = svdU( :, 1 : nMode );
        svdS = svdS( 1 : nMode, 1 : nMode );
        svdV = svdV( :, 1 : nMode );
        % Calculate modes and coefficients
        KE_mode = diag( svdS ).^2 / 2;                                      % Kinetic energy for each mode
        Mode = svdU;                                                        % Modes, Rows: Locations, Columns: Modes
        Coeff = Mode' * svdInput;                                           % Mode coefficients, Rows: Modes, Columns: Snapshots; can be also obtained using svdS * svdV'
        
    case 'Snapshot'                                                         % Method of Snapshots, building a correlation matrix for svd calculation (not recommended)
        % Check "skinny" matrix assumption
        if nLoc < nSnap
            error( 'nLocations must be no less than nSnapshots while using snapshot method.' )
        end
        % svd
        C = svdInput' * svdInput;                                           % Correlation matrix = Velocity' * Velocity (nSnap * nSnap)
        [ svdU, svdS, svdV ] = svd( C, 0 );                                 % C = U * S * V' ( U = V )
        nMode = rank( svdS );                                               % No. of modes
        % Eliminate extra rows/columns according to No. of modes
        svdU = svdU( :, 1 : nMode );
        svdS = svdS( 1 : nMode, 1 : nMode );
        svdV = svdV( :, 1 : nMode );
        % Calculate modes and coefficients
        KE_mode = diag( svdS ) / 2;                                         % Kinetic energy for each mode
        Coeff = svdS.^( 1/2 ) * svdU';                                      % Mode coefficients from svd result, Rows: Modes, Columns: Snapshots
        Mode = svdInput * svdU * diag ( ( diag ( svdS ) ).^( -1/2 ) );      % Modes based on svd result, Rows: Locations, Columns: Modes
        
    otherwise
        error( 'Wrong ''svd_method'' index. Please type ''help Perform_POD'' for detail.' )
end

%% Post-processing
% Step 1: Re-format mode for quiver plot: Rows: Locations, Columns: Modes(in complex form)
Mode = complex( Mode( 1 : nLoc , : ), ...
            Mode( ( nLoc + 1 ) : end , : ) );                            % Express real number modes in complex form
        
% Step 2: Keep most of mode 1 coefficients to be positive
if length( find( Coeff( 1, : ) < 0 ) ) > 1/2 * nSnap
    Coeff( 1, : ) =  - Coeff( 1, : );
    Mode( :, 1 ) = - Mode( :, 1 );
end

%% Outputs
Result.nLoc = nLoc;                                                         % No. of locations
Result.nSnap = nSnap;                                                       % No. of snapshots
Result.nMode = nMode;                                                       % No. of modes
% Result.svdInput = svdInput;                                               % svd input Velocity, Rows: Locations, Columns: Snapshots
% Result.C = C;                                                               % Corrlation matrix that is sent into svd
% Result.svdU = svdU;                                                         % svd result U
% Result.svdS = svdS;                                                         % svd result S
% Result.svdV = svdV;                                                         % svd result V
% Result.svdMode = svdMode;                                                   % Modes based on svd result, before rotating
% Result.svdCoeff = svdCoeff;                                                 % Mode coefficients based on svd result, before rotating
Result.EnsembleMean = ensemble_mean;                                        % Ensemble mean velocity, nLoc * 1
Result.Mode = Mode;                                                         % Modes, Rows: Locations, Columns: Modes (in complex form)
Result.Coeff = Coeff;                                                       % Mode coefficients, Rows: Modes, Columns: Snapshots
Result.KE_mode = KE_mode;                                                   % Kinetic Energy for each mode, in J/kg (nMode*1)
Result.CenterIndex = POD_center_index;                                      % 'Centered' = Centered around the ensemble mean before conducting POD, 'NotCentered' = Not centered
end