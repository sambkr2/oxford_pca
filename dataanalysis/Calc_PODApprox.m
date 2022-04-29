function [ PODApprox ] = Calc_PODApprox( PODResult, nModes, CycleNo )
%Calc_PODApprox calculates the POD approximated results at a given order
%(number of POD modes).

switch PODResult.CenterIndex
    case 'Centered'
        if mod( nModes, 1 ) == 0 && nModes > 0 && nModes <= PODResult.nMode
            PODcumsum = PODResult.Mode( :, 1:nModes ) * PODResult.Coeff( 1:nModes, CycleNo ) + repmat( PODResult.EnsembleMean, 1, length( CycleNo ) );
        elseif nModes == 0 % only ensemble mean
            PODcumsum = PODResult.EnsembleMean;
            CycleNo = 1;
        else
            error( 'Invalid nModes' )
        end
    case 'NotCentered'
        if mod( nModes, 1 ) == 0 && nModes > 0 && nModes <= PODResult.nMode
            PODcumsum = PODResult.Mode( :, 1:nModes ) * PODResult.Coeff( 1:nModes, CycleNo );
        elseif nModes == 0
            error( 'nModes = 0 is not supported when the POD is performed on a non-centered data.' )
        else
            error( 'Invalid nModes' )
        end
end

PODApprox.U = nan( PODResult.nRowsInOriginalGrid, PODResult.nColsInOriginalGrid, 1, length( CycleNo ) );
PODApprox.V = nan( PODResult.nRowsInOriginalGrid, PODResult.nColsInOriginalGrid, 1, length( CycleNo ) );
for jj = 1 : length( CycleNo )
    PODApprox.U(:,:,1,jj) = Convert_PODFormat( real( PODcumsum(:,jj) ), 'POD2Original', PODResult.nRowsInOriginalGrid, PODResult.nColsInOriginalGrid, PODResult.IndexInOriginalGrid );
    PODApprox.V(:,:,1,jj) = Convert_PODFormat( imag( PODcumsum(:,jj) ), 'POD2Original', PODResult.nRowsInOriginalGrid, PODResult.nColsInOriginalGrid, PODResult.IndexInOriginalGrid );
end

PODApprox.X = PODResult.X_OriginalGrid;
PODApprox.Y = PODResult.Y_OriginalGrid;

end

