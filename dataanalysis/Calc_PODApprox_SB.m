function [ PODApprox ] = Calc_PODApprox_SB( PODResult, nModes, CycleNo, nRowsInOriginal, nColsInOriginal, IndexInOriginal, ensembleMean )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

switch PODResult.CenterIndex
    case 'Centered'
        PODcumsum = PODResult.Mode( :, 1:nModes ) * PODResult.Coeff( 1:nModes, : ) + ensembleMean;        
    case 'NotCentered'
        PODcumsum = PODResult.Mode( :, 1:nModes ) * PODResult.Coeff( 1:nModes, : );   
end

PODApprox.U = nan( nRowsInOriginal, nColsInOriginal, 1, length( CycleNo ) );
PODApprox.V = nan( nRowsInOriginal, nColsInOriginal, 1, length( CycleNo ) );
for jj = 1 : length( CycleNo )
    PODApprox.U(:,:,1,jj) = Convert_PODFormat( real( PODcumsum(:,jj) ), 'POD2Original', nRowsInOriginal, nColsInOriginal, IndexInOriginal );
    PODApprox.V(:,:,1,jj) = Convert_PODFormat( imag( PODcumsum(:,jj) ), 'POD2Original', nRowsInOriginal, nColsInOriginal, IndexInOriginal );
end

end

