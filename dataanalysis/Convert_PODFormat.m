function [ OutputData, OutputDataType ] = Convert_PODFormat( InputData, DataConvertDirection, nRowsinOriginalData, nColsinOriginalData, ConvertIndex )
%Convert_PODFormat converts the POD data grid (no NaNs) to the original
%grid (matrix that contains NaN out of FOV), or the reverse.

switch DataConvertDirection
    case 'Original2POD'
        OutputDataType = 'POD';
        OutputData = InputData( ConvertIndex );
        
    case 'POD2Original'
        OutputDataType = 'Original';
        OutputData = nan( nRowsinOriginalData, nColsinOriginalData );
        OutputData( ConvertIndex ) = InputData;
end
    
end

