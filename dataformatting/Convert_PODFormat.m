function [ OutputData, OutputDataType ] = Convert_PODFormat( InputData, DataConvertDirection, nRowsinOriginalData, nColsinOriginalData, ConvertIndex )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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

