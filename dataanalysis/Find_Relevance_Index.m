function [ RI ] = Find_Relevance_Index( input_vector1, input_vector2, OutputType )
%Find_Relevance_Index finds the relevance index between two inputs. The
%inputs can be array or matrix with the same numbers of elements. The
%elements can either be real or complex number. If either of the elements
%in the two inputs is NaN, this element will be ignored.
%   Input(s):
%       input_vector1:  first input for comparison
%       input_vector2:  second input for comparsion, must have the same numbers
%                       of elements as vector1
%       OutputType:     the type for outputs
%                           'abs_str' = absolute value of relevance index as string with 4 digits, ranging from [ 0, 1 ]
%                           'normal_str' = value of relevance index as string with 4 digits, ranging from [ -1, 1 ]
%                           'abs_num' = absolute value of relevance index as a number, ranging from [ 0, 1 ]
%                           'normal_num' = value of relevance index as a number, ranging from [ -1, 1 ]
%   Output:
%       Relevance index in the specified format

%% Reshape inputs into a colmun array
reshaped_vector1 = input_vector1(:);
reshaped_vector2 = input_vector2(:);

%% Check numbers of elements
if length( reshaped_vector1 ) ~= length( reshaped_vector2 )
    error( 'The numbers of elements for two inputs must be the same.' )
end

%% Switch complex number elements into real numbers
if isreal( reshaped_vector1 ) ~= 1
    reformatted_vector1 = [  real( reshaped_vector1 ) ; imag( reshaped_vector1 ) ];
else
    reformatted_vector1 = reshaped_vector1;
end

if isreal( reshaped_vector2 ) ~= 1
    reformatted_vector2 = [  real( reshaped_vector2 ) ; imag( reshaped_vector2 ) ];
else
    reformatted_vector2 = reshaped_vector2;
end

%% Ignore NaN
vector1 = reformatted_vector1( ~isnan( reformatted_vector1 ) & ~isnan( reformatted_vector2 ) );
vector2 = reformatted_vector2( ~isnan( reformatted_vector1 ) & ~isnan( reformatted_vector2 ) );

%% Calculate relevance index
switch OutputType
    case 'abs_str'
        RI = abs( vector1' * vector2 ) / ( norm( vector1 ) * norm( vector2 ) );
        RI = num2str( RI, '%10.4f' );
    case 'normal_str'
        RI = ( vector1' * vector2 ) / ( norm( vector1 ) * norm( vector2 ) );
        RI = num2str( RI, '%10.4f' );
    case 'abs_num'
        RI = abs( vector1' * vector2 ) / ( norm( vector1 ) * norm( vector2 ) );
    case 'normal_num'
        RI = ( vector1' * vector2 ) / ( norm( vector1 ) * norm( vector2 ) );
    otherwise
        error( 'Illegal Output Type.' )  
end

end

