function figure_output = ColourQuiver_SB( x, y, u, v, figureprop )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%%
if ~exist( 'figureprop','var' ) || isempty( figureprop )
    warning( 'No figureprop defined, using default values.' )
    figureprop.axes_lim = [];
    figureprop.xlabel = [];
    figureprop.ylabel = [];
    figureprop.sparse_vector = [];
    figureprop.Clim = [];
else
    temp_figureprop_FieldsDetection = { 'axes_lim', 'xlabel', 'ylabel', 'sparse_vector', 'Clim', 'velocity_normalisation' };
    FieldNotExist = temp_figureprop_FieldsDetection( ~isfield( figureprop, temp_figureprop_FieldsDetection ) );
    for ii = 1 : length( FieldNotExist )
        warning( ['figureprop.',FieldNotExist{ii},' is not defined, using default value.'] )
        eval( [ 'figureprop.', FieldNotExist{ii}, '=[];' ] )
    end
end



%%
quiver_length = 1.5;

speedMap = abs( complex( u, v ) );
angleMap = angle( complex( u, v ) );



%%
if x(1,1) == x(1,2) && x(1,1) ~= x(2,1) && y(1,1) ~= y(1,2) && y(1,1) == y(2,1) % ndgrid
    x_imagesc = x(:,1).';
    y_imagesc = y(1,:).';
    speedMap_imagesc = speedMap.';

elseif x(1,1) ~= x(1,2) && x(1,1) == x(2,1) && y(1,1) == y(1,2) && y(1,1) ~= y(2,1) % meshgrid
    x_imagesc = x(1,:);
    y_imagesc = y(:,1);
    speedMap_imagesc = speedMap;
    
else
    error( 'Inputs x and y must be ndgrid or meshgrid, please check.' )
end


if isempty( figureprop.sparse_vector )
    sparse_vector = 1;
else
    if mod( figureprop.sparse_vector, 1 ) == 0 && figureprop.sparse_vector > 0
        sparse_vector = figureprop.sparse_vector; % Plot every n vectors
    else
        warning( 'figureprop.sparse_vector must be a postive integer, using default value 1.' )
        sparse_vector = 1;
    end
end

x_quiver = x( 1:sparse_vector:end, 1:sparse_vector:end );
y_quiver = y( 1:sparse_vector:end, 1:sparse_vector:end );
angleMap_quiver = angleMap( 1:sparse_vector:end, 1:sparse_vector:end );
speedMap_quiver = speedMap( 1:sparse_vector:end, 1:sparse_vector:end );


%%
% figure_output = figure( 'WindowState', 'maximized', 'Color', [ 1 1 1 ] );

figure_output = figure( 'Color', [ 1 1 1 ] );
hold on
box on

imagesc( x_imagesc, y_imagesc, speedMap_imagesc, 'AlphaData', ~isnan( speedMap_imagesc ) );
axis xy

if ~isempty( figureprop.Clim )
    set( gca, 'Clim', figureprop.Clim ) % set colorbar limits if user has specified
else
    set( gca, 'Clim', [ 0 round( max( speedMap_imagesc, [], 'all' )/10 )*10 ] ) % otherwise, start from zero and rounded to the nearest tens as maximum
end
temp_Clim = get(gca, 'Clim');

hcb = colorbar;


if isempty( figureprop.velocity_normalisation ) || ( figureprop.velocity_normalisation == 1 )
    hcb.Label.String = 'Speed (m/s)';
elseif figureprop.velocity_normalisation ~= 1
    hcb.Ticks = linspace( temp_Clim(1), temp_Clim(2), 6 ); 
    hcb.TickLabels = cellstr( num2str( ( 0:0.2:1 )' ) );
    hcb.Label.String = 'Speed (n.u.)';
end



% Plot all the vectors in white first
quiver( x_quiver, y_quiver, quiver_length * cos( angleMap_quiver ), quiver_length * sin( angleMap_quiver ), 0, 'Color', 'w', 'LineWidth', 1 ) 

% Find vectors that have high magnitude
HighSpeedCoordinates_x = ( x_quiver( speedMap_quiver >= 0.7*temp_Clim(2) ) );
HighSpeedCoordinates_y = ( y_quiver( speedMap_quiver >= 0.7*temp_Clim(2) ) );
HighSpeedCoordinates_angle = ( angleMap_quiver( speedMap_quiver >= 0.7*temp_Clim(2) ) );

% Plot the high magnitude vectors in black, in order to increase their
% visibility on the bright background
quiver( HighSpeedCoordinates_x, HighSpeedCoordinates_y, quiver_length * cos( HighSpeedCoordinates_angle ), quiver_length * sin( HighSpeedCoordinates_angle ), 0, 'Color', 'k', 'LineWidth', 1.5 ) 

axis equal


set( gca, 'FontSize', 22, 'FontName', 'Serif', 'xtick', -100:10:100, 'ytick', -100:10:100 );


if ~isempty( figureprop.axes_lim )
    axis( figureprop.axes_lim );
end

if ~isempty( figureprop.xlabel )
    xlabel( figureprop.xlabel );
end

if ~isempty( figureprop.ylabel )
    ylabel( figureprop.ylabel );
end


% QuiverGrid = [ x(:), y(:) ];
% 
% HighSpeedIndex
% 
% % HighRVBCoordinate = [ CFDData.RVB.Data.x_CFDGrid( temp_RVB >= 0.8*temp_Clim(2) ), CFDData.RVB.Data.y_CFDGrid( temp_RVB >= 0.8*temp_Clim(2) )];
% CFDVectorGrid = [temp_CFD_x(:),temp_CFD_y(:)];
% HighRVBIndex_in_CFDVectorGrid = dsearchn( CFDVectorGrid, HighRVBCoordinate );
% 
% CFDVector = [temp_CFD_u(:),temp_CFD_v(:)];
% 
% HighRVB_CFDVectorGrid =  CFDVectorGrid( HighRVBIndex_in_CFDVectorGrid, : ); % [ x y ]
% HighRVB_CFDVector = CFDVector( HighRVBIndex_in_CFDVectorGrid, : ); % [ u v ]
% 
% 
% quiver( HighRVB_CFDVectorGrid(:,1), HighRVB_CFDVectorGrid(:,2), HighRVB_CFDVector(:,1), HighRVB_CFDVector(:,2), 0, 'k', 'Linewidth', 1 )



end

