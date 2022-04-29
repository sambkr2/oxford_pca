function [ u_interp, v_interp ] = Interpolate_PIVVectors( x_data, y_data, u_data, v_data, x_interp, y_interp )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


InterpObjectU = scatteredInterpolant( x_data, y_data, u_data, 'natural', 'none' );
u_interp = InterpObjectU( x_interp, y_interp );

InterpObjectV = scatteredInterpolant( x_data, y_data, v_data, 'natural', 'none' );
v_interp = InterpObjectV( x_interp, y_interp );

end

