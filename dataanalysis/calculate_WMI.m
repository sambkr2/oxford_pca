%% calculate_WMI.m 

%% Code for calculating the weighted magnitude index between two vector fields

% Weights the absolute difference betwee 2 velocity fields by the median velocity magnitude in either field 

% Edited CW 20180206 & BS 20190108
% University of Oxford

% Expects to compare two 2D velocity fields (defined by u & v matrices)
% using a mask.

% Inputs: 

% a_u = 2D matrix of horizontal components (hence u) for vector field a.
% a_v = 2D matrix of vertical components (hence v) for vector field a.

% b_u = 2D matrix of horizontal components (hence u) for vector field b.
% b_v = 2D matrix of vertical components (hence v) for vector field b.

% mask = A mask to select a region of interest/validity within the overall
%        rectangular input matrices. Elements must be "1" for valid, "NaN" for
%        invalid. 


% Outputs: 

% output = A matrix of WMI values, with NaN values specified by the mask &
%          input vector fields. 


function [output]=calculate_WMI(a_u,a_v,b_u,b_v,mask,varargin)

    %% check inputs
    
    % Checking inputs...
    if length(varargin) >= 1
        error('Too many inputs to calculate_WMI')
    end
    
    %% calculate WMI
    
    % Calculates magnitude of first velocity field
    V_Mag_aa=(a_u(:,:).^2+a_v(:,:).^2).^0.5;       
    
    % Calculates magnitude of second velocity field
    V_Mag_bb=(b_u(:,:).^2+b_v(:,:).^2).^0.5;         

    % Apply crank angle mask to both velocity fields
    V_Mag_aa= double( V_Mag_aa.*double(mask(:,:))); 
    V_Mag_bb= double( V_Mag_bb.*double(mask(:,:)));
    
    % calculate median velocity present in either field to be used for
    % normalisation
    V_Mag_med = median([V_Mag_aa(:) ; V_Mag_bb(:)],'omitnan');
    
    % Calculating the difference on an absolute scale within the image and
    % then normalising by the median of either field
    VMagDiff=(abs((V_Mag_aa-V_Mag_bb))./V_Mag_med); 
    
    % Apply crank angle mask
    VMagDiff= double( VMagDiff.*double(mask(:,:)));  

    % Assigning the output
    output = VMagDiff;
    
end