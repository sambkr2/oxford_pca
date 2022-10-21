%% calculate_WRI.m 

%% Code for calculating the weighted relevance index between two vector fields 

% Weighting suppresses bias to poor values of RI when comparing relatively
% small vectors

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

% output = A matrix of WRI values, with NaN values specified by the mask &
%          input vector fields. 


function [output] = calculate_WRI(a_u,a_v,b_u,b_v,mask,varargin)

    %% Check inputs
    % if weight power not specified, use 2 (weighted by square of magnitudes)
    if length(varargin) >= 1
        error('Too many inputs to calculate_WRI')
    end

    

    %% Step 1 - Calculate angular difference
    
    % Find angle between U and V for each vector of first field
    theta_aa=atan2(a_v(:,:),a_u(:,:));        

    % Find angle between U and V for each vector of second field
    theta_bb=atan2(b_v(:,:),b_u(:,:));        

    % Calculate the absolute difference between the two
    thetaDiff=abs(theta_aa-theta_bb);         

    %% Step 2 - Calculate weighting factor
    
    % Calculates magnitude of first velocity field
    V_Mag_aa=(a_u(:,:).^2+a_v(:,:).^2).^0.5;           
    
    % Calculates magnitude of second velocity field
    V_Mag_bb=(b_u(:,:).^2+b_v(:,:).^2).^0.5;          

    % Apply crank angle mask for both velocity fields
    V_Mag_aa= double( V_Mag_aa.*double(mask(:,:)));   
    V_Mag_bb= double( V_Mag_bb.*double(mask(:,:)));


    % calculate median velocity present in either field
    V_mag_med_aa = median(V_Mag_aa(:),'omitnan');
    V_mag_med_bb = median(V_Mag_bb(:),'omitnan');


    %% Step 3 - Apply weighting to alignment penalty
    
    % Weight by smallest of 2 vectors to 'forgive' misalignment if either vector is small
    penalty = (1/2)*(1-cos(thetaDiff)).*(V_Mag_aa./V_mag_med_aa).*(V_Mag_bb./V_mag_med_bb);
 
    % Applying the digital mask and assigning the output
    output = double(penalty.*double(mask(:,:)));   % Apply mask to output


end