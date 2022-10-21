%% Set-up
gamma_grd = 200;

% gammaval from flow type 3 (cylinder PIV)
gammaval = logspace(1,log10(1600),gamma_grd);

% set-up snapshot matrices
velo = [ DMDData.U; DMDData.V ];
X1 = velo(:,1:end-1);
X2 = velo(:,2:end);

% SVD
[U,S,V] = svd(X1,'econ');
% r = rank(S);
% 
% % Truncated versions of U, S, and V
% U = U(:,1:r);
% S = S(1:r,1:r);
% V = V(:,1:r);

% Determine matrix UstarX1
UstarX1 = U'*X1;

%% run_dmdsp
% Matrix Vstar 
Vstar = V';

% The number of snapshots
N = size(Vstar,2);

% Optimal DMD matrix resulting from Schmid's 2010 algorithm
% Fdmd = U'*X1*V*inv(S)
Fdmd = (UstarX1*V)/S;
% Determine the rank of Fdmd
r = rank(Fdmd); % set the number of modes
% E-value decomposition of Fdmd 
[Ydmd,Ddmd] = eig(Fdmd); 
Edmd = diag(Ddmd); % e-values of the discrete-time system

% Form Vandermonde matrix
Vand = zeros(r,N);
zdmd = Edmd;

for i = 1:N
    
    Vand(:,i) = zdmd.^(i-1);
    
end

% Determine optimal vector of amplitudes xdmd 
% Objective: minimize the least-squares deviation between 
% the matrix of snapshots X0 and the linear combination of the dmd modes
% Can be formulated as:
% minimize || G - L*diag(xdmd)*R ||_F^2
L = Ydmd;
R = Vand;
G = S*Vstar;

% Form matrix P, vector q, and scalar s
% J = x'*P*x - q'*x - x'*q + s
% x - optimization variable (i.e., the unknown vector of amplitudes)
P = (L'*L).*conj(R*R');
q = conj(diag(R*G'*L));
s = trace(G'*G);

% My P isn't quite semi-definite (rounding error), so add small number to 
% diagonal of P before calling Cholesky:
R = chol(P + 1e-11*eye(size(P)),'lower');
xdmd = (R')\(R\q);

% % Cholesky factorization of P
% Pl = chol(P,'lower');
% 
% % Optimal vector of amplitudes xdmd
% xdmd = (Pl')\(Pl\q);

% Sparsity-promoting algorithm starts here

% Set options for dmdsp
options = struct('rho',1,'maxiter',10000,'eps_abs',1.e-6,'eps_rel',1.e-4);

%% SP_DMD
answer = dmdsp(P,q,s,gammaval,options);

%% DMD spectrum
% rpm
rpm = 1500;

% time between snapshots. *2 factor accounts for 4 stroke engine (seconds
% per cycle)
dt = 1/(rpm/60)*2;

% get DMD frequency
dmdFreq = log(diag(L2))/dt/2/pi;

% get amplitude
amp = abs(b2(:));

ampsp = abs(xdmd);

% scatter b2
figure1 = figure( 'Color', [ 1 1 1 ] );
axes1 = axes( 'Parent', figure1, 'FontName', 'Times', 'FontSize', 24,'linewidth',1);
hold( axes1, 'on' )
box( axes1, 'on')
semilogy(imag(dmdFreq), amp,'o','MarkerFaceColor','blue','MarkerSize',8)
% semilogy(imag(dmdFreq(11)), amp(11),'o','MarkerFaceColor','red','MarkerSize',8)
% semilogy(imag(dmdFreq(14)), amp(14),'o','MarkerFaceColor','red','MarkerSize',8)
xlim([0 6.5])
xlabel('Frequency (Hz)')
ylabel('DMD mode amplitude')
% title({'Raw spectrum:','cross-tumble plane, -285 CAD'})
% title({'De-noised spectrum:','cross-tumble plane, -285 CAD'})

% scatter xdmd
figure2 = figure( 'Color', [ 1 1 1 ] );
axes2 = axes( 'Parent', figure2, 'FontName', 'Times', 'FontSize', 24,'linewidth',1);
hold( axes2, 'on' )
box( axes2, 'on')
semilogy(imag(dmdFreq), ampsp,'o','MarkerFaceColor','blue','MarkerSize',8)
% semilogy(imag(dmdFreq(11)), amp(11),'o','MarkerFaceColor','red','MarkerSize',8)
% semilogy(imag(dmdFreq(14)), amp(14),'o','MarkerFaceColor','red','MarkerSize',8)
xlim([0 6.5])
xlabel('Frequency (Hz)')
ylabel('DMD mode amplitude')
% title({'Raw spectrum:','cross-tumble plane, -285 CAD'})
% title({'De-noised spectrum:','cross-tumble plane, -285 CAD'})

% % scatter raw
figure
semilogy(imag(dmdFreq), amp,'o','MarkerFaceColor','blue')
% xlim([0 1])

figure
semilogy(imag(dmdFreq), ampsp,'o','MarkerFaceColor','blue')
% xlim([0 1])

%% Reconstruction plot parameters
pivgrid = CFDData.Data.y_PIVGrid;

figureprop.axes_lim = [ -25 25 -20 2 ];
figureprop.xlabel = '{\it y} (mm)';
figureprop.ylabel = '{\it z} (mm)';
figureprop.velocity_normalisation = 5;
figureprop.sparse_vector = 2;
figureprop.Clim = [ 0 50 ];

zeroind = imag(dmdFreq) <= 0.1 & imag(dmdFreq) >= 0;
zeroindf = find(zeroind == 1);

%% Inspect at 0 Hz freq

for jj = 1:length(zeroindf)
    
    % DMD mode * DMD eigenvalue * amplitude
%     recon = P2(:,posf_ind(kk))*L2(posf_ind(kk),posf_ind(kk))*b2(posf_ind(kk),1);
%     recon = P2(:,zeroindf(jj))*L2(zeroindf(jj),zeroindf(jj))*amp(zeroindf(jj),1);
    recon = P2(:,zeroindf(jj))*-amp(zeroindf(jj),1);
    
    % Separate U and V velocities (real)
    locs = size(recon);
    Ucol = real(recon(1:locs/2));
    Vcol = real(recon(locs/2+1:end));

    % Separate U and V velocities (imag)
    locs = size(recon);
    Ucoli = imag(recon(1:locs/2));
    Vcoli = imag(recon(locs/2+1:end));
    
    % Reshape into PIV grid format
    pivx = size(pivgrid,1);
    pivz = size(pivgrid,2);
    dmdplotu = NaN(pivx,pivz);
    dmdplotv = NaN(pivx,pivz);
    dmdplotu(originalIndex) = Ucol;
    dmdplotv(originalIndex) = Vcol;

    % Reshape into PIV grid format (imag)
    dmdplotui = NaN(pivx,pivz);
    dmdplotvi = NaN(pivx,pivz);
    dmdplotui(originalIndex) = Ucoli;
    dmdplotvi(originalIndex) = Vcoli;
    
    figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    
%     t = {['DMD mode ', num2str(zeroindf(jj)), ', ', num2str(round(imag(dmdFreq(zeroindf(jj))))),' Hz']};
%     title(t)
    figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotui, dmdplotvi, figureprop );    

%     set(gcf, 'position', [440 377 560 290])
%     dmdname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/DMD',num2str(zeroindf(jj)),'.png'];
%     exportgraphics(gcf,dmdname,'resolution',600)

end

%% SP-DMD

for jj = 1:length(zeroindf)
    
    % DMD mode * DMD eigenvalue * amplitude
%     recon = P2(:,posf_ind(kk))*L2(posf_ind(kk),posf_ind(kk))*b2(posf_ind(kk),1);
%     recon = P2(:,zeroindf(jj))*L2(zeroindf(jj),zeroindf(jj))*amp(zeroindf(jj),1);
    recon = P2(:,zeroindf(jj))*ampsp(zeroindf(jj),1)*L2(zeroindf(jj),zeroindf(jj));
    
    % Separate U and V velocities (real)
    locs = size(recon);
    Ucol = real(recon(1:locs/2));
    Vcol = real(recon(locs/2+1:end));

    % Separate U and V velocities (imag)
    locs = size(recon);
    Ucoli = imag(recon(1:locs/2));
    Vcoli = imag(recon(locs/2+1:end));
    
    % Reshape into PIV grid format
    pivx = size(pivgrid,1);
    pivz = size(pivgrid,2);
    dmdplotu = NaN(pivx,pivz);
    dmdplotv = NaN(pivx,pivz);
    dmdplotu(originalIndex) = Ucol;
    dmdplotv(originalIndex) = Vcol;

    % Reshape into PIV grid format (imag)
    dmdplotui = NaN(pivx,pivz);
    dmdplotvi = NaN(pivx,pivz);
    dmdplotui(originalIndex) = Ucoli;
    dmdplotvi(originalIndex) = Vcoli;
    
    figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu, dmdplotv, figureprop );    
%     t = {['DMD mode ', num2str(zeroindf(jj)), ', ', num2str(round(imag(dmdFreq(zeroindf(jj))))),' Hz']};
%     title(t)
%     figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotui, dmdplotvi, figureprop );    

%     set(gcf, 'position', [440 377 560 290])
%     dmdname = ['/Users/sambaker/Documents/Oxford-Uni/Papers/dmd/fig/CTP/DMD',num2str(zeroindf(jj)),'.png'];
%     exportgraphics(gcf,dmdname,'resolution',600)

end

%% Plot specific DMD modes

index = [11 12];
% first = P2(:,index)*L2(index,index)*b2(index,1);
% first = P2(:,index)*L2(index,index)*-amp(index,1);
first = P2(:,index)*L2(index,index)*ampsp(index,1);

% Separate U and V velocities (real)
Ucol = real(first(1:nValidPoints));
Vcol = real(first(nValidPoints+1:end));

% Separate U and V velocities (imag)
Ucoli = imag(first(1:nValidPoints));
Vcoli = imag(first(nValidPoints+1:end));

% Reshape into PIV grid format (real)
pivx = size(pivgrid,1);
pivz = size(pivgrid,2);
dmdplotu = NaN(pivx,pivz);
dmdplotv = NaN(pivx,pivz);
dmdplotu(originalIndex) = Ucol;
dmdplotv(originalIndex) = Vcol;

% Reshape into PIV grid format (imag)
dmdplotui = NaN(pivx,pivz);
dmdplotvi = NaN(pivx,pivz);
dmdplotui(originalIndex) = Ucoli;
dmdplotvi(originalIndex) = Vcoli;

figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotu/50, dmdplotv/50, figureprop );    

% figure_output = ColourQuiver( MaskedData.X, MaskedData.Y, dmdplotui, dmdplotvi, figureprop );    

% ylim([-30 10])
% t = {['DMD mode ', num2str(index), ', 0 Hz']};
% title(t)
set(gcf, 'position', [440 377 560 290])
