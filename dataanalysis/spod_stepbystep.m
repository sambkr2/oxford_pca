X = veloT;
%%
dim     = size(X);
nt      = dim(1);
nx      = prod(dim(2:end));
isrealx = isreal(X(1));
%%
window = []; weight = []; nOvlp = []; 
nDFT        = 2^floor(log2(nt/10));
N = nDFT;
window = 0.54-0.46*cos(2*pi*(0:N-1)/(N-1))';
window_name = 'Hamming';
nOvlp = floor(nDFT/2);
weight      = ones(nx,1);
weight_name = 'uniform';
nBlks = floor((nt-nOvlp)/(nDFT-nOvlp));

%%
winWeight   = 1/mean(window);
x_mean      = mean(X,1);
x_mean      = x_mean(:);
mean_name   = 'long-time (true) mean';

f = (0:nDFT-1)/dt/nDFT;    
f = (0:ceil(nDFT/2))/nDFT/dt;
nFreq = length(f);

%%
Q_blk       = zeros(nDFT,nx);
Q_hat = zeros(nFreq,nx,nBlks);
for iBlk    = 1:nBlks
        offset                  = min((iBlk-1)*(nDFT-nOvlp)+nDFT,nt)-nDFT;
        timeIdx                 = (1:nDFT) + offset;
        Q_blk   = bsxfun(@minus,X(timeIdx,:),x_mean.');
        Q_blk                   = bsxfun(@times,Q_blk,window);
        Q_blk_hat               = winWeight/nDFT*fft(Q_blk);
        Q_blk_hat               = Q_blk_hat(1:nFreq,:);
        Q_blk_hat(2:end-1,:)    = 2*Q_blk_hat(2:end-1,:);
        Q_hat(:,:,iBlk)         = Q_blk_hat;
end

%%
L   = zeros(nFreq,nBlks);
P   = zeros(nFreq,nx,nBlks);
for iFreq = 1:nFreq
    disp(['frequency ' num2str(iFreq) '/' num2str(nFreq) ' (f=' num2str(f(iFreq),'%.3g') ')'])
    Q_hat_f             = squeeze(Q_hat(iFreq,:,:));
    M                   = Q_hat_f'*bsxfun(@times,Q_hat_f,weight)/nBlks;
    [Theta,Lambda]      = eig(M);
    Lambda              = diag(Lambda);
    [Lambda,idx]        = sort(Lambda,'descend');
    Theta               = Theta(:,idx);
    Psi                 = Q_hat_f*Theta*diag(1./sqrt(Lambda)/sqrt(nBlks));
    % Theta               = Theta*diag(sqrt(Lambda)*sqrt(nBlks));
    P(iFreq,:,:)        = Psi;
    L(iFreq,:)          = abs(Lambda);
end
P   = reshape(P,[nFreq dim(2:end) nBlks]);

