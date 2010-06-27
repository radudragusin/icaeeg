function [Wr,Zr,WZ] = MCMC4eeg(V,ICs)
%
% MCMC4EEG
% takes as parameters: data (channels as rows, time series as columns)
%    ICs: number of independent components
% returns: weight matrix, activation time series, backprojected data
%
% Based on the algorithm from D.H. Hald, S.E. Grützmeier, S. Harder, C.I.
% Blücher: An ICA algorithm based on MCMC


%% Initialization

% Fix random generators for repetibility
rand( 'state', 0 );  %#ok
randn( 'state', 0 );  %#ok

X = V;
[d N] = size(X);

% Sampler structure
prs.m = ICs;          % Number of components
prs.chains = 1;     % Number of chains
prs.skip = 10;     % number of samples to burn
prs.nsamples = 15; % total number of steps
prs.stride = 1;     % stepsize
prs.std = 0;        % standardize with z-score
prs.ss = 100;       % parameter in gammadistribution for noise variance(c)
prs.sr = 1;         % parameter in gammadistribution for noise variance(d)
prs.uu = 1;         % 1 if you want to optimize the parameters of K
prs.upsilon = ones( prs.m, 1 ); % initial values for K
prs.jit = 1e-4;     % add to the diagonal of K for numerical stability
prs.tpts = (1:N)/N; % time points, it is normalized time right now
prs.us = 2;         % shape of the gamma distribution for upsilon
prs.ur = 100;       % rate of the gamma distribution for upsilon
prs.rr = 60;        % rate of the gamma distribution for tau
prs.rs = 10;        % shape of the gamma distribution for tau

%% MCMC algorithm

% array allocation
hlen = prs.nsamples*prs.chains;
hSigma = zeros(d, hlen);
hW = zeros(d*prs.m, hlen);
hupsilon = zeros(prs.m, hlen);
hZ = zeros(prs.m*N, hlen);

K = zeros(N, N, prs.m);
R = ones(d, prs.m);
accr = zeros(prs.m, 1);

% Calculate the squared distance
ds = sq_dist( prs.tpts, prs.tpts );

fprintf( 'starting loop\n' )

for r=1:prs.chains                                      % Chain loop
    
    % Set initial values
    W = 0.5*randn( d, prs.m );                              % Normal
    Sigma = 1./gamrnd( prs.ss, 1/prs.sr, [ d, 1 ] );        % Gamma
    Z = 0.5*randn( prs.m, N );                              % Normal
    tau = gamrnd( prs.rs, 1/prs.rr );                       % Gamma
    upsilon = prs.upsilon;
    D = ones([d prs.m]);
    
    % Creating the K matrix - use jitter for numerical stability
    for j=1:prs.m
        K(:,:,j) = exp( -upsilon(j)*ds ) + prs.jit*eye( N );
    end
    
    for k=1:prs.nsamples                                % Sample loop
        fprintf( '%d of %d\n',k,prs.nsamples)
        
        for l=1:prs.skip*( k == 1 ) + prs.stride*( k > 1 ) % Burn-in loop
            
% ========================= Gibbs for Z ===================================
            fprintf( 'Gibbs for Z\n') 
            for j = randperm(prs.m)       
                
                Z(j,:) = 0;
                e = X - W*Z;
                
                % Mean and covariance calculated from the Bayesian
                % framework
                
                xs = W(:,j)./Sigma;
                cc = xs'*W(:,j);
                
                % Cholesky  decomposition
                L = chol( K(:,:,j) + eye( N )/cc, 'lower' );
                V = L\K(:,:,j);
                Sj = K(:,:,j) - V'*V;
                L = chol( Sj, 'lower' );                
                
                mj = Sj*e'*xs;
                
                Z(j,:) = mj + L*randn( N, 1 );
            end
            
% ========================== Gibbs for W ==================================
            e = X - W*Z;
            fprintf( 'Gibbs for W\n')
            for j = randperm(prs.m)
                e = e + W(:,j)*Z(j,:);
                B = 1./( Z(j,:)*Z(j,:)' + D(:,j).*Sigma );
                mj = B.*( Z(j,:)*e' )';
                Sj = B.*Sigma;
                W(:,j) = normrnd( mj, sqrt( Sj ) ).*R(:,j);
                e = e - W(:,j)*Z(j,:);
            end
            
% ========================= Gibbs for Sigma ===============================
            D = ones(d, prs.m);
            e = X - W*Z;
            fprintf( 'Gibbs for Sigma\n')
            % sample from conditional for Sigma
            for i=1:d
                Sigma(i) = 1/gamrnd( prs.ss + 0.5*N, ...
                    1/( prs.sr + 0.5*( e(i,:)*e(i,:)' ) ) );
            end

% ============================ Update K ===================================
            if prs.uu
                fprintf( 'Update K\n')
                % Gamma distribution for upsilon
                upsilonp = gamrnd( prs.us, prs.ur, [ prs.m 1 ] ); 
                
                for j=1:prs.m
                    % New K
                    Kp = exp( -upsilonp(j)*ds ) + prs.jit*eye( N );
                    
                    num = log( mvnpdf( Z(j,:)', zeros( N, 1 ), Kp ) ); %New
                    den = log( mvnpdf( Z(j,:)', zeros( N, 1 ), K ) );  %Old
                    
                    ratio = exp( num - den );
                    
                    if rand() < min( ratio, 1 )     % Metropolis
                        accr(j) = accr(j) + 1;      % accuracy 
                        upsilon(j) = upsilonp(j);   % Update upsilon
                        K(:,:,j) = Kp;              % Update K
                    end 
                end 
            end
            
            
        end % end burn-in loop
        
        idx = k + ( r - 1 )*prs.nsamples;
        
        % save
        hSigma(:,idx) = Sigma;
        hW(:,idx) = W(:);
        hZ(:,idx) = Z(:);
        hupsilon(:,idx) = upsilon;
        
    end % end sample loop
end % end chain loop

%% Return results

fprintf( 'Done.\n' )
result.hSigma = hSigma;
result.hW = hW;
result.hZ = hZ;
result.hupsilon = hupsilon;
result.accr = accr;
save('MCMC4eeg_result.mat','prs','result');

%% Sort

N = size(result.hZ,1)/prs.m;
Z = reshape( median( result.hZ, 2 ), prs.m, N );
Wr = reshape( median( result.hW, 2 ), d, prs.m);

% Standardizing
Energy = zeros(1,prs.m);
for j=1:prs.m,
    amplZ = std(Z(j,:));                % standard deviation
    signZ = sign(mean(Z(j,:)));         % scale with the sign of Z
    Z(j,:) = signZ*Z(j,:)/amplZ;
    Wr(:,j) = signZ*Wr(:,j)*amplZ;
    Energy(j) = sum(Wr(:,j).*Wr(:,j));  % energy of W
end

% Sorting (arrange components by energy)
[Energy,idx] = sort(Energy,'descend');
SortZ = zeros(size(Z));
SortW = zeros(size(Wr));
for j = 1:prs.m
    SortZ(j,:) = Z(idx(j),:);
    SortW(:,j) = Wr(:,idx(j));
end

Zr = SortZ; % components activation (time course)
Wr = SortW; % components scalp map
WZ = Wr*Zr; % Project the components back
Wr = pinv(SortW);
