function [corr_runica,...
    activations_runica, activations_icams, activations_icamf, ...
    super] = compareICA(channels,frames)
% Test the ability of different ICA algorithms to separate synthetic source
% signals. Function based on testica() from EEGLAB

%% Generate artificial data

% Default values for parameters
if nargin<2
    if nargin==0
        channels = 31;
    end
    frames = 1000;
end

% Default values
sources = channels; % independent components
exppow = -0.05;
shape = 1.2;

% Generate artificial super-Gaussian sources
exppowers = zeros(1,channels);
exppowers(1) = 1.0;
for s=1:sources
  exppowers(s) = exp(exppow*(s-1));
end

% Synthesize random source activations
super=randn(sources,frames).*(exppowers'*ones(1,frames));  
super=sign(super).*abs(super.^shape); % make super-Gaussian if shape > 1

% Mixing the simulated sources into channels
forward = randn(channels,sources);   % random forward mixing matrix
data = forward*super; % these are the simulated observed data

%% MCMC ICA
fprintf('Decomposing the resulting simulated data using MCMC ICA \n');
[weights_mcmc,timeseries_mcmc,backproj_mcmc] = MCMC4eeg(data,channels);

testid = inv(weights_mcmc)*forward(:,1:channels);
testresult = eyelike(testid);
figure; surf(testresult); view(-52,50); rotate3d;

%% Extended Infomax ICA Algorithm (EEGLAB)
fprintf('Decomposing the resulting simulated data using runica() \n');
[weights_runica,sphere_runica,~,~,~,~,activations_runica] = runica(data);

testid = weights_runica*sphere_runica*forward(:,1:channels);
testresult = eyelike(testid);
figure; surf(testresult); view(-52,50); rotate3d;

corr_runica = correlate(activations_runica,super);

%% The Molgedey and Schuster Algorithm (ICA:DTU Toolbox)
fprintf('Decomposing the resulting simulated data using icaMS() \n');
[activations_icams,weights_icams] = icaMS(data);

testid = inv(weights_icams)*forward(:,1:channels);
testresult = eyelike(testid);
figure; surf(testresult); view(-52,50); rotate3d;

%% The Mean Field ICA Algorithm (ICA:DTU Toolbox)
fprintf('Decomposing the resulting simulated data using icaMS() \n');
par = setupMFica(channels);
[activations_icamf,weights_icamf] = icaMF(data,par);

end

function [correlation] = correlate(activations,super)
% Test the correlation between the synthetic sources and the results of ICA
% algorithms.

    [corr,ind1,ind2] = matcorr(activations,super);
    fprintf('Absolute correlation between best-matching source and \n');
    fprintf('activation component pairs range from %5.4f to %5.4f\n', ...
        abs(corr(1)),abs(corr(length(corr))));
    correlation = [corr,ind1,ind2];

end

function [par] = setupMFica(channels)
% Set up MF ICA

    par.sources=channels;                  % number of sources
    par.optimizer = 'aem' ;                % optimizer
    par.solver = 'ec'; %variational' ;     % solver
    %par.method = 'positive' ;               % positive ICA
    % possible methods: {'positive','neg_kurtosis','pos_kurtosis','fa','ppca'}
    par.Sigmaprior = 'isotropic' ;     % noise variance
    par.S_tol = 1e-16 ;                % error tolerance E-step
    par.S_max_ite = 100 ;              % maximum number of iterations E-step                 
    par.tol = 1e-5;                    % error tolerance M-step for bfgs and conjgrad 
    par.max_ite = 50;                  % maximum number of iterations M-step
    par.draw = 1;                      % plot run time information

end