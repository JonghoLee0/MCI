clear all



%% Parameters
trialN = 1000;                      % number of trials
d = 1;                              % scene depth(m)


% Lighting conditions
r_a = 1.0;                          % relative ambient light strength: r_a = e_a/e_s
r_i = 1.0;                          % relative interfering signal strength: r_i = e_i/e_s
N = 5;                              % number of interfering cameras


% ToF parameters
f_mod = 20e6;                       % modulation frequency(Hz)
T = 10e-3;                          % total integration time(s)


% Number of generated electrons
e_s = 1e7;                          % average number of signal photons
e_a = r_a*e_s;                      % average number of ambient photons
e_i = r_i*e_s;                      % average number of interfering photons


% Misc
c = 3e8;                            % Light speed(m/s)
d_max = c/(2*f_mod);                % Maximum measurable distance(m)


% Parameters for our approaches
A = 7;                              % Peak power amplification
M = 500;                            % Number of slots


% Parameters for PN
stageN = 7;                         % Number of LFSR stages
bitN = 2^stageN - 1;                % Number of bits
sampleNperBit = 1000;               % Number of samples per bit


% Photon energy
h = 6.62607015e-34;                 % Planck constant
lambda = 900e-9;                    % source wavelength
E_photon = h*c/lambda;              % unit photon energy



%% Repeat depth estimation and get depth standard deviations


% dummy variables for various parameters
dummyVariables = 1 : 1 : 10;            % N
% dummyVariables = [10 : 10 : 90]*1e-3;   % T
% dummyVariables = [10 : 10 : 90]*1e6;    % f_mod



% buffers for depth standard deviations
stdSet_OCA = zeros(size(dummyVariables, 2), 1);
stdSet_CRA = zeros(size(dummyVariables, 2), 1);



% loop over dummy variables
idx = 1;

for N = dummyVariables
% for T = dummyVariables
% for f_mod = dummyVariables
    
    idx
    
    
    % optimal slot ON probability
    p_ON = min(1/(2*N+1), 1/A);
    
    
    % adjust parameters for fair comparisons
    T_CRA = T/(A*p_ON);
    
    
    % buffers for estimated depths
    dSet_OCA = zeros(trialN, 1);
    dSet_CRA = zeros(trialN, 1);
    
    for trial = 1 : trialN 

    
        % ACO
        d_hat = estimateDepth_OCA(d, c, N, e_s, e_a, e_i, f_mod, T);
        dSet_OCA(trial, 1) = d_hat;
        

        
        % SEC_sync
        d_hat = estimateDepth_CRA(d, c, p_ON, N, M, A, e_s, e_a, e_i, f_mod, T_CRA);
        dSet_CRA(trial, 1) = d_hat;


    end
    
    
    % save depth standard deviations by simulations
    stdSet_OCA(idx, 1) = std(dSet_OCA);
    stdSet_CRA(idx, 1) = std(dSet_CRA);
    
           
    idx = idx + 1;
end



%% Comparisons (lower is better)
colors = [
  1 0.4 0
  0 0.8 0.2
];


figure; hold on; grid on;
plot(dummyVariables', stdSet_OCA, 'color', colors(1, :), 'lineWidth', 4);
plot(dummyVariables', stdSet_CRA, 'color', colors(2, :), 'lineWidth', 4);
ylabel('Depth standard deviation (m)')
xlabel('Number of interfering cameras')
xlim([dummyVariables(1), dummyVariables(end)])

legend('OCA', 'CRA')




