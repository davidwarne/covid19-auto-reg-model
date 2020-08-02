function [Data_s] = simuldata(Data,theta)
%% SIMULDATA
% This function simulates stochastic epidemic curve with latent infectous population
% and regulatory mechnisms. Model is a discrete state coninuous time Markov process
% simulated using a Tau-Leaping method.
%
% Parameters:
%
% Data - observed data, use to initialise the time series.
%
% theta - model parameters: 
%         alpha_0 = theta(1); % leakage
%         alpha = theta(2); % infection rate
%         beta = theta(3);         % case recovery rate
%         gamma = theta(4);        % case detection rate
%         delta = theta(5);        % case death rate
%         eta = theta(6);          % unobserved removal rate
%         n = theta(7);            % inibitive strength ... like a hill constant
%         kappa = theta(8);        % initial under-estimate factor
%         logp = theta(9);         % log of reporting effect
%
% Returns:
%     Data_s - simulated data using the model with given parameters.
%
% Author:
%     David J. Warne (david.warne@qut.edu.au)
%           School of Mathematical Sciences
%           Science and Engineering Faculty 
%           Queensland University of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract parameters
alpha_0 = theta(1)/Data.P; % leakage
alpha = theta(2)/Data.P; % infection rate
beta = theta(3);         % case recovery rate
gamma = theta(4);        % case detection rate
delta = theta(5);        % case death rate
eta = theta(6);          % unobserved removal rate
n = theta(7);            % inibitive strength ... like a hill constant
kappa = theta(8);        % initial under-estimate factor
pC = 10^theta(9);
pR = 10^theta(10);
pD = 10^theta(11);

% intialise from first 100 confirmed cases
J = find((Data.C + Data.R + Data.D) >= 100);
start = J(1);

%copy initial condition
Data_s.C = Data.C(start:end);
Data_s.R = Data.R(start:end);
Data_s.D = Data.D(start:end);
Data_s.P = Data.P;
Data_s.I = zeros(size(Data_s.C)); % for additional inference
Data_s.Ru = zeros(size(Data_s.C));

% clear evolution data except initial condition
Data_s.C(2:end) = 0;
Data_s.R(2:end) = 0;
Data_s.D(2:end) = 0;


% set up model
f = @(C,R,D,n,pC,pR,PD) 1/(1+(C/pC + R/pR + D/pD)^n);
model = struct();

% dimensions 
model.k = [alpha_0;alpha;beta;gamma;delta;eta;n;pC;pR;pD]; % parameters
model.M = 5; % number of interactions
model.N = 6; % number of sub-populations

% state update matrices (M x N)
model.nu_minus = [1,1,0,0,0,0; % S + I -> 2I rate alpha_0 + alpha f(C)
                  0,1,0,0,0,0; % I -> C rate \gamma
                  0,0,1,0,0,0; % C -> R rate \beta
                  0,0,1,0,0,0; % C -> D rate \delta
                  0,1,0,0,0,0];% I -> R^u rate \eta*\beta
model.nu_plus = [0,2,0,0,0,0;
                 0,0,1,0,0,0;
                 0,0,0,0,1,0;
                 0,0,0,1,0,0;
                 0,0,0,0,0,1];
model.nu = model.nu_plus - model.nu_minus;

% initial state
model.X0 = [Data_s.P - Data_s.C(1) - ceil(kappa*Data_s.C(1)) - Data_s.R(1)-Data_s.D(1);
           ceil(kappa*Data_s.C(1));
           Data_s.C(1);
           Data_s.D(1);
           Data_s.R(1);
           0];

% hazard functions
model.a = @(X,k) [(k(1) + k(2)*f(X(3),X(4),X(5),k(7),k(8),k(9),k(10)))*X(1)*X(2);
                  k(4)*X(2);
                  k(3)*X(3);
                  k(5)*X(3);
                  k(3)*k(6)*X(2)];

% number of simulated days
T = length(Data_s.C)-1;

% run simulation
tau = 1.0; % timestep (half day resolution) 
[Z,t] = TauLeapingMethod(model,T,tau);
for i =1:T
    [J] = find(t <= i);
    Data_s.C(1+i) = Z(3,J(end));
    Data_s.R(1+i) = Z(5,J(end));
    Data_s.D(1+i) = Z(4,J(end));
    Data_s.I(1+i) = Z(2,J(end));
    Data_s.Ru(1+i) = Z(6,J(end));
end
Data_s.I(1) = ceil(kappa*Data_s.C(1));
%append initial zeros (to correctly match full time-series)
Data_s.C = [zeros(start-1,1); Data_s.C];
Data_s.R = [zeros(start-1,1); Data_s.R];
Data_s.D = [zeros(start-1,1); Data_s.D];
Data_s.I = [zeros(start-1,1); Data_s.I];
Data_s.Ru = [zeros(start-1,1); Data_s.Ru];
