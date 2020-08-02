function [Z,t] = TauLeapingMethod(bcrn,T,tau)
%% The Tau Leaping Method
% An approximate stochastic simulation algorithm for a discrete-state 
% continuous-time Markov process.
%
% Inputs:
%    bcrn - a interaction network struct
%    T    - the end time of the simulation
%    tau   - the timestep
% Outputs:
%    Z    -  time series of copy number vectors
%    t    -  vector of times
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Science and Engineering Faculty
%         Queensland University of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise
Nt = floor(T/tau);
Z = zeros(length(bcrn.X0),Nt+1);
t = zeros(1,Nt+1);
Z(:,1) = bcrn.X0;
for i=1:Nt
    % compute propensities
    a = bcrn.a(Z(:,i),bcrn.k);
    % generate poisson variates
    Y = poissrnd(a*tau);
    % update copy numbers
    Z(:,i+1) = Z(:,i) + (bcrn.nu') * Y;
    Z(Z(:,i+1) < 0,i+1) = 0;
    t(i+1) = t(i) + tau;
end
