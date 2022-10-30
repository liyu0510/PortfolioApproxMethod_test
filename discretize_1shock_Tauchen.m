function [a_grid,P]=discretize_1shock_Tauchen(mew,rho,sigma,znum,Tauchen_q, tauchenoptions)
% Create states vector, a_grid, and transition matrix, P, for the discrete markov process approximation 
%    AR(1) process z'=mew(1-rho)+rho*z+e, e~N(0,sigma^2), by Tauchen method
%
% Inputs
%   mew(1-rho)     - constant term coefficient, corrected for a mean of mew
%   rho            - autocorrelation coefficient
%   sigma          - standard deviation of (gaussion) innovations
%   znum           - number of states in discretization of z (must be an odd number)
%   Tauchen_q      - (Hyperparameter) Defines max/min grid points as mew+-nSigmas*sigmaz (I suggest 2 or 3)
% Optional Inputs (tauchenoptions)
%   parallel:      - set equal to 2 to use GPU, 0 to use CPU
%   dshift:        - allows approximating 'trend-reverting' process around a deterministic trend (not part of standard Tauchen method)
% Outputs
%   z_grid         - column vector containing the znum states of the discrete approximation of z
%   P              - transition matrix of the discrete approximation of z;
%                    transmatrix(i,j) is the probability of transitioning from state i to state j
% Helpful info:
%   Var(z)=(sigma^2)/(1-rho^2). So sigmaz=sigma/sqrt(1-rho^2);   sigma=sigmaz*sqrt(1-rho^2)
%                                  where sigmaz= standard deviation of z
%     E(z)=mew/(1-rho)
%%%%%%%%%%%%%%%
% Original paper:
% Tauchen (1986) - "Finite state Markov-chain approximations to univariate and vector autoregressions"
% Original code: Kirkby Toolbox

















