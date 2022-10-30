% Compute approximated portfolio solution given the parameters

function[Ewprime,Varwprime,ki,kj] = Code_Approximate_Portfolio_func(para,w_in,pi_in,rf_in)

% define the financial intermediary parameters
omega = para.omega;
delta = para.omega;

% define the interest rate pass-through parameters
theta = para.theta;


% define the management parameters
zeta   = para.zeta;         % governs overall management efficiency
psi    = para.psi;          % governs management efficiency w.r.t wealth
kappad = para.kappad;       % pledgeability of domestic assets, affects rho_tilde
kappaf = para.kappaf;       % pledgeability of foreign asset, affects sigma2j_tilde

% define the TFP shock parameters
mui     = para.mui;
muj     = para.muj;
sigma2i = para.sigma2i;
sigma2j = para.sigma2j;
rho     = para.rho;

% define the uncertainty shock parameter
sigma2epsilon = para.sigma2epsilon;

% input variables
rf = rf_in;    % risk-free rate

% redefine the augmented parameters
sigma2j_tilde  = sigma2j - zeta*psi*kappaf*theta*(muj-rf)*sigma2epsilon;
rho_tilde      = (rho*sqrt(sigma2i)*sqrt(sigma2j)-0.5*zeta*psi*kappad*theta*(mui-rf))/(sqrt(sigma2i)*sqrt(sigma2j_tilde));


%%
% wealth, consumption, portfolio
w  = w_in;     % total asset level w = e+d
pi = pi_in;    % consumption
p  = w-pi;     % total portfolio size

% parameter to iterate
alpha = 1;

currdist    = Inf;
Tolerance   = 10^(-5);

% loop for alpha_new

while currdist > Tolerance
    
alpha_new = alpha;
    

% define the uncertainty risks under wealth effect
sigma2epsilon_tilde = sigma2epsilon*(1-psi*(1-omega*delta)*w*rf);

% define the portfolio solution
ki = 1/(1-rho_tilde^2)*((mui-rf))/(2*alpha_new*sigma2i) ...
    - 1/(1-rho_tilde^2)*rho_tilde*(muj-rf)/(2*alpha_new*sqrt(sigma2i)*sqrt(sigma2j_tilde)) ...
    + 1/(1-rho_tilde^2)*rho_tilde*zeta*sigma2epsilon_tilde/(2*alpha_new*sqrt(sigma2i)*sqrt(sigma2j_tilde));
kj = 1/(1-rho_tilde^2)*((muj-rf))/(2*alpha_new*sigma2j_tilde) ...
    - 1/(1-rho_tilde^2)*rho_tilde*(mui-rf)/(2*alpha_new*sqrt(sigma2i)*sqrt(sigma2j_tilde)) ...
    - 1/(1-rho_tilde^2)*zeta*sigma2epsilon_tilde/(2*alpha_new*sqrt(sigma2i)*sqrt(sigma2j_tilde));

% define new portfolio moments
Ewprime   = (theta*(mui-rf))*ki + (theta*(muj-rf))*kj + (1+rf)*(p);
Varwprime = sigma2i*ki^2 + sigma2j_tilde*kj^2 + 2*rho_tilde*sqrt(sigma2i)*sqrt(sigma2j_tilde)*ki*kj+zeta*sigma2epsilon_tilde*kj;
Ewprimeratio = Ewprime/w;

% define new alpha
alpha_new_num = (theta*(muj-rf))/(2*rho_tilde*sqrt(sigma2i)*sqrt(sigma2j_tilde)) ...
    -(theta*(mui-rf))/(2*sqrt(sigma2i)*sqrt(sigma2j_tilde))+1+rf;

alpha_new_denom = ki*(1-(sqrt(sigma2i))/sqrt(sigma2j_tilde)) ...
    + kj*(sqrt(sigma2j_tilde)/(rho_tilde*sqrt(sigma2i))-rho_tilde) ...
    + zeta*sigma2epsilon_tilde/(2*rho_tilde*sqrt(sigma2i)*sqrt(sigma2j_tilde));

alpha = alpha_new_num/alpha_new_denom;

currdist = abs(alpha_new - alpha)
end

b=p-ki-kj;



