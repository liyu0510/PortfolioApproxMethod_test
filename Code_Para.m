para.beta = 0.99;

para.omega = 0.80;
para.delta = 0.70;

% define the interest rate pass-through parameters
para.theta = 1;

% define the management parameters
para.zeta   = 0.90;       % governs overall management efficiency
para.psi    = 0.55;    % governs management efficiency w.r.t wealth
para.kappad = 0;       % pledgeability of domestic assets, affects rho_tilde
para.kappaf = 0;       % pledgeability of foreign asset, affects sigma2j_tilde

% define the TFP shock parameters
para.mui     = 0.25;
para.muj     = 0.25;
para.sigma2i = 0.015;
para.sigma2j = 0.015;
para.rho     = -0.50;

% define the uncertainty shock parameter
para.sigma2epsilon = 1.5*para.sigma2j;

