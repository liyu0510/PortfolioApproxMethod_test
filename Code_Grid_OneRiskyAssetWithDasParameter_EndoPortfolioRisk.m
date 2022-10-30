%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！%
% Date: 20220813
% Constraints:
%   e     = pi+ki'+b-d'        % current period equity
%   e'    = Ri*ki'+Rf*b-Rd*d'  
%   delta = d'/(e+d')
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！%

clear all
close all
clc



%% Model parameters
Code_Para;

% discounting
beta = para.beta;

% define the financial intermediary parameters
omega = para.omega;
delta = para.delta;

% define the interest rate pass-through parameters
theta = para.theta;


% define the management parameters
zeta   = para.zeta;       % governs overall management efficiency
psi    = para.psi;    % governs management efficiency w.r.t wealth
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


%% Interest rates
R_f   = 1.02;
rf    = R_f-1;
R_d   = 1+rf*omega;

% redefine the augmented parameters
sigma2j_tilde  = sigma2j - zeta*psi*kappaf*theta*(muj-rf)*sigma2epsilon;
rho_tilde      = (rho*sqrt(sigma2i)*sqrt(sigma2j)-0.5*zeta*psi*kappad*theta*(mui-rf))/(sqrt(sigma2i)*sqrt(sigma2j_tilde));


%% Endogenous state variables
% a denotes endogenous states
% a = (e)

% grid min max
grid_min = 0;
grid_max = 2;

n_a   = [9];
num_a = length(n_a);

grid_a.e_grid  = linspace(grid_min,grid_max,n_a(1))';
grid_a.e_grid  = grid_a.e_grid(2:end);


% grid combined
state_key = zeros(num_a, length(grid_a.e_grid));

for i = 1:numel(fieldnames(grid_a))
    fn = fieldnames(grid_a);
    state_key(i,1:length(grid_a.(fn{i}))) =  grid_a.(fn{i});
end

n_grid = size(grid_a.e_grid,1);

N_a = n_grid;

grid_value_a = cell2mat(struct2cell(grid_a))';

grid_value_w = (1+delta/(1-delta))*grid_value_a;


%% Control variables
% d denotes control variables
% d = (ki+kj, b)

n_d   = [2*N_a+1,2*N_a+1];
num_d = length(n_d);
grid_size = grid_value_a(2)- grid_value_a(1);

grid_d.grid_ki  = (1+delta/(1-delta))*linspace(min(grid_value_a-grid_size), max(grid_value_a), n_d(1))';
grid_d.grid_b   = (1+delta/(1-delta))*linspace(min(grid_value_a-grid_size)*0.99, max(grid_value_a)*0.99, n_d(1))';

% length
control_key = zeros(num_d, max([length(grid_d.grid_ki), length(grid_d.grid_b)]));

for i = 1:numel(fieldnames(grid_d))
    fn = fieldnames(grid_d);
    control_key(i,1:length(grid_d.(fn{i}))) =  grid_d.(fn{i});
end

n_grid = size(grid_d.grid_ki,1);

grid2 = [];
for i = 1:length(grid_d.grid_b)
    grid2 = [grid2;[grid_d.grid_ki ones(n_grid,1)*grid_d.grid_b(i,1)]];
end


grid_value_d = grid2';

check_nonneg = grid_value_d(1,:) + grid_value_d(2,:);
% grid_value_d = grid_value_d_all;
grid_value_d = grid_value_d(:,find(check_nonneg < (1+delta/(1-delta))*grid_max));

s                       = sum(grid_value_d,1);
[sortedSums, sortOrder] = sort(s, 'Ascend');
grid_value_d            = grid_value_d(:, sortOrder);

n_grid = size(grid_value_d,2);

N_d = n_grid;

%% Current period utility
% consum gives the matrix for current period consumption
    % given this period endogenous state (w) ---> N_a
    % given the choice variable (p,b) ---> N_d
    % how much one gets to consume in this period?

alpha = 3;

% c = linspace(0,5,10);
% utility function: CARA/exponential function
% u(c)= 1-1/alpha*exp(-alpha*c)
% x = 1-1/alpha*exp(-alpha.*c);

% utility function: CRRA/power function
% u(c)=c^(1-alpha)/(1-alpha)-1
% x = (c.^(1-alpha))/(1-alpha)+1;


consum = zeros(N_a, N_d);

for c_a=1:N_a
    for c_d=1:N_d
        consum(c_a, c_d) = grid_value_a(1,c_a) - grid_value_d(1,c_d) - grid_value_d(2,c_d) + delta/(1-delta)*grid_value_a(1,c_a);
    end
end
consum(consum < 0) = 0;



%% Exogenous state variables

% set inital discretization parameters to pin down loop
% parameters of the process

N            = [5]';
random_draws = 1000;
method       = 1;   
N_z          = N; % number of exogenous states

% discretize the portfolio return
portfolio_discretized = cell(N_a, N_d);

tic;
tempcounter = 1;
for a_c=1:N_a
    for d_c=1:N_d
                % given a level of total wealth and consumption, 
                % there is an associated level of portfolio mean and variance
                    % (ki, kj) is computed directly
                    % b is obtained from budget constraint
                w = grid_value_w(a_c);
                pi = consum(a_c,d_c);
                [Ewprime,Varwprime,~,~] = Code_Approximate_Portfolio_func(para, w, pi, rf); 
                
                A0           = 1;
                A1           = Ewprime;
                A2           = 0;  % no persistence
                SIGMA        = Varwprime;
                
                [Pr_mat,Pr_mat_key,zbar] = fn_var_to_markov(A0,A1,A2,SIGMA,N_z,random_draws,method);
                
                portfolio_discretized{a_c,d_c} = {Pr_mat,Pr_mat_key,zbar};
                
                if rem(tempcounter,100)==0
                    disp(tempcounter)
                end
                tempcounter=tempcounter+1;
    end
end
toc;

%% Next period state variables
% Phi_aprime gives the matrix for next period state variable w
    % given next period exogenous state (Return on P) ---> N_z
    % given this period endogenous state (w) ---> N_a (deposit = delta*w)
    % given the choice variable (p,b) ---> N_d
    % what is the wealth for next period?


Phi_aprime = zeros(N_z, N_d, N_a);
    for j = 1:N_d
        for k = 1:N_a
        % For each given level of (w,pi), there is an associated level of
        % risks with the portfolio
        pfolio     = portfolio_discretized{k,j};
        Pr_mat     = pfolio{1};
        Pr_mat_key = pfolio{2};
            
            for i = 1:N_z 
                    % investment is less than total wealth = equity+deposit
                    Phi_aprime(i,j,k) = Pr_mat_key(1,i)*(grid_value_a(1,k) - consum(k,j))...
                    - R_d*grid_value_w(1,k);
            end
        end
    end    

openvar("Phi_aprime");




%% Loop for value function
currdist    = Inf;
Tolerance   = 10^(-3);


% define empty value and policy function
VKron = zeros(N_a,N_z);
PolicyIndexesKron=zeros(N_a,N_z); %indexes the optimal choice for d given rest of dimensions a,z

tic;
tempcounter = 1; 
while currdist>Tolerance

    VKronold=VKron;

    for c_a=1:N_a
        
        for c_z=1:N_z
           
            pfolio     = portfolio_discretized{c_a,c_d};
            Pr_mat     = pfolio{1};
            
            RHSpart2=zeros(N_d,1);
            for zprime_c=1:N_z
                if Pr_mat(c_z,zprime_c)~=0 %multilications of -Inf with 0 gives NaN, this replaces them with zeros (as the zeros come from the transition probabilites)
                %                     Phi_of_d=Phi_aprimeKron(:,c_a,c_z,zprime_c);
                %                     RHSpart2=RHSpart2+VKronold(Phi_of_d,zprime_c)*pi_z(c_z,zprime_c);
                    for c_d=1:N_d
                        if isempty(find(grid_value_a <= Phi_aprime(zprime_c,c_d,c_a), 1, 'last'))
                            RHSpart2(c_d) = 0;
                        else
                            locate_w_grid = find(grid_value_a <= Phi_aprime(zprime_c,c_d,c_a), 1, 'last');
                            RHSpart2(c_d) = RHSpart2(c_d)+VKronold(locate_w_grid,zprime_c)*Pr_mat(c_z,zprime_c);
                        end
                    end
                end
            end
            % entireRHS = (consum(c_a,:)'.^(1-alpha))/(1-alpha) + beta*RHSpart2; %d by 1
            entireRHS = 1-1/(alpha)*exp(-alpha.*consum(c_a,:)') + beta*RHSpart2; %d by 1
            
            % then maximizing d indexes
            [VKron(c_a,c_z),PolicyIndexesKron(c_a,c_z)]=max(entireRHS,[],1);
        end
    end

    VKrondist=reshape(VKron-VKronold,[N_a*N_z,1]); 
    VKrondist(isnan(VKrondist))=0;
    currdist=max(abs(VKrondist));
        if rem(tempcounter,100)==0
           disp(tempcounter)
           disp(currdist)
        end

    tempcounter=tempcounter+1;
end
toc;

openvar("PolicyIndexesKron");


%% Stationary Distribution

simoptions.parallel = 0;
Case2_Type          = 1;

StationaryDist=StationaryDist_Case2_EndoRisk(PolicyIndexesKron,Phi_aprime,Case2_Type,N_d,N_a,N_z,portfolio_discretized,simoptions,grid_value_a)


%% General Equilibrium

% call Dynare package to compute steady state

% compute aggregate deposit supply

% compute the difference between deposit supply and demand

% if larger than 0, increase R^f. If smaller than 0, decrease R^f

%%




