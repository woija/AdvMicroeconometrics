%  Week1xa.m: Ex ante code for Week 1
%  Purpose: To introduce Gauss and to study the properties of OLS estimator. 
%  Agenda:
%       EXERCISE 1: Generate a synthetic data set.
%       EXERCISE 2: Estimate the parameters of a linear regression model by OLS.
%       EXERCISE 3: Carry out a Monte Carlo Study to verify the unbiasedness of the OLS Estimator.


%*********************************************
%* EXERCISE 1: GENERATE A SYNTHETIC DATA SET *
%*********************************************

% CLEAR WORKSPACE; CLEAR COMMAND WINDOW 
clear all; clc;

% GENERATE DATA
my_seed_variable = 1;                                   % set seed to 1
rng(my_seed_variable);                      % seed the random number generator --- produces same sequence every time. 
n=100;                                      % define n=100 - size of data set    
b=[1;-0.5;2];                               % define parameter values            
x1=normrnd(0,2,n,1);                        % generate x1~N(0,4) -> nx1          
x2=normrnd(5,3,n,1);                        % generate x2~N(5,9) -> nx1          
eps=normrnd(0,1,n,1);                       % generate eps~N(0,1)-> nx1          
x=[ones(n,1),x1,x2];                        % collect the expl. variables      
y= b*x + eps;                                % generates y -> nx1... recall that y = b_0 + b_1*x_1 + b_2*x_2 + b_3*x*3 + eps.


%*****************************************************************************
%  EXERCISE 2: ESTIMATE THE PARAMETERS OF A LINEAR REGRESSION MODEL BY OLS   *    
%*****************************************************************************

% ESTIMATION
b_hat  = ;                         % estimate kx1 vector of OLS parameters            
res    = y-x*b_hat;                         % OLS residuals                                    
sigma  = [FILL IN];                         % calculate sigma^2                                
b_var  = [FILL IN];                         % calculate cover. matrix of b_ols (kxk)           
b_stde = sqrt(diag(b_var));                 % extract standard errors from covar-matrix        
t      = b_hat./b_stde;                     % calculate t-values                               

% PRINT RESULTS
disp('           OLS-estimates                  ');
disp('==========================================');
disp('      b       b_hat    stde_b    t-value  ');
disp([b, b_hat, b_stde, t]);
disp('==========================================');



%******************************************
%*  EXERCISE 3: SIMPLE MONTE CARLO STUDY  *    
%******************************************
s=200;                                      % # replications                                      
n=100;                                      % # observations (as before)         
m1=zeros(s,size(b,1));                      % (Sx3) vector to save b estimates     
m2=zeros(s,size(b,1));                      % (Sx3) vector to save var estimates     
for i=1:s                                   % start loop                         

    % GENERATE DATA
    x1=normrnd(0,2,n,1);                    % generate x1~N(0,4) -> nx1  
    x2=normrnd(5,3,n,1);                    % generate x2~N(5,9) -> nx1  
    eps=normrnd(0,1,n,1);                   % generate eps~N(0,1)-> nx1  
    x=[FILL IN];                            % expl. variables    -> nx3  
    y=[FILL IN];                            % dep. variable      -> nx1      

    % ESTIMATION
    b_hat=[FILL IN];                        % OLS estimates -> 3x1  
    res=y-x*b_hat;                          % OLS residuals -> nx1  
    var_b=[FILL IN];                        % Variance-Covariance matrix of beta -> 3x3 
    
    % STORE ESTIMATES IN m1 AND m2
    m1(i,:)=b_hat;                          % put OLS estimate i into row i of m1         
    m2(i,:)=diag(var_b);                    % put estimate i of est. var into row i of m2 
 
end;                                        % end loop  

% DO MEAN AND STANDARD ERRORS
E_b_hat=mean(m1,1);                                 % Mean of OLS estimates                       
E_se_b=[FILL IN];                                   % Mean of estimated standar errors            
Mcse=sqrt(1/(s-1)*sum((m1-ones(s,1)*mean(m1,1)).^2,1));% Standard deviation of Monte Carlo Estimates 
% alternatively use the built-in function 
%Mcse = std(m1);

% PRINT RESULTS
disp('');
disp('           Monte Carlo Experiment                 ');
disp('==================================================');
disp('Mean estimates   Mean est. s.e.   Monte Carlo s.e.');
disp([E_b_hat', E_se_b', Mcse']);
disp('--------------------------------------------------');
disp(['Number of obs:          ' num2str(n)]);
disp(['Number of replications: ' num2str(s)]);
disp('==================================================');

% HISTOGRAM
hist(m1(:,2),20);                                           % hist(200x1 vector est., number of bins) 


