function y = Rastrigin(x)
% 
% Rastrigin function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 2; 
s = 0;
for j = 1:n
    s = s+(x(j)^2-10*cos(2*pi*x(j))); 
end
y = 10*n+s;
%{ 
---------------------------------------------------------------------------
% RUN symDIRECT      
fun = 'Rastrigin'; x_L = [-5.0 -5.0]'; x_U = [6.0 6.0]'; 
GLOBAL.tolerance = 0.01; GLOBAL.optimalvalue = 0; PriLev = 2;
Result = symDIRECT(fun,x_L,x_U,GLOBAL,PriLev)
---------------------------------------------------------------------------  
% RUN Direct.m (Finkel implementation)
clear all;
% 1. Establish bounds for variables
bounds = [-5 6;-5 6];

% 2. Send options to Direct
% We tell DIRECT that the globalmin = 0.
% It will stop within 0.01% of solution
options.testflag  = 1; options.globalmin = 0.; 
options.showits   = 1; options.tol       = 0.01;

% 2a. NEW!
% Pass Function as part of a Matlab Structure
Problem.f = 'Rastrigin';

% 3. Call DIRECT
[fmin,xmin,hist] = Direct(Problem,bounds,options);

% 4. Plot iteration statistics
plot(hist(:,2),hist(:,3))
xlabel('Fcn Evals');
ylabel('f_{min}');
title('Iteration Statistics for Rastrigin test Function');
---------------------------------------------------------------------------
%}