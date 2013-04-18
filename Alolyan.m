function f = Alolyan(x)
f = x(1)*x(2)^2 + x(2)*x(1)^2 - x(1)^3 - x(2)^3;
%{ 
---------------------------------------------------------------------------
% RUN symDIRECT   
fun = 'Alolyan'; x_L = [-1 -1]'; x_U = [1 1]'; 
GLOBAL.tolerance = 0.01; GLOBAL.optimalvalue = -1.18519; PriLev = 2;
Result = symDIRECT(fun,x_L,x_U,GLOBAL,PriLev)
---------------------------------------------------------------------------  
% RUN Direct.m (Finkel implementation)
clear all;
% 1. Establish bounds for variables
bounds = [-1 1;-1 1];

% 2. Send options to Direct
% We tell DIRECT that the globalmin = -1.18519
% It will stop within 0.01% of solution
options.testflag  = 1; options.globalmin = -1.18519; 
options.showits   = 1; options.tol       = 0.01;

% 2a. NEW!
% Pass Function as part of a Matlab Structure
Problem.f = 'Alolyan';

% 3. Call DIRECT
[fmin,xmin,hist] = Direct(Problem,bounds,options);

% 4. Plot iteration statistics
plot(hist(:,2),hist(:,3))
xlabel('Fcn Evals');
ylabel('f_{min}');
title('Iteration Statistics for Branin test Function');
---------------------------------------------------------------------------
%}








