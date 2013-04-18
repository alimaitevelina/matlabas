function f = Easom(x)
f = -cos(x(1))*cos(x(2))*exp(-( (x(1)-pi)^2 + (x(2)-pi)^2 ));

%{ 
---------------------------------------------------------------------------
% RUN symDIRECT   
fun = 'Easom'; x_L = [-100 -100]'; x_U = [100 100]'; 
GLOBAL.tolerance = 0.01; GLOBAL.optimalvalue = -1.0; PriLev = 2;
Result = symDIRECT(fun,x_L,x_U,GLOBAL,PriLev)
---------------------------------------------------------------------------  
% RUN Direct.m (Finkel implementation)
clear all;
% 1. Establish bounds for variables
bounds = [-100 100;-100 100];

% 2. Send options to Direct
% We tell DIRECT that the globalmin = -1.0
% It will stop within 0.01% of solution
options.testflag  = 1; options.globalmin = -1.0; 
options.showits   = 1; options.tol       = 0.01;

% 2a. NEW!
% Pass Function as part of a Matlab Structure
Problem.f = 'Easom';

% 3. Call DIRECT
[fmin,xmin,hist] = Direct(Problem,bounds,options);

% 4. Plot iteration statistics
plot(hist(:,2),hist(:,3))
xlabel('Fcn Evals');
ylabel('f_{min}');
title('Iteration Statistics for Branin test Function');
---------------------------------------------------------------------------
%}






























































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































