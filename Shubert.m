function f = Shubert(x)
s1 = 0; 
s2 = 0;
for i = 1:5;   
    s1 = s1+i*cos((i+1)*x(1)+i);
    s2 = s2+i*cos((i+1)*x(2)+i);
end
f = s1*s2;

%{
---------------------------------------------------------------------------
% RUN symDIRECT (gblSolve.m by Bjorkman. Slighly modified)  
fun = 'Shubert'; x_L = [-10 -10]'; x_U = [10 10]'; 
GLOBAL.tolerance = 0.01;GLOBAL.optimalvalue = -186.730908831024;PriLev = 2;
Result = symDIRECT(fun,x_L,x_U,GLOBAL,PriLev)
---------------------------------------------------------------------------  
% RUN Direct.m (Finkel implementation)
clear all;
% 1. Establish bounds for variables
bounds = [-10 10;-10 10];

% 2. Send options to Direct
% We tell DIRECT that the globalmin = -186.730908831024
% It will stop within 0.01% of solution
options.testflag  = 1; options.globalmin = -186.730908831024; 
options.showits   = 1; options.tol       = 0.01;

% 2a. NEW!
% Pass Function as part of a Matlab Structure
Problem.f = 'Shubert';

% 3. Call DIRECT
[fmin,xmin,hist] = Direct(Problem,bounds,options);

% 4. Plot iteration statistics
plot(hist(:,2),hist(:,3))
xlabel('Fcn Evals');
ylabel('f_{min}');
title('Iteration Statistics for Shubert test Function');
---------------------------------------------------------------------------
%}