function f = Branin(x)
f = (x(2)-5*x(1)^2/(4*pi^2)+5*x(1)/pi-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
%{ 
---------------------------------------------------------------------------
% RUN gblSolve (gblSolve.m by Bjorkman)   
fun = 'Branin'; x_L = [-5 0]'; x_U = [10 15]'; GLOBAL.tolerance = 0.01; 
GLOBAL.optimalvalue = 0.397887357729739; PriLev = 2;
Result = gblSolve(fun,x_L,x_U,GLOBAL,PriLev)
---------------------------------------------------------------------------
% RUN gblSolveMod (gblSolveMod.m by Bjorkman. Slighly modified)   
fun = 'Branin'; x_L = [-5 0]'; x_U = [10 15]'; GLOBAL.tolerance = 0.01; 
GLOBAL.optimalvalue = 0.397887357729739; PriLev = 2;
Result = gblSolveMod(fun,x_L,x_U,GLOBAL,PriLev)

% SCATTER PLOT WITH CONTOUR LINE 
C=Result.GLOBAL.C;
plot(C(1,:),C(2,:),'.');

[X,Y] = meshgrid(-5:.1:10,0:.1:15);
Z = (Y-5*X.^2/(4*pi^2)+5*X./pi-6).^2+10*(1-1/(8*pi))*cos(X)+10;
contour(X,Y,Z,10);

V=Result.GLOBAL.V;
plot(V(1,:),V(2,:),'.');
---------------------------------------------------------------------------  
% RUN Direct.m (Finkel implementation)
clear all;
% 1. Establish bounds for variables
bounds = [-5 10;0 15];

% 2. Send options to Direct
% We tell DIRECT that the globalmin = 0.397887357729739
% It will stop within 0.01% of solution
options.testflag  = 1; options.globalmin = 0.397887357729739; 
options.showits   = 1; options.tol       = 0.01;

% 2a. NEW!
% Pass Function as part of a Matlab Structure
Problem.f = 'Branin';

% 3. Call DIRECT
[fmin,xmin,hist] = Direct(Problem,bounds,options);

% 4. Plot iteration statistics
plot(hist(:,2),hist(:,3))
xlabel('Fcn Evals');
ylabel('f_{min}');
title('Iteration Statistics for Branin test Function');
---------------------------------------------------------------------------
%}





