[X,Y] = meshgrid(-5:.5:10,0:.5:15);
Z = (Y-5*X.^2/(4*pi^2)+5*X./pi-6).^2+10*(1-1/(8*pi))*cos(X)+10;
subplot(1,2,1);
surf(X,Y,Z);
subplot(1,2,2);
contour(X,Y,Z,100);
% Issaugojimas
set(gcf,'PaperPositionMode','auto');
saveas(gcf,'Branin.eps','psc2');