function [fMin2,xMin2]=AdaptyvusMonteCarlo2(funkcija,a1,b1)
% Adaptyvi Monte Carlo (Random Search) realizacija
%1. Sugeneruojame 100 atsitiktiniu tasku int [-10.10]^2
%2. Surandame, kuriame funkcija igyja didziausia (maziausia) reiksme
%a1=-10; %pradine sritis
%b1=10; %pradine sritis
%PALEIDIMAS
%Pvz: a1=-10; b1=10; funkcija=@sincos2;
% [fMin2,xMin2]=AdaptyvusMonteCarlo2(funkcija,a1,b1)
n=2;% Dimensija (matavimas)
k1=50;% tasku (vektoriu) skaicius
x1=a1 + (b1-a1).*rand(k1,n);% Perdaryti kad generuotu dvimacius
f=[];
for i=1:k1
  f(i)=funkcija(x1(i,:));
end
[fMin1,indMin1]=min(f);
%[fMax1,indMax1]=max(f);
xMin1=x1(indMin1,:);
%xMax1=x(indMax1,:);
fprintf('Surastas min=%6.4f taske x=(%6.4f,%6.4f)\n',fMin1,xMin1(1),xMin1(2));
%fprintf('Surastas max=%6.4f taske x=(%6.4f,%6.4f)\n',fMax1,xMax1(1),xMax1(2));
hold on;
scatter(x1(:,1),x1(:,2),'b.');
scatter(xMin1(1),xMin1(2),'r*');
text(xMin1(1)+0.3,xMin1(2),num2str(fMin1));
rectangle('Position',[-10.0,-10.0,20.0,20.0],...
    'LineWidth',5,'LineStyle','--')

%sumazinta sritis

a2=xMin1(1)-2; %sumazinta sritis
b2=xMin1(1)+2; %sumazinta sritis
a3=xMin1(2)-2;
b3=xMin1(2)+2;
n=1;
k2=50;% tasku (vektoriu) skaicius
if (a2<-10) a2=-10; b2=-6;
    if (a3<-10) a3=-10; b3=-6;
        if (b3>10) b3=10; a3=6;
            if (b2>10) b2=10; a2=6;
x2(:,1)=a2 + (b2-a2).*rand(k2,n);% Perdaryti kad generuotu dvimacius
x2(:,2)=a3 + (b3-a3).*rand(k2,n);
f2=[];
for i=1:k2
  f2(i)=funkcija(x2(i,:));
end
[fMin2,indMin2]=min(f2);
%[fMax2,indMax2]=max(f2);
xMin2=x2(indMin2,:);
%xMax2=x2(indMax,:);
fprintf('Surastas min=%6.4f taske x2=(%6.4f,%6.4f)\n',fMin2,xMin2(1),xMin2(2));
%fprintf('Surastas max=%6.4f taske x=(%6.4f,%6.4f)\n',fMax2,xMax2(1),xMax2(2));
hold on;
scatter(x2(:,1),x2(:,2),'b.');
scatter(xMin2(1),xMin2(2),'r*');
text(xMin2(1)+0.3,xMin2(2),num2str(fMin1));
rectangle('Position',[a2,a3,4.0,4.0],...
    'LineWidth',5,'LineStyle','--')
%grafikas2 
%}
