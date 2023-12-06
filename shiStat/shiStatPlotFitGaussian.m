function [PSE,Spread,Func,PLOT] = shiStatPlotFitGaussian(X,Y,Color)

% fits the Y-X scatter plot with a gaussian function (normcdf)
%

if nargin < 3
    Color = 'b';
end

f = @(beta,x)((1+erf((x-beta(1))/(sqrt(2)*beta(2))))/2);

b0 = [0.5 0.8];
b = nlinfit(X,Y,f,b0);
PSE = b(1);
Spread = b(2);

PLOT = plot(100*(0:0.02:1),100*f(b,0:0.02:1),'Color',Color,'LineWidth',2);
hold on;
plot(100*X,100*Y,'o','Color',Color,'LineWidth',2,'MarkerFaceColor',Color);
plot(100*[0,PSE],100*[0.5,0.5],'--','Color','k','LineWidth',1.5);
plot(100*[PSE,PSE],100*[0,0.5],'--','Color',Color,'LineWidth',1.5);
hold off;

ax = gca;
ax.XAxis.TickLabelFormat = '%g%%';
ax.YAxis.TickLabelFormat = '%g%%';

Func = @(x)((1+erf((x-b(1))/(sqrt(2)*b(2))))/2);

