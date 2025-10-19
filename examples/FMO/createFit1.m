function [xx,YFIT] = createFit1()
%% Fit: 'FMO fit'.
load fmo.csv;
x=fmo(:,1);
Jor=(x.^2).*fmo(:,2);
J=Jor';
M=500
x1=x(1:M)';
J1=J(1:M)';
Jsm=smooth(x1,J1,0.1,'loess');
Jsm=transpose(Jsm);
J(1:M)=Jsm(1:M);
J=smooth(x,J,0.01,'loess');

[xData, yData] = prepareCurveData( x, J );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.99997;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
%figure( 'Name', 'FMO fit' );
%h = plot( fitresult, xData, yData );
%legend( h, 'J vs. x', 'FMO fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
%xlabel( 'x', 'Interpreter', 'none' );
%ylabel( 'J', 'Interpreter', 'none' );

xx=linspace(x(1),x(length(x)),2000);
YFIT=ppval(fitresult.p,xx);
plot(xx,YFIT);
hold on;
plot(x,Jor);
hold off;

end

