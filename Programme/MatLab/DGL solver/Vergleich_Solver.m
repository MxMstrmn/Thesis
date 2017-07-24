fcn = @(t,x) [x(2);-x(1)] ; 
initial = [0;1]; 
string = {'Einfacher Euler','Verbesserter Euler','Prädiktor-Korrektor','Runge-Kutta 4.Ord.'};
T = linspace(0,2*pi,100);
for i=1:4
    [y,t] = ode_solver(fcn,initial,0,2*pi,10,i);
    subplot(2,2,i); plot(t,y(1,:),T,sin(T)); title(string(i)); xlim([0 2*pi]); ylim([-1.5 1.5]);
end