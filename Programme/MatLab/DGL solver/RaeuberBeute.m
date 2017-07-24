e1 = 0.1; %Reproduktion der Beute 
e2 = 0.1; %Sterberate der Räuber ohne Beute
y1 = 0.4; %Fressrate der Räuber der Beute pro Lebewesen
y2 = 0.7; %Reproduktion der Räuber pro Beutetier

ode445  = @(t,x) [x(1)*(e1-y1*x(2)); -x(2)*(e2-y2*x(1))] ; 
ode     = @(t,x) [x(1)*(e1-y1*x(2)); -x(2)*(e2-y2*x(1))] ; 
initial = [0.1;0.4];
[y,t]   = ode_solver(ode,initial,0,200,150,6); 
%[t2,y2] = ode45(ode445,[0;50],initial); 
plot(t,y)
legend('Beute','Räuber')