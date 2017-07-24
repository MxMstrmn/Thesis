a = 0.246;
tAB = 3;
res = 150;

%==========================================================% 
% Plot über ausgeählten Weg in der 1.BZ                    %
%==========================================================%

G = [0;0];
K = pi/a*[2/sqrt(3) ; 2/3];
M = pi/a*[1/sqrt(3) ; 1];

PATH = [G M K G];
n    = [res 0.6*res res];
ki   = vecspace(PATH,n);

f =@(ki,a) exp(-1i*ki(1,:)*a/sqrt(3)) + exp(1i*ki(1,:)*a/(2*sqrt(3))).*(exp(-1i*ki(2,:)*a/2)+exp(1i*ki(2,:)*a/2));
f_ki = f(ki,a);

E = zeros(size(ki));
E_as = E;
for j=1:5
for i=1:length(f_ki)
    H = [0 -tAB*f_ki(i);
        -tAB*conj(f_ki(i))  0];
    % Symmetriebrechung 
    delta = 0.1*j^2.5;
    H_as = H + delta*[1 0; 0 -1];
    
    E_as(:,i) = real(eig(H_as));
    E(:,i) = real(eig(H));
end


x = linspace(1,length(E(1,:)),length(E(1,:)));
subplot(2,1,1);
plot(x,E,'b-')
x_tick = [1, n(1), sum(n(1:2))-1, sum(n)-2];
set(gca, 'xtick', x_tick)
set(gca, 'xticklabel', {'\Gamma', 'M', 'K', '\Gamma'})

subplot(2,1,2);
plot(x,E_as,'b-')
x_tick = [1, n(1), sum(n(1:2))-1, sum(n)-2];
set(gca, 'xtick', x_tick)
set(gca, 'xticklabel', {'\Gamma', 'M', 'K', '\Gamma'})
hold on
end



%==========================================================% 
% Plot über Gebiet                                         %
%==========================================================%


[kx,ky] = meshgrid(-1.2*norm(K):norm(K)/100:1.2*norm(K));
f_ki = f([kx(:)';ky(:)'],a);

for i=1:length(f_ki)
    H = [0 -tAB*f_ki(i);
        -tAB*conj(f_ki(i))  0];
    E(:,i) = real(eig(H));
end

figure
E1 = reshape(E(1,:),size(kx));
s = surf(kx,ky,E1);
set(s,'edgecolor','none'); 

