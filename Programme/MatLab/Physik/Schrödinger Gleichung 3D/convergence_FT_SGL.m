%================================================% 
% Test convergence with weighting function den
%================================================%
% Results shwon in calc file 08/05/17

y           = [];
dist        = [];
int_pnts    = 40;
temp_state  = 0; 
Int         = [];
steps       = 20;
N_min       = 3; 

for i=1:5
% A N M E R K U N G
%====================================
% Gewichtung mit der Funktion, um selbst die Prioritäten zu setzen, aber
% hier gehen vorteile von Gauß/ GTS u.U. verloren! 

% c       = get_coefficient(10,5,100,1);
% density = @(x)x.^c(1).*exp(c(2)-c(3).*x);    
% 
% I = [0 20 100 500]; 
% for j=1:length(I)-1; Int(j) = I(j) ; end %+ (I(j+1)-I(j))/2
% N = density(Int);
% N = round(int_pnts/sum(N)*N,0);
% N(1) = N(2);
% N(N<N_min) = N_min;

%====================================

I = [0 2.5 5 10 25 100 500];
N = [50 50 50 50 5 5];

[EW, states, k]    = SGL_FT_Hebung(I,N);
bound_states       = EW(EW<0);
states             = states(:,1:length(bound_states));

y(i)        = bound_states(1);
dist(i)     = abs(bound_states(1)-temp_state);
temp_state  = bound_states(1);
end

%%
%===============================================================%
%                    P L O T T I N G                            %         
%===============================================================%

% Diese Plots funktionieren für die Konvergenzuntersuchung

% subplot(3,1,1)
% bar(Int,N)
% legend(sprintf('Gewichtete Verteilung für I=[0,%d], pro step: N+=%d', max(I),int_pnts)) %
% xlabel('k-Intervalle')
% ylabel('Stützstellen')
% 
% subplot(3,1,2)
% plot(1:i-1,dist(2:end),'b-x')
% legend(sprintf('I=[0:5:5*i]')); %[0:10:100 120:20:300 350:50:500 700:200:1500]
% xlabel('Laufindex i')
% ylabel('Änderung zum vorigen Eigenwert')
% 
% subplot(3,1,3)
% plot(1:i,abs(y/13600),'b-x')
% legend(sprintf('Energieeigenwert'))
% xlabel('Laufindex i')
% ylabel('Relative Genauigkeit')


f = figure;
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = ('Wasserstoff: Wellenfunktionen und Spektrum'); 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';


[k,idx] = sort(k);
states = states(idx,:);
rows = 4;
for i=1:rows-1
subplot(rows,3,[3*i-2 3*i-1],'Parent', p)
% numerical results 
PSI = (states(:,i));
numerical = plot(k,PSI,'r.-');hold on;
plot(k,0*k,'k-')
% analytical results
if i < 4
    [ka,PSIa(i,:)] = FT_schlecht([0,40],50,i);
    analytical = plot(ka,PSIa(i,:),'b.');
    legend([numerical(1) analytical],sprintf('Eigenfunktion zum EW=%4.1f meV', bound_states(i)),'Analytisches Ergebnis')
end
xlabel('k')
xlim([0,40-5*i]);
ylabel('')
end

subplot(rows,3,[10 11])
for i=1:3
    y  = states(:,i);
    y2 = PSIa(i,:);
    numerical2 = plot(k,y.^2,'.');hold on
    %analytical2 = plot(ka,y2.^2,'bo');
    xlim([0,7]);
end

subplot(rows,3,[3 6 9 12],'Parent', p)
analytisch = -13.6*1./[1:length(bound_states)].^2;
a = plot([0 1],[bound_states bound_states]*1e-3,'r-'); hold on
b = plot([1 2],[analytisch' analytisch'],'b-');
legend([a(1) b(1)], 'numerisch','analytisch','Location','best')
xlabel('Rydbergspektrum')
ylabel('Energiewerte in eV')

