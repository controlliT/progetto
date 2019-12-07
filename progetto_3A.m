%%Introduzione

%Progetto di Controlli Automatici - T
%Tipologia III variante A: Controllo di un sistema idroelettrico con condotta forzata
%Specifiche vecchie.

%Gruppo:
%Andrea Proia 0000825784
%Federico Maria Macchiavelli 0000825621
%Mattia Innocenti 0000825046
%Luca Bartolomei 0000825005


%%
%Caratteristiche impianto

%Si tratta di un sistema SISO.

%Tabella delle caratteristiche:
tab.a = pi^2;
tab.a_inv = 1/tab.a;
tab.C_d = 2*pi^2;
tab.R_0 = 45;
tab.eta = 0.6;
tab.W = 30;
tab.omega_n = 1000;
tab.A_n = 0.02;
tab.B_n = 30;
tab.T_a_1 = 0.3;
tab.T_a_0 = 0.033;
tab.x_equilibrio_1 = 10;
tab.x_equilibrio_2 = 6;
tab.x_equilibrio_3 = 6;
tab.u_equilibrio = (tab.x_equilibrio_1/(tab.x_equilibrio_2*abs(tab.x_equilibrio_2))-tab.R_0)/(tab.C_d);

%La costante g non è presente nella tabella:
tab.g_const = 9.80655;

%Metto insieme la condizione di equilibrio dello stato:
% tab_x_equilibrio = [tab_x_equilibrio_1, tab_x_equilibrio_2, tab_x_equilibrio_3];

%Definizione del sistema:
%Il sistema è formato da tre variabili di stato:

%x_1: la pressione dell'acqua sul fondo del bacino.
%x_2: la portata in uscita dal bacino.
%x_3: la portata in ingresso allo stesso.

%L'uscita del sistema viene rappresentata dalla variabile y:
%y: rappresenta l'energia elettrica generata attraverso la turbina.
%y= -eta*x_1*x_2

%Forma di stato:
%x_dot_1 = (g/a)*(-x_2+x_3)
%x_dot_2 = -C_d*u*x_2*|x_2|-R_0*x_2*|x_2|+x_1
%x_dot_3 = 0
 
%%
%Specifiche di progetto

%1 Errore a regime nullo con riferimento a gradino w(t) = W1(t)

%Non posso soddisfare questa specifica: il polo all'origine verrà
%influenzato dallo zero a parte reale positiva, rendendo instabile l'anello
%chiuso.

%Quindi imposto un errore a regime minore di 0.1.
E_regime = 0.1;

%2 Per garantire una certa robustezza del sistema si deve avere un margine di fase M_f* = 45 

Mf=45;

%3 Il sistema può accettare un sovraelongazione percentuale al massimo del 5% : S_% <= 5%

s_perc = 0.05;

%4 Il tempo di assestamento all'1% deve essere tenuto relativamente basso T_a_1 = 0.3 s
%5 Devo abbattere i rumori n di 30 volte che in db =
B_n_db = 20*log10(tab.B_n);

%Specifiche opzionali:

%1_opz Il tempo di assestamento all'1% deve essere tenuto a 0.033 s

%Introduco i vincoli indiretti:

%Dato dal rumore di misura.
%Dato che non ho vincoli di moderazione, questo rimane il vincolo più forte
omega_c_max=tab.omega_n; % 1000 rad/s
%Limite superiore per la frequenza di attraversamento.

%Viene ricavato dalla sovraelongazione con la formula classica.
xi=sqrt(log10(s_perc)^2/(pi^2+log10(s_perc)^2)); %0.6901

%Mf>38.2618° richiesta meno limitatnte delle specifiche (Mf>45°).
%PROCEDIMENTO: calcolo xi con la formula inversa, poi calcolo il Mf e
%valuto la condizione più restrittiva
Mf_s_perc=xi*100; %38.2618

%Dato che le specifiche sono più stringenti uso quelle.
%Mf = Mf_s_perc

%Calcolo la frequenza di attraversamento minimo attraverso la formula:
%460/(Mf* T*) cioè 460/(45 * 0.3)
%Questo limite inferiore è dettato dal tempo di assestamento.
omega_c_min=460/(Mf * tab.T_a_1); % 34.0741 rad/s

%Abbiamo individuato l'intervallo per la pulsazione di attraversamento 
%omega_c_star [34.0741 rad/s, 1000 rad/s]

%%
%Linearizzazione del sistema

%Ridefinisco la prima equazione di stato:
%x_dot_1 = -(g/a) x_2 + (g/a) x_3

%Ridefinisco la seconda equazione di stato:
%x_dot_2 = x_1 -(C_d*u+R_0) x_2*|x_2|
%x_dot_2 = x_1 -(R_0*|x_2|) x_2 - (C_d*x_2*|x_2|) u

%La derivata di x_2 in x_dot_2 non è immediata:
%d/dx_2 (x_1 -(C_d*u+R_0) x_2*|x_2|) = -(2(x_2)^2(C_d*u+R_0))/|x_2|

%Ridefinisco y
%y= -(eta*x_2) x_1
%y= -(eta*x_1) x_2

%Punto di equilibrio dove linearizzare il sistema:
%x_equilibrio = (10,6,6)
%u_equilibrio = (-2.265654245335609)

%Applico la serie di Taylor per linearizzare il sistema:
%A = 0 + d/dx (x_dot)|x=x_equilibrio, u=u_equilibrio
%B = 0 + d/du (x_dot)|x=x_equilibrio, u=u_equilibrio
%C = 0 + d/dx (y)|x=x_equilibrio, u=u_equilibrio
%D = 0 + d/du (y)|x=x_equilibrio, u=u_equilibrio

%Definiamo le matrici A,B,C,D derivabili dalla forma di stato e dall'uscita

%A è una 3x3 perchè devo moltiplicare per le tre equazioni di stato 3x1 e 
%deve saltare fuori un 3x1 quindi 3x3 * 3x1 = 3x1
A = [0, -(tab.g_const/tab.a),                                                                   (tab.g_const/tab.a);
     1, -(2*(tab.x_equilibrio_2^2)*(tab.C_d*tab.u_equilibrio+tab.R_0)/abs(tab.x_equilibrio_2)), 0;
     0, 0,                                                                                      0];

%B è una 3x1 perchè deve moltiplicare per l'ingresso 1x1 e deve saltare
%fuori una 3x1 quindi 3x1 * 1x1 = 3x1
B = [0;
     -(tab.C_d*tab.x_equilibrio_2*abs(tab.x_equilibrio_2));
     0];

%C è 1x3 perchè 1x3 * 3x1 = 1x1
C = [-(tab.eta*tab.x_equilibrio_2), -(tab.eta*tab.x_equilibrio_1), 0];

D = 0;

%Per verificare il contenuto delle variabili globali uso la funzione disp
%disp(A);

%Definisco la funzione di traferimento
s=tf('s');
[Num,Den]=ss2tf(A,B,C,D);
G=tf(Num,Den);
%Stampo G
zpk(G)

%Plot del diagramma di bode

%Margini di visualizzazione del diagramma.
omega_plot_min=10^(-2);
omega_plot_max=10^5;

%Ricavo il bode della G da plottare
[Mag_G,phase_G,omega_G]=bode(G,{omega_plot_min,omega_plot_max});

%Nuova finestra grafica
figure();

%Vincolo sulla omega_c_min
patch([omega_plot_min,omega_c_min,omega_c_min,omega_plot_min],[-200,-200,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Indico la frequenza di attraversamento minima
hold on;
text(omega_plot_min*10,-100, sprintf('w_c^*>=%.2f rad/sec', omega_c_min));

%Vincolo sulla omega_c_max
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[200,200,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Indico la frequenza di attraversamento massima
hold on;
text(omega_c_max*5,60, sprintf('w_c^*<=%.2f rad/sec', omega_c_max));

%Vincolo sull'attenuazione di n
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[-B_n_db,-B_n_db,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%plotto G
hold on;
margin(Mag_G,phase_G,omega_G);

%Vincolo sul margine di fase: -180° + arg(L(jw_c))
hold on;
%Coppie di punti (omega_c_min, -180+Mf), (omega_c_max, -180+Mf), 
%(omega_c_max, -270), (omega_c_min, -270)
patch([omega_c_min,omega_c_max,omega_c_max,omega_c_min],[-180+Mf,-180+Mf,-180,-180],'green','FaceAlpha',0.2,'EdgeAlpha',0); 
grid on;

%%
%Progettazione della rete regolatrice statica

%Guardo il luogo delle radici
figure();
rlocus(G);
 
% Ho un limite molto forte sul guadagno:  circa 0.00039, -68,18 db

%Imposto un regolatore senza polo per evitare instabilità.
R_s = 0.00039;
G_e = R_s*G;

%Stampo gli zeri e i poli di G_e 
zpk(G_e)

%Ricavo i dati sulla G_e
[Mag_G_e,phase_G_e,omega_G_e]=bode(G_e,{omega_plot_min,omega_plot_max});

%Nuova finestra grafica.
figure();

%Vincolo sulla omega_c_min
patch([omega_plot_min,omega_c_min,omega_c_min,omega_plot_min],[-200,-200,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Indico la frequenza di attraversamento minima
hold on;
text(omega_plot_min*10,-100, sprintf('w_c^*>=%.2f rad/sec', omega_c_min));

%Vincolo sulla omega_c_max
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[200,200,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Indico la frequenza di attraversamento massima
hold on;
text(omega_c_max*5,60, sprintf('w_c^*<=%.2f rad/sec', omega_c_max));

%Vincolo sull'attenuazione di n
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[-B_n_db,-B_n_db,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%Plotto G_e
hold on;
margin(Mag_G_e,phase_G_e,omega_G_e);

%Vincolo sul margine di fase: -180° + arg(L(jw_c))
hold on;
%Coppie di punti (w_c_min, -180+Mf), (w_c_max, -180+Mf), (w_c_max, -270),
%(w_c_min, -270)
patch([omega_c_min,omega_c_max,omega_c_max,omega_c_min],[-180+Mf,-180+Mf,-180,-180],'green','FaceAlpha',0.2,'EdgeAlpha',0); 


%Dal grafico si vede che  G_e sta sotto lo 0db, non ho più un
%attraversamento. Non so come progettare la rete dinamica.
%Però posso vedere cosa succede attraverso il luogo delle radici.

figure();
rlocus(G_e);

%Calcolo la retta di attraversamento del tempo di assestamento.
retta_T = 4.6 * tab.T_a_1;

%Il polo vicino all'origine non soddisfa la specifica, dovrei spostarlo a
%sinistra ma c'è l'influenza dello zero reale positivo.


%%

%Calcolo R
R = R_s;

%Dati usati in simulink per la R.
% [NumR, DenR] = tfdata(R);
% NumR = NumR{1,1};
% DenR = DenR{1,1};

%Calcolo L
L = R * G;

%Dati usati in simulink per la L.
[NumL, DenL] = tfdata(L);
NumL = NumL{1,1};
DenL = DenL{1,1};

%Stampo gli zeri e i poli di L

%Ricavo i dati sulla L
[Mag_L,phase_L,omega_L]=bode(L,{omega_plot_min,omega_plot_max});

%Plotto L sul grafico precedente e la confronto con G_e
% figure();
% margin(Mag_L,phase_L,omega_L);
% grid on;

%%
%Chiusura del loop

%Chiudo il loop: F

F=L/(1+L);

%Ricavo informazioni
[Mag_F,phase_F,w_F]=bode(F,{omega_plot_min,omega_plot_max});

%Stampo gli zeri e i poli di F
zpk(F)

%Plotto F
figure();
margin(Mag_F,phase_F,w_F);
grid on;

%Impostazioni per il gradino: imposto un impulso di ampiezza W
stepOption = stepDataOptions('StepAmplitude', tab.W);

%Plotto la risposta a gradino
figure();
step(F, stepOption);

%Per muovere il cursore:
% datacursormode on

%Informazioni sullo step
%Simulo di nuovo lo step ma in questo caso non plotto ma ricavo i dati:
[Y_F,T_F] = step(F, stepOption);
%Imposto un vincolo dell'1% sul tempo di assestamento e ricavo le info:
F_stepinfo = stepinfo(Y_F, T_F,'SettlingTimeThreshold',0.01);
disp(F_stepinfo);

%Stampo il luogo delle radici di F
figure();
rlocus(F);
