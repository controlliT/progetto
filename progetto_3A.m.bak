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
xi=sqrt(log(s_perc)^2/(pi^2+log(s_perc)^2)); %0.6901

%Mf>69.01 gradi richiesta meno limitatnte delle specifiche (Mf>45 gradi).
%PROCEDIMENTO: calcolo xi con la formula inversa, poi calcolo il Mf e
%valuto la condizione più restrittiva
Mf_s_perc=xi*100; %69.01

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

%Noto uno zero a parte reale positiva

%Imposto un regolatore con un polo nell'origine per rispettare il vincolo
%e_inf=0
R_s = 1/s;
G_e = R_s*G;

%Stampo gli zeri e i poli di G_e 
zpk(G_e)

%Guardo il luogo delle radici
figure(1);
rlocus(G_e);
title("Luogo delle radici, regolatore Rs=1/s)");

%Calcolo la retta di attraversamento del tempo di assestamento.
retta_T = 4.6 * tab.T_a_1;

hold on;
plot([-retta_T, -retta_T], [-10, 10]);

grid on;

%Commento finale:
%Per rispettare il vincolo di tempo di assestamento dovrei impostare un
%guadagno dinamico di almeno 0.000289 ma così facendo sposto il polo
%nell'origine nella parte a destra dell'asse immaginario, rendendo
%instabile il sistema. Dovrei violare il vincolo e_inf = 0 per avere un
%margine molto basso sul guadagno: ho comunque un polo vicino allo zero.
%Quindi non è detto che riesca a regolare comunque la rete con il livello
%di conoscenza attuale.

%Se imposto un regolatore statico senza polo:

R_s = 1;
G_e = R_s*G;

%Stampo gli zeri e i poli di G_e 
zpk(G_e)

%Guardo il luogo delle radici
figure(2);
rlocus(G_e);
title("Luogo delle radici, regolatore Rs=1)");

%Calcolo la retta di attraversamento del tempo di assestamento.
retta_T = 4.6 * tab.T_a_1;

hold on;
plot([-retta_T, -retta_T], [-10, 10]);

grid on;

%Ho bisogno di una rete anticipatrice per cancellare il polo in -0.33

alpha_rete_ant = 0.001;
tau_rete_ant = 1/0.33;
tau_alpha_rete_ant = tau_rete_ant*alpha_rete_ant;

%Ricavo il regolatore dinamico (anticipatore).
R_d_ant = (1+s*tau_rete_ant)/(1+tau_alpha_rete_ant*s);

G_e_1 = R_d_ant * G_e;


%Guardo il luogo delle radici
figure(3);
rlocus(G_e_1);
title("Luogo delle radici G_e = Rd * Rs * G");

%Calcolo la retta di attraversamento del tempo di assestamento.
retta_T = 4.6 * tab.T_a_1;

hold on;
plot([-retta_T, -retta_T], [-10, 10]);

grid on;

%Devo abbassare di parecchio il guadagno.
mu_d = 0.000055;

L = mu_d * G_e_1;

%Vendo come si comporta F
F=L/(1+L);

zpk(F)

figure(4);
pzmap(F);
title("Zeri e poli di F");
grid on;

%Margini di visualizzazione del diagramma.
omega_plot_min=10^(-2);
omega_plot_max=10^5;

%Ricavo il bode della G da plottare
[Mag_F,phase_F,omega_F]=bode(F,{omega_plot_min,omega_plot_max});
[Mag_L,phase_L,omega_L]=bode(L,{omega_plot_min,omega_plot_max});

figure(5);


% %Vincolo sulla omega_c_min
% patch([omega_plot_min,omega_c_min,omega_c_min,omega_plot_min],[-200,-200,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 
% 
% %Indico la frequenza di attraversamento minima
% hold on;
% text(omega_plot_min*10,-100, sprintf('w_c^*>=%.2f rad/sec', omega_c_min));

%Vincolo sulla omega_c_max
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[200,200,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%Indico la frequenza di attraversamento massima
hold on;
text(omega_c_max*5,60, sprintf('w_c^*<=%.2f rad/sec', omega_c_max));

%Vincolo sull'attenuazione di n
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[-B_n_db,-B_n_db,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%plotto F
hold on;
margin(Mag_F,phase_F,omega_F);

%plotto L
hold on;
margin(Mag_L,phase_L,omega_L);

title("Bode di L e F con vincoli");

% %Vincolo sul margine di fase: -180gradi + arg(L(jw_c))
% hold on;
% %Coppie di punti (omega_c_min, -180+Mf), (omega_c_max, -180+Mf), 
% %(omega_c_max, -270), (omega_c_min, -270)
% patch([omega_c_min,omega_c_max,omega_c_max,omega_c_min],[-180+Mf,-180+Mf,-180,-180],'red','FaceAlpha',0.2,'EdgeAlpha',0); 
% grid on;

%Impostazioni per il gradino: imposto un impulso di ampiezza W
stepOption = stepDataOptions('StepAmplitude', tab.W);
figure(6);
step(F, stepDataOptions);
title("Risposta al gradino (W) di F");
grid on;
