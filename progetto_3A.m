%%Introduzione

%Progetto di Controlli Automatici - T
%Tipologia III variante A: Controllo di un sistema idroelettrico con condotta forzata

%Gruppo:
%Andrea Proia 0000825784
%Federico Maria Macchiavelli 0000825621
%Mattia Innocenti 0000825046
%Luca Bartolomei 0000825005


%%
%Caratteristiche impianto

%Si tratta di un sistema SISO.

%Tabella delle caratteristiche:
tab_a = pi^2;
tab_a_inv = 1/(pi^2);
tab_C_d = 2*pi^2;
tab_R_0 = 45;
tab_eta = 0.6;
tab_W = 30;
tab_omega_n = 1000;
tab_A_n = 0.02;
tab_B_n = 30;
tab_h_perc = 1;
tab_T_a_h_perc = 0.3;
tab_T_a_0 = 0.033;
tab_x_equilibrio_1 = 10;
tab_x_equilibrio_2 = 6;
tab_x_equilibrio_3 = 6;
tab_u_equilibrio = (tab_x_equilibrio_1/(tab_x_equilibrio_2*abs(tab_x_equilibrio_2))-tab_R_0)/(tab_C_d);

%La costante g non è presente nella tabella:
tab_g_const = 9.80655;

%Metto insieme la condizione di equilibrio dello stato:
tab_x_equilibrio = [tab_x_equilibrio_1, tab_x_equilibrio_2, tab_x_equilibrio_3];

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
%2 Per garantire una certa robustezza del sistema si deve avere un margine di fase M_f* = 45°
%3 Il sistema può accettare un sovraelongazione percentuale al massimo del 5% : S_% <= 5%

s_perc = 0.05;

%4 Il tempo di assestamento all'1% deve essere tenuto relativamente basso T_a_1 = 0.3 s
%5 Devo abbattere i rumori n di 30 volte che in db =
B_n_db = 20*log(tab_B_n);

%Specifiche opzionali:

%1_opz Il tempo di assestamento all'1% deve essere tenuto a 0.033 s

%Introduco i vincoli indiretti:

%Dato dal rumore di misura.
w_c_max=tab_omega_n; % 1000 rad/s

%Viene ricavato dalla sovraelongazione con la formula classica.
xi=sqrt(log(s_perc)^2/(pi^2+log(s_perc)^2));

%Mf>69.011° richiesta più limitatnte delle specifiche (Mf>45°),
Mf=xi*100;

%Calcolo la frequenza di attraversamento minimo attraverso la formula:
%460/(Mf* T*) cioè 460/(69.011 * 0.3)
w_c_min=460/(Mf * tab_T_a_h_perc); % 22.2188 rad/s

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
A = [0, -(tab_g_const/tab_a),                                                                   (tab_g_const/tab_a);
     1, -(2*(tab_x_equilibrio_2^2)*(tab_C_d*tab_u_equilibrio+tab_R_0)/abs(tab_x_equilibrio_2)), 0;
     0, 0,                                                                                      0];

%B è una 3x1 perchè deve moltiplicare per l'ingresso 1x1 e deve saltare
%fuori una 3x1 quindi 3x1 * 1x1 = 3x1
B = [0;
     -(tab_C_d*tab_x_equilibrio_2*abs(tab_x_equilibrio_2));
     0];

%C è 1x3 perchè 1x3 * 3x1 = 1x1
C = [-(tab_eta*tab_x_equilibrio_2), -(tab_eta*tab_x_equilibrio_1),0];

D = 0;

%Definisco la funzione di traferimento
s=tf('s');
[Num,Den]=ss2tf(A,B,C,D);
G=tf(Num,Den);
%Stampo G
zpk(G)

%Si evidenzia uno zero che diverge.
%Si evidenzia uno zero nell'origine e un polo nell'origine.
%Si evidenzia due poli non nell'origine.

%Plot del diagramma di bode

%Margini di visualizzazione del diagramma.
w_plot_min=10^(-2);
w_plot_max=10^5;

%Ricavo il bode della G da plottare
[Mag,phase,w]=bode(G,{w_plot_min,w_plot_max});


%Nuova finestra grafica
% figure();
% 
% %Vincolo sulla w_c_min
% patch([w_plot_min,w_c_min,w_c_min,w_plot_min],[-200,-200,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 
% 
% %Vincolo sulla w_c_max
% hold on;
% patch([w_plot_max,w_c_max,w_c_max,w_plot_max],[120,120,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 
% 
% %Vincolo sull'attenuazione di n
% hold on;
% patch([w_plot_max,w_c_max,w_c_max,w_plot_max],[-B_n_db,-B_n_db,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 
% 
% %plotto G
% hold on;
% margin(Mag,phase,w);
% 
% %Vincolo sul margine di fase: -180° + arg(L(jw_c))
% hold on;
% %Coppie di punti (w_c_min, -180+Mf), (w_c_max, -180+Mf), (w_c_max, -270),
% %(w_c_min, -270)
% patch([w_c_min,w_c_max,w_c_max,w_c_min],[-180+Mf,-180+Mf,-270,-270],'green','FaceAlpha',0.2,'EdgeAlpha',0); 
% grid on;

%%
%Progettazione della rete regolatrice statica

%Ho bisogno di un polo per il vincolo numero 1.
%Il guadagno statico resta libero: verrà modificato se necessario.
R_s = 1/s;
G_e = R_s*G;

%Stampo gli zeri e i poli di G_e 
zpk(G_e)

%Ricavo i dati sulla G_e
[Mag_e,phase_e,w_e]=bode(G_e,{w_plot_min,w_plot_max});

%Nuova finestra grafica.
figure();

%Vincolo sulla w_c_min
patch([w_plot_min,w_c_min,w_c_min,w_plot_min],[-200,-200,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Vincolo sulla w_c_max
hold on;
patch([w_plot_max,w_c_max,w_c_max,w_plot_max],[120,120,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Vincolo sull'attenuazione di n
hold on;
patch([w_plot_max,w_c_max,w_c_max,w_plot_max],[-B_n_db,-B_n_db,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%Plotto G_e
hold on;
bodeplot(G_e, {w_plot_min,w_plot_max});
% margin(Mag_e,phase_e,w_e);

%Vincolo sul margine di fase: -180° + arg(L(jw_c))
hold on;
%Coppie di punti (w_c_min, -180+Mf), (w_c_max, -180+Mf), (w_c_max, -270),
%(w_c_min, -270)
patch([w_c_min,w_c_max,w_c_max,w_c_min],[-180+Mf,-180+Mf,-270,-270],'green','FaceAlpha',0.2,'EdgeAlpha',0); 

%%
%Progettazione della rete regolatrice dinamica

%Si tratta della frequenza di attraversamento scelta. NON FUNZIONA BOH
omega_c_star = 45;

%Ricavo i dati di attraversamento di G_e
[Mag_G_e_omega_c_star,phase_G_e_omega_c_star,omega_c_star]=bode(G_e, omega_c_star);

Mag_G_e_omega_c_star_db = 20*log(Mag_G_e_omega_c_star);

M_star = 1/10^(Mag_G_e_omega_c_star_db/20);

phi_star = Mf-180-phase_G_e_omega_c_star;

tau_rete_anticipatrice = (M_star-cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
tau_alpha_rete_anticipatrice = (cos(phi_star*pi/180)-1/M_star)/(omega_c_star*sin(phi_star*pi/180));

alpha_rete_anticipatrice = tau_alpha_rete_anticipatrice/tau_rete_anticipatrice;

%Ricavo il regolatore dinamico.
R_d = (1+s*tau_rete_anticipatrice)/(1+tau_alpha_rete_anticipatrice*s);

%Creo la funzione ad anello aperto (L).
L = R_d*G_e;

%Stampo gli zeri e i poli di L
zpk(L)
%Ricavo i dati sulla L
[Mag_L,phase_L,w_L]=bode(L,{w_plot_min,w_plot_max});

%Plotto L sul grafico precedente e la confronto con G_e
hold on;
margin(Mag_L,phase_L,w_L);
grid on;

%Il grafico sembra reggere.

%%
%Regolatore finale
R=R_s*R_d;

%%
%Chiusura del loop

%Chiudo il loop: F

F=L/(1+L);

%Ricavo informazioni
[Mag_F,phase_F,w_F]=bode(F,{w_plot_min,w_plot_max});

%Stampo gli zeri e i poli di F
zpk(F)

%Plotto F
figure();
margin(Mag_F,phase_F,w_F);
grid on;

figure();
step(F);

