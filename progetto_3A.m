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
tab.x_equilibrio_3 = 7;
tab.u_equilibrio = (tab.x_equilibrio_1/(tab.x_equilibrio_2*abs(tab.x_equilibrio_2))-tab.R_0)/(tab.C_d);

%La costante g non � presente nella tabella:
tab.g_const = 9.80655;

%Metto insieme la condizione di equilibrio dello stato:
% tab_x_equilibrio = [tab_x_equilibrio_1, tab_x_equilibrio_2, tab_x_equilibrio_3];

%Definizione del sistema:
%Il sistema � formato da tre variabili di stato:

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
%2 Per garantire una certa robustezza del sistema si deve avere un margine di fase M_f* = 45 
%3 Il sistema pu� accettare un sovraelongazione percentuale al massimo del 5% : S_% <= 5%

s_perc = 0.05;

%4 Il tempo di assestamento all'1% deve essere tenuto relativamente basso T_a_1 = 0.3 s
%5 Devo abbattere i rumori n di 30 volte che in db =
B_n_db = 20*log10(tab.B_n);

%Specifiche opzionali:

%1_opz Il tempo di assestamento all'1% deve essere tenuto a 0.033 s

%Introduco i vincoli indiretti:

%Dato dal rumore di misura.
%Dato che non ho vincoli di moderazione, questo rimane il vincolo pi� forte
omega_c_max=tab.omega_n; % 1000 rad/s
%Limite superiore per la frequenza di attraversamento.

%Viene ricavato dalla sovraelongazione con la formula classica.
xi=sqrt(log(s_perc)^2/(pi^2+log(s_perc)^2)); %0.6901

%Mf>69.011� richiesta pi� limitatnte delle specifiche (Mf>45�).
%PROCEDIMENTO: calcolo xi con la formula inversa, poi calcolo il Mf e
%valuto la condizione pi� restrittiva
Mf=xi*100; %69.0107 

%Calcolo la frequenza di attraversamento minimo attraverso la formula:
%460/(Mf* T*) cio� 460/(69.011 * 0.3)
%Questo limite inferiore � dettato dal tempo di assestamento.
omega_c_min=460/(Mf * tab.T_a_1); % 22.2188 rad/s

%Abbiamo individuato l'intervallo per la pulsazione di attraversamento 
%w_c* [22.2188 rad/s, 1000 rad/s]

%%
%Linearizzazione del sistema

%Ridefinisco la prima equazione di stato:
%x_dot_1 = -(g/a) x_2 + (g/a) x_3

%Ridefinisco la seconda equazione di stato:
%x_dot_2 = x_1 -(C_d*u+R_0) x_2*|x_2|
%x_dot_2 = x_1 -(R_0*|x_2|) x_2 - (C_d*x_2*|x_2|) u

%La derivata di x_2 in x_dot_2 non � immediata:
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

%A � una 3x3 perch� devo moltiplicare per le tre equazioni di stato 3x1 e 
%deve saltare fuori un 3x1 quindi 3x3 * 3x1 = 3x1
A = [0, -(tab.g_const/tab.a),                                                                   (tab.g_const/tab.a);
     1, -(2*(tab.x_equilibrio_2^2)*(tab.C_d*tab.u_equilibrio+tab.R_0)/abs(tab.x_equilibrio_2)), 0;
     0, 0,                                                                                      0];

%Consideriamo la pressione sul fondo costante.
% A = [0, 0,                                                                                      0;
%      1, -(2*(tab.x_equilibrio_2^2)*(tab.C_d*tab.u_equilibrio+tab.R_0)/abs(tab.x_equilibrio_2)), 0;
%      0, 0,                                                                                      0];

%B � una 3x1 perch� deve moltiplicare per l'ingresso 1x1 e deve saltare
%fuori una 3x1 quindi 3x1 * 1x1 = 3x1
B = [0;
     -(tab.C_d*tab.x_equilibrio_2*abs(tab.x_equilibrio_2));
     0];

%C � 1x3 perch� 1x3 * 3x1 = 1x1
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

%Si evidenzia uno zero che diverge.
%Si evidenzia uno zero nell'origine e un polo nell'origine.
%Si evidenzia due poli non nell'origine.
%Si evidenzia uno zero non nell'origine (quello divergente).

%Plot del diagramma di bode

%Margini di visualizzazione del diagramma.
omega_plot_min=10^(-2);
omega_plot_max=10^5;

%Ricavo il bode della G da plottare
[Mag_G,phase_G,omega_G]=bode(G,{omega_plot_min,omega_plot_max});

%Nuova finestra grafica
figure();

%Vincolo sulla w_c_min
patch([omega_plot_min,omega_c_min,omega_c_min,omega_plot_min],[-200,-200,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0);
text(omega_c_min-22,-100,'w_c^* >= 22.2188 rad/sec');

%Vincolo sulla w_c_max
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[200,200,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Vincolo sull'attenuazione di n
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[-B_n_db,-B_n_db,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%plotto G
hold on;
margin(Mag_G,phase_G,omega_G);

%Vincolo sul margine di fase: -180� + arg(L(jw_c))
hold on;
%Coppie di punti (w_c_min, -180+Mf), (w_c_max, -180+Mf), (w_c_max, -270),
%(w_c_min, -270)
patch([omega_c_min,omega_c_max,omega_c_max,omega_c_min],[-180+Mf,-180+Mf,-270,-270],'green','FaceAlpha',0.2,'EdgeAlpha',0); 
grid on;



%-------------------------------------------------------------------------------------------------------------------------------------


%%
%Progettazione della rete regolatrice statica

%Ho bisogno di un polo per il vincolo numero 1.
%Volgio un e_inf = 0
%Il guadagno statico resta libero: verr� modificato se necessario.
R_s = 1/s;
G_e = R_s*G;

%Stampo gli zeri e i poli di G_e 
zpk(G_e)

%Ricavo i dati sulla G_e
[Mag_G_e,phase_G_e,omega_G_e]=bode(G_e,{omega_plot_min,omega_plot_max});

%Nuova finestra grafica.
figure();

%Vincolo sulla w_c_min
patch([omega_plot_min,omega_c_min,omega_c_min,omega_plot_min],[-200,-200,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Vincolo sulla w_c_max
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[120,120,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Vincolo sull'attenuazione di n
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[-B_n_db,-B_n_db,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 


%Plotto G_e
hold on;
bodeplot(G_e, {omega_plot_min,omega_plot_max});
% margin(Mag_G_e,phase_G_e,omega_G_e);

%Vincolo sul margine di fase: -180� + arg(L(jw_c))
hold on;
%Coppie di punti (w_c_min, -180+Mf), (w_c_max, -180+Mf), (w_c_max, -270),
%(w_c_min, -270)
patch([omega_c_min,omega_c_max,omega_c_max,omega_c_min],[-180+Mf,-180+Mf,-180,-180],'green','FaceAlpha',0.2,'EdgeAlpha',0); 

%Dal grafico si deduce che siamo caduti in uno scenario B

%%
%Progettazione della rete regolatrice dinamica

%Si tratta della frequenza di attraversamento scelta.
omega_c_star_ant = 150;

%Alzo il margine di fase per contrastare i due poli della rete ritardatrice
%da progettare successivamente per rispettare il vincolo su n.
Mf_ant=Mf+10;

%Ricavo i dati di attraversamento di G_e
[Mag_G_e_omega_c_star_ant,phase_G_e_omega_c_star_ant,omega_c_star_ant]=bode(G_e, omega_c_star_ant);

Mag_G_e_omega_c_star_ant_db = 20*log(Mag_G_e_omega_c_star_ant);

M_star_ant = 10^-(Mag_G_e_omega_c_star_ant_db/20);
phi_star_ant = Mf_ant-180-phase_G_e_omega_c_star_ant;

tau_rete_ant = (M_star_ant-cos(phi_star_ant*pi/180))/(omega_c_star_ant*sin(phi_star_ant*pi/180));
tau_alpha_rete_ant = (cos(phi_star_ant*pi/180)-1/M_star_ant)/(omega_c_star_ant*sin(phi_star_ant*pi/180));

alpha_rete_ant = tau_alpha_rete_ant/tau_rete_ant;

%Ricavo il regolatore dinamico (anticipatore).
R_d_ant = (1+s*tau_rete_ant)/(1+tau_alpha_rete_ant*s);

%Stampo il grafico per vedere in che scenario siamo caduti.
G_e_1 = R_d_ant  * G_e;

%Stampo gli zeri e i poli di G_e_1
% zpk(G_e_1)

%Ricavo i dati sulla G_e_1
[Mag_G_e_1,phase_G_e_1,omega_G_e_1]=bode(G_e_1,{omega_plot_min,omega_plot_max});

figure();
%Vincolo sulla w_c_min
patch([omega_plot_min,omega_c_min,omega_c_min,omega_plot_min],[-200,-200,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Vincolo sulla w_c_max
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[120,120,0,0],'yellow','FaceAlpha',0.3,'EdgeAlpha',0); 

%Vincolo sull'attenuazione di n
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[-B_n_db,-B_n_db,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%Plotto G_e_1
hold on;
% bodeplot(G_e_1, {omega_plot_min,omega_plot_max});
margin(Mag_G_e_1,phase_G_e_1,omega_G_e_1);

%Vincolo sul margine di fase: -180� + arg(L(jw_c))
hold on;
%Coppie di punti (w_c_min, -180+Mf), (w_c_max, -180+Mf), (w_c_max, -270),
%(w_c_min, -270)
patch([omega_c_min,omega_c_max,omega_c_max,omega_c_min],[-180+Mf,-180+Mf,-180,-180],'green','FaceAlpha',0.2,'EdgeAlpha',0);

%Siamo caduti in uno scenario A:
%Dato che il guadagno � libero lo uso per correggere la rete.
%Abbiamo calcolato il valore di mu attraverso il grafico di bode.

omega_c_star_mu = 45;

[Mag_G_e_1_omega_c_star_ant,phase_G_e_1_omega_c_star_mu,omega_c_star_mu]=bode(G_e_1, omega_c_star_mu);

mu_d = 1/Mag_G_e_1_omega_c_star_ant;

G_e_2 = mu_d * G_e_1;

%Plotto G_e_2
[Mag_G_e_2,phase_G_e_2,omega_G_e_2]=bode(G_e_2,{omega_plot_min,omega_plot_max});
hold on;
% bodeplot(G_e_2, {omega_plot_min,omega_plot_max});
margin(Mag_G_e_2,phase_G_e_2,omega_G_e_2);

%La specifica di attraversamento e margine di fase � rispettata, manca
%solo il vincolo sulla n: uso una rete ritardatrice.
%%
%RETE RITARDATRICE PER RISOLUZIONE DI n(t)
%Si tratta della frequenza di inizio del polo.
omega_c_star_rit = 150;
alpha_rete_rit = 0.01;

tau_rete_ritardatrice = 1/omega_c_star_rit;
tau_alpha_rete_ritardatrice = alpha_rete_rit * tau_rete_ritardatrice;

R_d_rit = (1+s*tau_alpha_rete_ritardatrice)^2/(1+tau_rete_ritardatrice*s)^2;


%%
%Calcolo L con aggiunta della rete ritardatrice
L = mu_d * R_d_ant * G_e;

[NumL, DenL] = tfdata(L);
NumL = NumL{1,1};
DenL = DenL{1,1};

%Stampo gli zeri e i poli di L

%Ricavo i dati sulla L
[Mag_L,phase_L,omega_L]=bode(L,{omega_plot_min,omega_plot_max});

%Plotto L sul grafico precedente e la confronto con G_e
hold on;
margin(Mag_L,phase_L,omega_L);
grid on;

%Il grafico sembra reggere.

%%
%Chiusura del loop

%Chiudo il loop: F

F=L/(1+L);

%Ricavo informazioni
[Mag_F,phase_F,w_F]=bode(F,{omega_plot_min,omega_plot_max});

%Stampo gli zeri e i poli di F
zpk(F)

%Plotto F
% figure();
% margin(Mag_F,phase_F,w_F);
% grid on;

%Informazioni sullo step
stepinfo(F)

figure();
step(F);


