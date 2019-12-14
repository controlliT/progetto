%%
%--Introduzione--

%Progetto di Controlli Automatici - T
%Tipologia III variante A: Controllo di un sistema idroelettrico con condotta forzata
%Nuova versione con le specifiche nuove.

%Gruppo:
%Andrea Proia 0000825784
%Federico Maria Macchiavelli 0000825621
%Mattia Innocenti 0000825046
%Luca Bartolomei 0000825005


%%
%--Caratteristiche impianto--

%Si tratta di un sistema SISO.

%Tabella delle caratteristiche:
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
tab.u_equilibrio = (tab.x_equilibrio_1/(tab.x_equilibrio_2*abs(tab.x_equilibrio_2))-tab.R_0)/(tab.C_d);

%Definizione del sistema:
%Il sistema è formato da due variabili di stato:

%x_1: la pressione dell'acqua sul fondo del bacino. Consideriamo la 
%pressione sul fondo costante.
%x_2: la portata in uscita dal bacino.

%L'uscita del sistema viene rappresentata dalla variabile y:
%y: rappresenta l'energia elettrica generata attraverso la turbina.
%y= -eta*x_1*x_2

%Forma di stato:
%x_dot_1 = 0
%x_dot_2 = -C_d*u*x_2*|x_2|-R_0*x_2*|x_2|+x_1
 
%%
%--Specifiche di progetto--

%1 Errore a regime nullo con riferimento a gradino w(t) = W1(t)
%2 Per garantire una certa robustezza del sistema si deve avere un margine di fase M_f* = 45 

% Mf=45;

%3 Il sistema può accettare un sovraelongazione percentuale al massimo del 5% : S_% <= 5%

s_perc = 0.05;

%4 Il tempo di assestamento all'1% deve essere tenuto relativamente basso T_a_1 = 0.3 s
%5 Devo abbattere i rumori n di 30 volte che in db =
B_n_db = 20*log10(tab.B_n);

%Specifiche opzionali:

%1 Il tempo di assestamento all'1% deve essere tenuto a 0.033 s

%Introduco i vincoli indiretti:

%Dato dal rumore di misura.
%Dato che non ho vincoli di moderazione, questo rimane il vincolo più forte
omega_c_max=tab.omega_n; % 1000 rad/s
%Limite superiore per la frequenza di attraversamento.

%Viene ricavato dalla sovraelongazione con la formula classica.
%Ocio: qui va lasciato il log normale
xi=sqrt(log(s_perc)^2/(pi^2+log(s_perc)^2)); %0.6901

%PROCEDIMENTO: calcolo xi con la formula inversa, poi calcolo il Mf e
%valuto la condizione più restrittiva
Mf_s_perc=xi*100; %69.01 

%Mf>69.01° richiesta più limitatnte delle specifiche (Mf>45°).
Mf = Mf_s_perc;

%Calcolo la frequenza di attraversamento minima attraverso la formula:
%460/(Mf* T*) cioè 460/(45 * 0.3)
%Questo limite inferiore è dettato dal tempo di assestamento.

%ATTENZIONE: può subire variazioni in quanto va tenuto in conto Mf del
%sistema in retroazione (F).
omega_c_min = 460/(Mf * tab.T_a_1); 

%Calcolo la frequenza di attraversamento minima opzionale:
%ATTENZIONE: può subire variazioni in quanto va tenuto in conto Mf del
%sistema in retroazione (F).
omega_c_min_fac = 460/(Mf * tab.T_a_0); 

%Calcolo anche l'asse di specifica del tempo di assestamento all'1%
T_axis = 4.6 / tab.T_a_1;
T_axis_fac = 4.6 / tab.T_a_0;

%Abbiamo individuato l'intervallo per la pulsazione di attraversamento 
fprintf('Range di attraversamento obbligatorio [%.2f rad/s, %.2f rad/s]\n', omega_c_min, omega_c_max);
%Versione opzionale:
fprintf('Range di attraversamento opzionale [%.2f rad/s, %.2f rad/s]\n', omega_c_min_fac, omega_c_max);

%%
%Impostazioni per il gradino: imposto un impulso di ampiezza W
stepOption = stepDataOptions('StepAmplitude', tab.W);

%Margini di visualizzazione del diagramma.
omega_plot_min=10^(-2);
omega_plot_max=10^5;

%%
%--Linearizzazione del sistema--

%Ridefinisco la seconda equazione di stato:
%x_dot_2 = x_1 -(C_d*u+R_0) x_2*|x_2|
%x_dot_2 = x_1 -(R_0*|x_2|) x_2 - (C_d*x_2*|x_2|) u

%La derivata di x_2 in x_dot_2 non è immediata:
%d/dx_2 (x_1 -(C_d*u+R_0) x_2*|x_2|) = -(2(x_2)^2(C_d*u+R_0))/|x_2|

%Ridefinisco y
%y= -(eta*x_2) x_1
%y= -(eta*x_1) x_2

%Punto di equilibrio dove linearizzare il sistema:
%x_equilibrio = (10,6)
%u_equilibrio = (-2.265654245335609)

%Applico la serie di Taylor per linearizzare il sistema:
%A = 0 + d/dx (x_dot)|x=x_equilibrio, u=u_equilibrio
%B = 0 + d/du (x_dot)|x=x_equilibrio, u=u_equilibrio
%C = 0 + d/dx (y)|x=x_equilibrio, u=u_equilibrio
%D = 0 + d/du (y)|x=x_equilibrio, u=u_equilibrio

%Definiamo le matrici A,B,C,D derivabili dalla forma di stato e dall'uscita

%A è una 3x3 perchè devo moltiplicare per le tre equazioni di stato 3x1 e 
%deve saltare fuori un 3x1 quindi 3x3 * 3x1 = 3x1
A = [0, 0,                                                                                    ;
     1, -(2*(tab.x_equilibrio_2^2)*(tab.C_d*tab.u_equilibrio+tab.R_0)/abs(tab.x_equilibrio_2))];

%B è una 3x1 perchè deve moltiplicare per l'ingresso 1x1 e deve saltare
%fuori una 3x1 quindi 3x1 * 1x1 = 3x1
B = [0;
     -(tab.C_d*tab.x_equilibrio_2*abs(tab.x_equilibrio_2))];

%C è 1x3 perchè 1x3 * 3x1 = 1x1
C = [-(tab.eta*tab.x_equilibrio_2), -(tab.eta*tab.x_equilibrio_1)];

D = 0;

%Per verificare il contenuto delle variabili globali uso la funzione disp
%disp(A);

%Definisco la funzione di traferimento
s=tf('s');
[NumG,DenG]=ss2tf(A,B,C,D);
G=tf(NumG,DenG);

%Ricavo il bode della G da plottare
[mag_G,phase_G,omega_G]=bode(G,{omega_plot_min,omega_plot_max});

%Stampo G
zpk(G)

%Chiusura di G senza il regolatore
CG = G/(1+G);

%Ricavo il bode della G chiusa da plottare
[mag_CG,phase_CG,omega_CG]=bode(CG,{omega_plot_min,omega_plot_max});

%Stampo la risposta al gradino di G in closed loop
figure(1);
step(CG, stepOption);
title(sprintf("Risposta al gradino (W=%d) di G in anello chiuso", tab.W));
legend("CG");
grid on;

%Informazioni sullo step
%Simulo di nuovo lo step ma in questo caso non plotto ma ricavo i dati:
[Y_CG,T_CG] = step(CG, stepOption);
%Imposto un vincolo dell'1% sul tempo di assestamento e ricavo le info:
CG_stepinfo = stepinfo(Y_CG, T_CG,'SettlingTimeThreshold',0.01);
disp(CG_stepinfo);

%Si intravede un piccolo errore a regime trascurabile. Non sappiamo se
%attenua il rumore di misura, non presente nella funzione step.

%Plot del diagramma di bode
%Nuova finestra grafica
figure(2);

%Vincolo sulla omega_c_min
patch([omega_plot_min,omega_c_min,omega_c_min,omega_plot_min],[-200,-200,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%Indico la frequenza di attraversamento minima
hold on;
text(omega_plot_min*10,-100, sprintf('w_c^*>=%.2f rad/sec', omega_c_min));

% %Vincolo sulla omega_c_min opzionale
% hold on;
% patch([omega_c_min,omega_c_min_fac,omega_c_min_fac,omega_c_min],[-200,-200,0,0], [0.9100, 0.4100, 0.1700] ,'FaceAlpha',0.3,'EdgeAlpha',0); 

%Vincolo sulla omega_c_max
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[200,200,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%Indico la frequenza di attraversamento massima
hold on;
text(omega_c_max*5,60, sprintf('w_c^*<=%.2f rad/sec', omega_c_max));

%Vincolo sull'attenuazione di n
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[-B_n_db,-B_n_db,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%plotto G
hold on;
margin(mag_G,phase_G,omega_G);

%Vincolo sul margine di fase: -180° + arg(L(jw_c))
hold on;
%Coppie di punti (omega_c_min, -180+Mf), (omega_c_max, -180+Mf), 
%(omega_c_max, -270), (omega_c_min, -270)
patch([omega_c_min,omega_c_max,omega_c_max,omega_c_min],[-180+Mf,-180+Mf,-180,-180],'red','FaceAlpha',0.2,'EdgeAlpha',0); 


%%
%--Progettazione della rete regolatrice statica--

%Ho bisogno di un polo per il vincolo numero 1 (e_inf = 0)
%Il guadagno statico resta libero: verrà modificato se necessario.
%Aggiungo un altro polo perchè ho bisogno di più pendenza
R_s = 1/s;
G_e = R_s*G;

%Stampo gli zeri e i poli di G_e 
zpk(G_e)

%Ricavo i dati sulla G_e
[mag_G_e,phase_G_e,omega_G_e]=bode(G_e,{omega_plot_min,omega_plot_max});

%Finestra precedente
figure(2);
hold on;
margin(mag_G_e,phase_G_e,omega_G_e);
legend("G",  "Vincoli", "G con Rs");
grid on;
title("Bode di G e di G con regolatore statico");

%Dal grafico si deduce che siamo caduti in uno scenario B

%%
%--Progettazione della rete regolatrice dinamica--

%Prendo una omega nel range di attraversamento
omega_c_star_ant = 800;
%Richiedo un margine di fase maggiore per robustezza
Mf_ant = Mf+10;

%Ricavo i dati di attraversamento di G_e
[Mag_G_e_omega_c_star_ant,phase_G_e_omega_c_star_ant,omega_c_star_ant]=bode(G_e, omega_c_star_ant);

M_star_ant = 1/Mag_G_e_omega_c_star_ant;
phi_star_ant = Mf_ant-180-phase_G_e_omega_c_star_ant;

tau_rete_ant = (M_star_ant-cos(phi_star_ant*pi/180))/(omega_c_star_ant*sin(phi_star_ant*pi/180));
tau_alpha_rete_ant = (cos(phi_star_ant*pi/180)-1/M_star_ant)/(omega_c_star_ant*sin(phi_star_ant*pi/180));
alpha_rete_ant = tau_alpha_rete_ant/tau_rete_ant;

%Checks
%M_star_ant > 1, 0 < phi_star_ant < 90 gradi, 
%cos(phi_star_ant*pi/180) > 1/M_star_ant

if M_star_ant <= 1
    error("M* <= 1 (%f)", M_star_ant);
elseif phi_star_ant < 0
    error("phi* < 0 (%f)", phi_star_ant);
elseif phi_star_ant > 90
    error("phi* > 90 (%f)", phi_star_ant);
elseif cos(phi_star_ant*pi/180) <= 1/M_star_ant
    error("cos(phi*) <= 1/M* (%f <= %f)", cos(phi_star_ant*pi/180), 1/M_star_ant);
end


%Ricavo il regolatore dinamico (anticipatore).
R_d_ant = (1+s*tau_rete_ant)/(1+tau_alpha_rete_ant*s);

G_e_1 = R_d_ant * G_e;


%Stampo il grafico per vedere in che scenario siamo caduti.
%Ricavo i dati sulla G_e_1
[mag_G_e_1,phase_G_e_1,omega_G_e_1]=bode(G_e_1,{omega_plot_min,omega_plot_max});

%Nuova finestra grafica
figure(3);

%Vincolo sulla omega_c_min
patch([omega_plot_min,omega_c_min,omega_c_min,omega_plot_min],[-200,-200,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%Indico la frequenza di attraversamento minima
hold on;
text(omega_plot_min*10,-100, sprintf('w_c^*>=%.2f rad/sec', omega_c_min));

% %Vincolo sulla omega_c_min opzionale
% hold on;
% patch([omega_c_min,omega_c_min_fac,omega_c_min_fac,omega_c_min],[-200,-200,0,0], [0.9100, 0.4100, 0.1700] ,'FaceAlpha',0.3,'EdgeAlpha',0); 

%Vincolo sulla omega_c_max
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[200,200,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%Indico la frequenza di attraversamento massima
hold on;
text(omega_c_max*5,60, sprintf('w_c^*<=%.2f rad/sec', omega_c_max));

%Vincolo sull'attenuazione di n
hold on;
patch([omega_plot_max,omega_c_max,omega_c_max,omega_plot_max],[-B_n_db,-B_n_db,0,0],'red','FaceAlpha',0.3,'EdgeAlpha',0); 

%Plotto G_e
hold on;
margin(mag_G_e,phase_G_e,omega_G_e);

%Vincolo sul margine di fase: -180° + arg(L(jw_c))
hold on;
%Coppie di punti (w_c_min, -180+Mf), (w_c_max, -180+Mf), (w_c_max, -180),
%(w_c_min, -180)
patch([omega_c_min,omega_c_max,omega_c_max,omega_c_min],[-180+Mf,-180+Mf,-180,-180],'red','FaceAlpha',0.2,'EdgeAlpha',0);

%Plotto G_e_1
hold on;
margin(mag_G_e_1,phase_G_e_1,omega_G_e_1);

legend("G con Rs",  "Vincoli", "G con Rs*Rd");
grid on;
title("Bode di G con i vari regolatori");

%Il margine di fase può essere recuperato:
%Siamo caduti in uno scenario A.

%COSA SUCCEDE IN ANELLO APERTO:
zpk(G_e_1)
%Il polo -3.333 viene cancellato con lo zero della rete (30% di scarto con
%relativa coda di assestamento), mentre il polo della rete è trascurabile.
%Rimane il polo nell'origine. Se collegassi G_e_1 in anello aperto
%funzionerebbe come integratore, almeno per frequenze basse.

%COSA SUCCEDE IN ANELLO CHIUSO (G_e_1/(1+G_e_1)): 
zpk(G_e_1/(1+G_e_1))
%Mi rimane solo una coppia di poli complessi coniugati, oltre ad un
%guadagno molto elevato.

%--------------------------------------------------------------------------
%Soluzione: modifico il guadagno.

%Settare questa opzione per aggiungere un margine di attenuazione ulteriore
%Si tratta di un offset negativo.
mu_d_offset_db = 1;

%Ricavo il valore del modulo di G_e_1 in omega_n
[mag_G_e_1_omega_n,phase_G_e_1_omega_n,omega_n]=bode(G_e_1, tab.omega_n);
mag_G_e_1_omega_n_db = 20*log10(mag_G_e_1_omega_n);

%Devo imporre un valore di mu tale che:
%Mag_G_e_1_omega_n_db+mu_d_db=-B_n_db-mu_d_offset_db:

mu_d_db=-B_n_db-mu_d_offset_db-mag_G_e_1_omega_n_db;
mu_d_1 = 10^(mu_d_db/20);

mu_d = mu_d_1;

L = mu_d * G_e_1;

%Plotto L
figure(3);
[mag_L,phase_L,omega_L]=bode(L,{omega_plot_min,omega_plot_max});
hold on;
margin(mag_L,phase_L,omega_L);

legend("G con Rs", "Vincoli", "G con Rs*Rd", "L");
grid on;
title("Bode di G con i vari regolatori fino a L");

%Il grafico sembra soddisfare le specifiche.

%NOTA IMPORTANTE: sembrerebbe che dal grafico non vengano soddisfatte le
%specifiche di tempo di assestamento facoltative. In realtà bisogna tenere
%in conto la chiusura in retroazione di L:
%il margine di fase nell'attraversamento assume un altro valore e di
%conseguenza anche omega_c_min si sposta: andando a plottare i vincoli
%sulla F si nota che viene rispettato anche questo vincolo.

%Calcolo R
R = mu_d * R_d_ant * R_s;

%Dati usati in simulink per la R.
[NumR, DenR] = tfdata(R);
NumR = NumR{1,1};
DenR = DenR{1,1};

%L già calcolata

%Dati usati in simulink per la L.
[NumL, DenL] = tfdata(L);
NumL = NumL{1,1};
DenL = DenL{1,1};


%%
%Chiusura dei loop

F=L/(1+L);


%Ricavo informazioni
[mag_F,phase_F,omega_F]=bode(F,{omega_plot_min,omega_plot_max});

%Stampo gli zeri e i poli di F
zpk(F)

%Plotto la risposta a gradino
figure(5);
step(F, stepOption);
title(sprintf("Risposta al gradino (W=%d) di L in anello chiuso", tab.W));
legend("F");
grid on;

%Confronto la risposta a gradino con quella della closed loop di G.
%Avevo già plottato in figura 1.

%Per muovere il cursore:
%datacursormode on

%Informazioni sullo step
%Simulo di nuovo lo step ma in questo caso non plotto ma ricavo i dati:
[Y_F,T_F] = step(F, stepOption);
%Imposto un vincolo dell'1% sul tempo di assestamento e ricavo le info:
F_stepinfo = stepinfo(Y_F, T_F,'SettlingTimeThreshold',0.01);
disp(F_stepinfo);

%Dallo stepinfo abbiamo un Tempo di assestamento di 0.0321 sec e una
%sovraelongazione percentuale dello 0%, valori al di sotto dei vincoli:
%inoltre la specifica opzionale di tempo di assestamento è stata risolta.

%Rappresento con il diagramma di bode F e CG
figure(6);

margin(mag_F, phase_F, omega_F); 
hold on;
margin(mag_CG, phase_CG, omega_CG); 

title("Confronto tra L e G in anello chiuso con Bode");
legend("L closed loop", "G closed loop");
grid on;


%Dal grafico si osserva che:
%Il vincolo di misura viene rispettato: -30db circa in 1000 rad/s
%Il vincolo di sovraelongazione viene rispettato: alla pulsazione di taglio
%ho un valore di -70 gradi.
%Non riesco a capire se il vincolo di tempo di assestamento viene
%rispettato, lo verifico dalla risposta a gradino nel grafico precedente.
%L'errore a regime è nullo, sempre verificato dal grafico precedente.

%%
%CONCLUSIONI

%Specifiche F:
%Tempo assestamento: 0.0321 sec
%Sovraelongazione percentuale: 0%
%Attenuazione errore di misura: 30db
%Errore e_inf = 0
%Margine di fase: 70 gradi

%Specifiche G in anello chiuso:
%Tempo assestamento: 0.000515 sec
%Sovraelongazione percentuale: 0%
%Attenuazione errore di misura: -0.2db
%Errore e_inf non nullo (ma molto basso: e<0.03 verificato con simulink)
%Margine di fase: 135 gradi

%In conclusione il sistema G senza regolatore chiuso in retroazione
%funziona molto più velocemente in quanto è approssimabile ad un sistema
%del tipo: G(s)=mu nel range di frequenze fino a 4000 rad/s circa. L'unica
%pecca è la violazione sui vincoli di e_inf=0 e attenuazione di misura.
%Se non ci fosse stato il problema di questi vincoli, allora il sistema G
%in anello chiuso sarebbe bastato, a meno di una correzione del guadagno.

%%
%PROVA IN SIMULINK
%open('progetto_simulink_new.slx')

%Il regolatore sul sistema linearizzato risponde in modo impeccabile
%Il regolatore sul sistema non lineare funziona piuttosto bene in un
%intorno della coppia di equilibrio.

