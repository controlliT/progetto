%%Introduzione

%Progetto di Controlli Automatici - T
%Tipologia III variante A: Controllo di un sistema idroelettrico con condotta forzata

%Gruppo:
%Andrea Proia 0000825784
%Federico Maria Macchiavelli 0000825621
%Mattia Innocenti 0000825046
%Luca Bartolomei 0000825005


%%Caratteristiche impianto

%Si tratta di un sistema SISO.

%Tabella delle caratteristiche:
tab_a = pi^2;
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

%Definizione del sistema:
%Il sistema è formato da tre variabili di stato:

%x_1: la pressione dell’acqua sul fondo del bacino.
%x_2: la portata in uscita dal bacino.
%x_3: la portata in ingresso allo stesso.

%L'uscita del sistema viene rappresentata dalla variabile y:
%y: rappresenta l’energia elettrica generata attraverso la turbina.
%y=?eta*x_1*x_2

%Forma di stato:
%x_dot_1 = (g/a)*(-x_2+x_3)
%x_dot_2 = -C_d*u*x_2*|x_2|-R_0*x_2*|x_2|+x_1
%x_dot_3 = 0
 
%%Specifiche di progetto


%%Linearizzazione del sistema

%Ridefinisco la seconda equazione:
%x_dot_2 = x_1 -(C_d*u*|x_2|+R_0*|x_2|) x_2
%x_dot_2 = x_1 -(R_0*|x_2|) x_2 - (C_d*x_2*|x_2|) u

%Ridefinisco y
%y= ?(eta*x_2) x_1
%y= ?(eta*x_1) x_2

%Punto di equilibrio dove linearizzare il sistema:
%x_equilibrio = (10,6,6)
%u_equilibrio = (-2.265654245335609)

%Definiamo le matrici A,B,C,D derivabili dalla forma di stato e dall'uscita
%Non so se sono giuste: bisogna usare simulink per esserne sicuri.

A = [0, -(tab_g_const/tab_a),                                                                (tab_g_const/tab_a);
     1, -(tab_C_d*tab_u_equilibrio*abs(tab_x_equilibrio_2)+tab_R_0*abs(tab_x_equilibrio_2)), 0;
     0, 0,                                                                                   0];
 
B = [0;
     -(tab_C_d*tab_x_equilibrio_2*abs(tab_x_equilibrio_2));
     0];

C = [-(tab_eta*tab_x_equilibrio_2), -(tab_eta*tab_x_equilibrio_1),0];

D = 0;

%Definisco la funzione di traferimento
s=tf('s');
[N,D]=ss2tf(A,B,C,D);
G=tf(N,D);
%Stampa la forma canonica di bode.
zpk(G)

w_plot_min=10^(-2);
w_plot_max=10^5;

[Mag,phase,w]=bode(G,{w_plot_min,w_plot_max});
figure()
margin(Mag,phase,w);
grid on;


