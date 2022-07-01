function [massa_pod_irr, dati_trave, massa_emp, dati_longherone, Massa_ala, Massa_strutt_complessiva] = struct_eval(dist_payload, guess_tail, Sh, diagV, alfa_rombo, co, y0, max_sp, S)
% INPUT
% dist_payload: distribuzione del payload, 9 sacchette da 300g si
% descrivono con un vettore [3 1 9] (3 lungo il body, 9 in altezza)
% guess_tail: lunghezza della trave di coda
% Sh: superficie dello stabilizzatore
% diagV: diagonale verticae del rombo (lungo la fusoliera)
% co: vettore con corde dell'ala (importa solo il primo elemento, ovvero
% la corda di radice
% y0: posizione delle corde specificate in co (importa solo il terzo, in
% questo caso ultimo, elemento). es. [0 1 1.3] in questo caso il dato
% importante è l'1.3 da cui lo script ricava l'apertura
% max_sp: massimo spessore percentuale del profilo dell'ala (0.30 valore
% sensato)
% S superficie alare

% OUTPUT
% valori importanti:
% massa_pod (o fusoliera)
% massa_emp: massa impennaggi (orizzontali + verticali)
% Massa_ala
% Dovrebbero essere tutti in grammi


% Lo script procende dimensionando consequenzialmente "pod", tubo di coda
% e ala. Le ipotesi di fondo sono che la massa del pod tiene conto solo
% dell'ingombro del payload, che il carico sollecitante il tubo di coda
% tenga conto solo della massa di payload e pod e che il carico
% sollecitante l'ala tenga conto di masa di payload, pod e tubo di coda (da
% qui il dimensionamento "a cascata")

%_________FUSOLIERA(POD)___________

% Stima della massa della fusoliera
% Questo script serve a dare una prima approssimazione della massa della
% fusoliera a partire dal numero di sacchette trasportate e dalla loro
% disposizione. I dati in unput sono l, b e h che sono rispettivamente il
% numero di sacchette ripetute nella direzione della lunghezza, della
% larghezza e dell'altezza. La disposizione identificata con la terna 3,2,4
% ad esempio indica un totale di 24 sacchette disposte in 2 file parallele
% da 3 "vagoni" ripetute su 4 livelli in altezza. Si assume di utilizzare
% le sacchette da 300g. Viene assunta una configurazione costruttiva a
% guscio in sandwich con facce in textreme da 160g/m2 e depron da 40g/m2

d_tex= 160 ; %g/m2
d_depr= 40 ; %g/m2
d_sand= d_depr + 2*2*d_tex ;
d_stiff= 120 ;  %irrigidimenti longitudinali= 120g/m
massa_sacchetta= 100;  % [grammi]
l_s=250 ;  %lunghezza sacchetta
b_s=115 ;  %larghezza sacchetta
h_s=15  ;  %altezza acchetta

v= dist_payload;            %vettore del numero e disposizione delle sacchette - INPUT
n_sacchette= v(1)*v(2)*v(3); % numero delle sacchette
sup= v(2)*b_s*v(1)*l_s*2 + (v(2)*b_s+v(1)*l_s)*2*v(3)*h_s ;
vol=v(1)*l_s * v(2)*b_s  *v(3)*h_s ;

l_dim= v(1)*l_s ; %lunghezza (dimensionale) [mm]
b_dim= v(2)*b_s ; %largehzza (dimensionale) [mm]
h_dim= v(3)*h_s ; %altezza (dimensionale)   [mm]
sez_front= b_dim*h_dim; %mm2
sup_vol_ratio= sup/vol;

massa_pod= sup * 10^-6 * d_sand + v(1)*l_s * 10^-3 * d_stiff; %risultato in grammi
massa_pod_irr= massa_pod + d_stiff * l_dim * 0.001;
Wcarr = 0.256 + 0.2 + 0.011 + 0.02;
Weln = 0.367 + 0.0109 + 0.010*4 + 0.177;
massa_pod = massa_pod_irr + (Wcarr + Weln)*1000;





%% ____________ TRAVE DI CODA _________________


Res_Materiale=1100000000;                 %[N/m^2]
Dens_Materiale=1750;                    %[Kg/m^3]
E=300000000000;                         %[N/m^2]
Sf_t=3;
L=guess_tail ;                 %[m] lunghezza trave - INPUT
s_t=0.0008:0.0008:0.0032;                 %[m] spessore parete della trave
d=0.01:0.001:0.2;                       %[m] diametro trave                               %[Kg]
P= 3.5* 9.81* 0.001* (n_sacchette*massa_sacchetta + massa_pod);                           %[N]
I_lim_t=10*P*L^2/(3*E);

Mf_t=L*P;                                 %[Nm]
%figure(1)
for i=1:4
      Sigma(i,:)=Mf_t.*d.*32./(pi.*(d.^4-(d-2*s_t(i)).^4));
      %plot(d,Sigma(i,:))  %sigma=f(diametro) con vari spessori
      %hold all
end
%plot(d,Res_Materiale/Sf_t*ones(size(d))) %valore di soglia (resistenza del materiale)
%hold all

for i=1:4
    f=1;
    for j=1:191
        if (Sigma(i,j)<Res_Materiale/Sf_t && f==1)
            d_res(i)=d(j);                       %[m]
            f=0;
        end
    end    
end
I_t=(pi.*(d_res.^4-(d_res-2*s_t(i)).^4))/64;                %[m^4]


Massa_t=L.*pi.*(d_res.^2-(d_res-2*s_t).^2).*Dens_Materiale/4; %[kg]

for i=1:4
    if (I_t(i)<I_lim_t)
        Massa_t(i)=1000;
    end
end

d_final_trave=d_res(Massa_t==min(Massa_t)); %OUTPUT
s_final_trave=s_t(Massa_t==min(Massa_t));
Massa_trave=min(Massa_t); % OUTPUT
Massa_trave_g= Massa_trave*1000; %[grammi]
dati_trave = [d_final_trave s_final_trave Massa_trave]; %m m kg



%% __________ IMPENNAGGI _________
%Per renderequesta sezione universale per tutte le configurazioni di
%impennaggi è necessario inserire in sup_imp_wet la somma delle sperfici
%bagnate e la somma delle aperture (ad esempio somma dell'apertura di
%imennaggio orizzontale e verticale o somma dell'apertur dei due impennaggi
%di una coda a V)


dens_skin_emp= 360; % densità superficiale rivestimento depron+vetro [g/m2]
dens_spar= 60; % densità lineare tondinio in carbonio pultruso diam esterno 10mm, interno 7mm [g/m]
sup_imp_wet=2*Sh; % superficie bagnata totale degli impennaggi [m^2] -   INPUT
delta_x = diagV/2 - guess_tail;
beta = 180 - alfa_rombo;
bh = 3*delta_x*tand(beta/2);
L_imp_tot= bh; %somma delle aperture degli impennaggi [m]

massa_emp= dens_skin_emp * sup_imp_wet + dens_spar * L_imp_tot; %[grammi] % OUTPUT





%% ________ LONGHERONE ________



Res_Materiale=1100000000;           %[N/m^2]
Sf_l=1.5;
L_ala=y0(3)*2;                            %Apertura alare [m]    % INPUT
Dens_Materiale=1750;                %[Kg/m^3]
Dens_balsa= 160;                    %[kg/m^3]


Mf_l= Sf_l* 9.81* 0.001* (n_sacchette*massa_sacchetta + massa_pod+ Massa_trave + massa_emp) * (L_ala * 0.5 );  %[Nm]
Mod_Res=Mf_l/(Res_Materiale*Sf_l);      %[m^3]
I_lim_l=10*Mf_l^2/(48*E);

c=co(1);                              %[m]      % INPUT (radice)
h=max_sp;                             %[m]      % INPUT (spessore massimo radice)  
b=0.001:0.001:0.02;                 %[m]
s_l=0.001:0.001:h/2;
h_real=h*(1-b/c);
I=Mod_Res.*h_real;

%figure(2)
for i=1:20
    Sol(i,:)=2/3.*s_l.^3.*b(i)-s_l.^2.*b(i).*h_real(i)+s_l.*b(i).*h_real(i).^2/2-I(i);
    %plot(s_l,Sol(i,:))
     % hold all
end
for i=1:20
    f=1;
    for j=1:length(s_l)
        if (Sol(i,j)>0 && f==1)
            s_res(i)=s_l(j);                       %[m]
            f=0;
        end
    end    
end
I_l=2.*(b.*s_res.^3./12)+2.*(b.*s_res).*(h_real./2-s_res./2).^2; 
Massa_l=2.*s_res.*b.*L_ala.*Dens_Materiale;

for i=1:length(h_real)
    if (I_l(i)<I_lim_l)
        Massa_l(i)=1000;
    end
end

h_final_longherone=h_real(Massa_l==min(Massa_l)); % OUTPUT
b_final_longherone=b(Massa_l==min(Massa_l)); % OUTPUT
s_final_longherone=s_res(Massa_l==min(Massa_l)); % OUTPUT

Massa_longherone=min(Massa_l);
Massa_longherone_g=Massa_longherone*1000; %[g]

vol_core=L_ala*b_final_longherone(1)*(h_final_longherone(1)-2*s_final_longherone(1)); %[m^3] % OUTPUT
Massa_longherone_con_core=Massa_longherone + (Dens_balsa*vol_core)*1000;  %[grammi] % OUTPUT

dati_longherone = [h_final_longherone(1) b_final_longherone(1) s_final_longherone(1) vol_core Massa_longherone_con_core];


%% __________ MASSA ALA ____________



sup_ala_wet= 2*S;   %superficie bagnata ala [m^2]                   % INPUT
dens_skin=320;      %densità rivestimento ala depron + vetro [g/m^2]


Massa_ala= Massa_longherone_con_core + dens_skin*sup_ala_wet;  %[g] % OUTPUT


%% __________ MASSA COMPLESSIVA______

Massa_strutt_complessiva= Massa_ala + Massa_trave + massa_emp + massa_pod_irr;  %[grammi] %  % OUTPUT 

end
