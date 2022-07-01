clear all
close all

%% Geometric data based on current regulations
lato_rombo = 1.5; %metri
alfa_rombo = 60:3:80; %deg

Xcg_pos = 0.30; % 30% corda

wings_per_angle = 4;    % How many wings to generate
payload_per_wing = 6;   % How many payloads to generate


%% Caricamento file profili
% Airfoil data extraction 
% (see xfoilfetch.mat)
foils = [load('NACA0012.mat') load('clarkY.mat')  load('S1210.mat') load('NACA2412.mat') load('SD7037.mat')];
foil_name = ["NACA0012"; "clarkY"; "FX-63_100"; "FX-63_137"; "S1210"; "NACA2412"; "SD7037";];

% Max thickness
v_spess = [0.12 0.10 0.12 0.12 0.09];

% Fuselage drag
load('fitCDpod.mat') 


%% Main loop

l = 1;

for a = alfa_rombo 
    
    % Generate wings
    for t = 1:wings_per_angle
        Vs = [5 50]; % Speed interval to consider
        % Generate a random wing
        [c, y0, S] = wing_gen(lato_rombo, a);
        for z = 2:length(foils)
            % Loop through the airfoils
            
            airfoil_data = foils(z).airfoil_data;
            
            % Lifting Line Theory analysis setup
            Vrange = Vs(1):5:Vs(end);
            for i = 1:length(Vrange)
                cs = [c(end:-1:2) c]';      % Chord distribution
                ys = -y0(end):0.1:y0(end);  % Corresponding Y coord.
                cd = interp1(y0,c,ys);      % Interpolation of chords
                MAC = mean(cs);
                
                Re = round((1.225 * MAC * Vrange(i))/1.81e-5, -5); % Elininates some non convergence problems
                % Re in airfoil data is distretized in multiples of 10000
                
                if Re>900000
                    % Not much change abouve 900000
                    Re = 900000;
                end
                if Re == 0
                    % If Re = 0, round it to the minimum one
                    Re = 100000;
                end
                
                % Pre-allocation of memory based on deisred number of
                % planes
                alfas = airfoil_data{Re/100000,1};
                L = zeros(1,length(alfas));
                D = zeros(1,length(alfas));
                CL = zeros(1,length(alfas)); 
                CD = zeros(1,length(alfas));
                
                % Lifting Line Theory analysis (see clean_LLT.m)
                for j = 1:length(alfas) 
                    [CL(j),CD(j),L(j),D(j),E(j),Mfoc] = clean_LLT(c,y0,alfas(j),Vrange(i),airfoil_data);
                end
                wing_data(i,:) = {Vrange(i) alfas CL CD L D E Mfoc};
            end
            
            
            % Assign payload and analyise plane
            for p = 1:payload_per_wing
                Vs = [10 50];
                 % Generates a payload disposition (see gen_payload.m)
                [dist_payload, area, Wp] = gen_payload();

                % Estimates of structures weight
                pod_length = dist_payload(1)*0.25; %m
                diagV = 2 * lato_rombo * sind(a/2);
                guess_tail = diagV/2 *0.8 - pod_length/2 + c(1)*0.30;
                guess_Sh = (0.20*S)/2 *3; %including fin
                
                % Structures weight by Mario and Nicholas (see struct_eval.m)
                spess_max = v_spess(z)*c(1);
                [massa_pod_irr, dati_trave, massa_emp, dati_longherone, Massa_ala, Mtot] = struct_eval(dist_payload, guess_tail, guess_Sh, diagV, a, c ,y0, spess_max, S);
                
                % Additional masses
                Wstrut = Mtot*0.001*9.81; % g -> kg
                Wcarr = (0.256 + 0.2 + 0.011 + 0.02)*9.81; 
                Weln = (0.367 + 0.0109 + 0.018*4 + 0.177)*9.81;
                Wtot = Wp + Wcarr + Weln + Wstrut;
                
                % Equilibrium analysis [deprecated]
                %T = @(v) -0.285*v + 19.8; % mia stima
                %T = @(v) -0.5164.*v + 23.5820; %vecchia elica
                %T = @(v) -0.216.*v + 11.23;
                %T = @(v) -0.285.*v + 11.23; %misto motocalc e mia stima
                %T = @(v) -0.285.*v + 11.23;
                %T = @(v) -1.959.*v + 22;
                
                % Thrust curve from motor model
                T = @(V) -0.0037.*V.^2 -0.3683.*V + 16.9020;

                % [deprecated]
                %[Tnec,eq,wing_data] = eq_eval(c, y0, S, Wtot, Vs, airfoil_data);
                
                % First equilibrium analysis
                [Tnec,eq] = eq_eval_light(c, y0, S, Wtot, Vs, airfoil_data, wing_data, zeros(1,length(Vrange)), zeros(1,length(Vrange)));

                % Add drag of fuselage and calculate necessary thrust
                % Fit some data calculated with OpenVSP
                CD_pod = @(vv,area) sfitCD(vv,area);
                Dspod = 0.5*1.225.*(Vs(1):5:Vs(end)).^2 .*CD_pod(Vs(1):5:Vs(end),area) *area;

                % Necessary thrust
                ctnec = polyfit(Vs(1):5:Vs(2), Tnec+Dspod,2);
                tnec = @(v) ctnec(1).*v.^2 + ctnec(2).*v + ctnec(3);
                
                % Speed and CL of equilibrium by crossing necessary and
                % available power
                d = @(x) tnec(x)-T(x);
                Veq = fzero(d,40);
                CLeq = Wtot/(0.5*1.225*Veq^2 *S);
                

                % Tail evaluation
                
                VH = 0.5; % Tail volume
                
                % Tail Geometry (see find_tail.m)
                [tail_arm_store, Sh, i_tail, bh] = find_tail(VH, S, c, y0, a, lato_rombo, wing_data, Veq, Wtot);
                ch = Sh/bh;
                dts = 1;
                for dt = Vs(1):5:Vs(end)
                    [CLh,CDh,Lh(dts),Dtoth(dts),E,Mfoch] = clean_LLT([ch ch ch], [0 bh/4 bh/2], i_tail, dt, foils(1).airfoil_data);
                    dts=dts+1;
                end
                
                %% Redo analysis with "real" weights
                
                % Same as before
                [massa_pod_irr, dati_trave, massa_emp, dati_longherone, Massa_ala, Mtot] = struct_eval(dist_payload, tail_arm_store, Sh, diagV, a, c ,y0, spess_max, S);

                Wstrut = Mtot*0.001*9.81;
                Wcarr = (0.256 + 0.2 + 0.011 + 0.02)*9.81;
                Weln = (0.367 + 0.0109 + 0.018*4 + 0.177)*9.81;
                Wtot = Wp + Wcarr + Weln + Wstrut;
                
                [Tnec,eq] = eq_eval_light(c, y0, S, Wtot, Vs, airfoil_data, wing_data, Dtoth, Lh);

                CD_pod = @(vv,area) sfitCD(vv,area);
                Dspod = 0.5*1.225.*(Vs(1):5:Vs(end)).^2 .*CD_pod(Vs(1):5:Vs(end),area) *area;

                %tnec = @(v) interp1(Vs(1):5:Vs(2),Tnec+Dspod,v);
                ctnec = polyfit(Vs(1):5:Vs(2), Tnec+Dspod,2);
                tnec = @(v) ctnec(1).*v.^2 + ctnec(2).*v + ctnec(3);

                % New equilibrium
                d = @(x) tnec(x)-T(x);
                Veq = fzero(d,40);
                CLeq = Wtot/(0.5*1.225*Veq^2 *S);
                
                % New tail
                VH = 0.5;
                [tail_arm_store, Sh, i_tail, bh] = find_tail(VH, S, c, y0, a, lato_rombo, wing_data, Veq, Wtot);
                ch = Sh/bh;
                
                %% Evaluation of competition flight

                Vs = Vs(1):5:Vs(2);
                Clseq = Wtot./(0.5*1.225.*Vs.^2 *S);
                Cdseq = Tnec./(0.5*1.225.*Vs.^2 *S);
                E = Clseq./Cdseq;
                E(E(:)==Inf) = 0; % Removes non-equilibrium cases
                [Emax,n] = max(E);
                VEmax = Vs(n);
                
                % Power for climb evaluation
                Pns =Tnec.*Vs;
                Pd = @(v) -0.5164.*v.^2 + 23.5820.*v;
                cpn = polyfit(Vs,Pns,2);
                Pn = @(v) cpn(1).*v.^2 + cpn(2).*v + cpn(3);
                Vs = Vs(1):1:Vs(end);
                dPmax =  max(Pd(Vs)-Pn(Vs));

                RofC = dPmax/Wtot; % Rate of climb
                height = RofC*60;
                
                % Maximum climb check
                if height > 100
                    height = 100;
                end
                
                % Takeoff distance
                CLto = 1.2;
                CD0 = min(wing_data{1,4});
                CDpodto = CD_pod(10,area);
                D=@(V) 0.5*1.225*S*CD0*V.^2 + 0.5*1.225*area*CDpodto*V.^2;
                mi = 0.05;
                K = Wtot*mi;
                m=Wtot/9.81;
                dS=@(V) V*m./(T(V)-D(V)-K);
                Vmin=1.2.*sqrt(2*Wtot/(1.225*CLto*S));
                
                % [DEPRECATED]

%                 Vi=linspace(1,15,100);
%                 lto=zeros(1,length(Vi));
%                 for i=1:length(Vi)
%                     Vto=Vi(i);
%                     lto(i)=integral(dS,0,Vto);
%                 end
%                 g = 9.81;
%                 T0 = T(0); %spinta statica             
%                 c = polyfit(0:20, T(0:20), 2);
%                 mi = 0.04;
%                 CD0 = min(wing_data{1,4});
%                 CLto = 1.2;
%                 Vto = sqrt((2.*Wtot)./(S.*1.225.*Cl));
%                 A = g.*(T0./Wtot - mi);
%                 B = g./Wtot .* (0.5*rho*S.*(CD0 - c_fric.*CLto) - c(1));
%     
%                 Vx = @(v) v./(A-B.*v.^2); %funzione da integrare
%                 x = integral(Vx, 0, Vto); %integrazione numerica

                
%% Build struct with generated data

                teams(l).ID = l;
                teams(l).alfa_rombo = a;
                teams(l).Surface = S;
                teams(l).S_h = Sh;
                teams(l).i_tail = i_tail;
                teams(l).tail_arm = tail_arm_store;
                teams(l).bh = bh;
                teams(l).c = c;
                teams(l).y0 = y0;
                teams(l).foil_name = foil_name(z);
                teams(l).Wpayload = Wp;
                teams(l).Wtot = Wtot;
                teams(l).Wstructures = Wstrut;
                teams(l).Wdry = Wtot-Wp;
                teams(l).CLeq = CLeq;
                teams(l).payload = dist_payload;
                teams(l).height = height;
                teams(l).PSh = -3.92e-5 * height^4 + 1.08e-2 * height^3 -1.156 * height^2 + 64.2 * height -537;
                teams(l).Vmax = Veq;
                teams(l).Emax = Emax;
                %teams(l).X_to = lto(round(Vmin));
                teams(l).distance = Veq * 120;
                teams(l).Sh = 0;
                teams(l).Sp = 0;
                teams(l).Sd = 0;
                teams(l).S = 0;
                teams(l).massa_pod_irr = massa_pod_irr;
                teams(l).dati_trave = dati_trave;
                teams(l).massa_emp = massa_emp;
                teams(l).dati_longherone = dati_longherone;
                teams(l).Massa_ala = Massa_ala;
                disp(l)
                l = l+1;
            end
        end
    end
end

%% SCORING
% Follows the algorithm in the regulations

Pmax = 0;
Dmax = 1;
PSmax = 1;

for i = 1:l-1
    if i==1
        Pmax = teams(i).Wpayload;
        Dmax = teams(i).distance;
        teams(i).Sp = 1000;
        teams(i).Sd = 1000;
        teams(i).height = 100;
        teams(i).Sh = 1000;
        teams(i).S = (1000 + 1000 + teams(i).Sh)/3;
    else
        if teams(i).Wpayload > Pmax
            Pmax = teams(i).Wpayload;
            teams = ricalcolo_p(teams, Pmax, i);
            teams(i).Sd = 1000*teams(i).distance/Dmax;
            teams(i).Sh = 1000*teams(i).PSh/PSmax;
            teams(i).S = (teams(i).Sp + teams(i).Sd + teams(i).Sh)/3;
        elseif teams(i).distance > Dmax
            Dmax = teams(i).distance;
            teams = ricalcolo_d(teams, Dmax, i);
            teams(i).Sp = 1000*teams(i).Wpayload/Pmax;
            teams(i).Sh = 1000*teams(i).PSh/PSmax;
            teams(i).S = (teams(i).Sp + teams(i).Sd + teams(i).Sh)/3;
        elseif teams(i).height > PSmax
            PSmax = teams(i).PSh;
            teams = ricalcolo_h(teams, PSmax, i);
            teams(i).Sd = 1000*teams(i).distance/Dmax;
            teams(i).Sp = 1000*teams(i).Wpayload/Pmax;
            teams(i).S = (teams(i).Sp + teams(i).Sd + teams(i).Sh)/3;
        else
            teams(i).Sp = 1000*(teams(i).Wpayload/Pmax);
            teams(i).Sd = 1000*(teams(i).distance/Dmax);
            teams(i).Sh = 1000*teams(i).PSh/PSmax;
            teams(i).S = (teams(i).Sp + teams(i).Sd + teams(i).Sh)/3;
        end
        teams(i).S = (teams(i).Sp + teams(i).Sd + teams(i).Sh)/3;
    end
end

% Convert to table
teams = struct2table(teams);

%% Support functions
function [t] = ricalcolo_p(teams, Pmax, n)
for i=1:n
    teams(i).Sp = 1000*teams(i).Wpayload/Pmax;
    teams(i).S = (teams(i).Sp + teams(i).Sd + teams(i).Sh)/3;
    
end
t=teams;
end

function [t] = ricalcolo_d(teams, Dmax, n)
for i=1:n
    teams(i).Sd = 1000*teams(i).distance/Dmax;
    teams(i).S = (teams(i).Sp + teams(i).Sd + teams(i).Sh)/3;
end
t=teams;
end

function [t] = ricalcolo_h(teams, PSmax, n)
for i=1:n
    teams(i).Sh = 1000*teams(i).PSh/PSmax;
    teams(i).S = (teams(i).Sp + teams(i).Sd + teams(i).Sh)/3;
end
t=teams;
end

