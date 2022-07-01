function [tail_arm_store, Sh, i_tail, bh] = find_tail(VH, S, c, y, alfa_romb, lato_rombo, wing_data, Veq, Wtot)

tol_AR = 0.05;
target_AR = 6;

ch = 0.12; % m

root = c(1);
c = [c(end:-1:2) c]'; % Chords on entire wing
c1=c;
y = [-y(end:-1:2) y]';
y0 = linspace(y(1),y(end),15);
c = interp1(y,c,y0); % Chords interpolation
MACw = mean(c);
diagV = 2 * lato_rombo * sind(alfa_romb/2); 
diagH = 2 * lato_rombo * cosd(alfa_romb/2);
tail_arm = (diagV/2)*0.20; % Initial guess
err_AR = 1;
cont = 1;
while abs(err_AR) > tol_AR && cont<50
    close all
    Sh = VH*MACw*S/tail_arm;
    V = Veq - mod(Veq,5);
    CLeq = Wtot/(0.5*1.225*Veq^2 * S);
    if V < 10
        V = 10;
    end
    M = wing_data{[wing_data{:,1}]==V,8};
    Cm = M/(0.5*1.225*V^2 * S*MACw);
    CLh = (Cm + CLeq*(0.30-0.25))/VH;
    CL_alfa_h = 0.11; % [1/deg] assuming 2*pi approximation
    
    % Horizontal tail aingle
    i_tail = CLh/CL_alfa_h;
    
    % Arms calculation
    beta = 180 - alfa_romb;
    x_TE_h = tail_arm - (0.30-0.25)*root + 0.75*ch;
    delta_x = diagV/2 - x_TE_h;
    
    % Tail span
    bh = 2*delta_x*tand(beta/2);
    AR = bh^2 / Sh;
    tail_arm_store = tail_arm;
    err_AR = (target_AR - AR)/target_AR; %errore relativo su AR
    tail_arm = tail_arm - tail_arm*err_AR*0.05;
    cont = cont+1;
end

% Plane drawing (optional)

% figure(2)
% rombo = polyshape([diagH/2 0 -diagH/2 0],[0 diagV/2 0 -diagV/2]);
% plot(rombo, 'FaceAlpha', 0.1, 'FaceColor', 'w');
% pbaspect([1 1 1]);
% axis([-1.5 1.5 -1.5 1.5]);
% hold on
% cs = [c1(1)/2 0.30*c1(2)+0.20*c1(1) 0.30*c1(3)+0.20*c1(1)];
% cs = [cs cs(2:-1:1) -c1(1)/2 -(c1(2)-(0.30*c1(2)+0.20*c1(1))) -(c1(3)-(0.30*c1(3)+0.20*c1(1))) -(c1(2)-(0.30*c1(2)+0.20*c1(1))) -c1(1)/2];
% ys = [y(1) y(2) y(3) y(4) y(5)];
% ys = [ys ys(end:-1:1)];
% ala = polyshape(ys,cs);
% plot(ala)
% plot(ys(1:5),0.20.*c1)
% 
% ch = Sh/bh;
% 
% stab = polyshape([-bh/2 bh/2 bh/2 -bh/2],[(-x_TE_h+ch) (-x_TE_h+ch) -x_TE_h -x_TE_h]);
% plot(stab)

end