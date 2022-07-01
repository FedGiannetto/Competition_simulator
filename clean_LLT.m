function [CL,CD,L,Dtot,E,Mfoc] = clean_LLT(c, yc, alfa, Vinf, airfoil_data)

n = 60; % stations

theta0 = linspace(0,pi,n);

%% GEOMETRY AND SETUP

b = yc(end) * 2;
S = 2*trapz(yc,c);
c = [c(end:-1:2) c]'; % Chords on entire wing
y = [-yc(end:-1:2) yc]';
y0 = -b/2 * cos(theta0);
c = interp1(y,c,y0); % Chords interpolation

% First evaluation of Reynolds to choose correct CL_alpha
MAC = mean(c);
Re = round((1.225 * MAC * Vinf)/1.81e-5, -5);
if Re>900000 % Doesn't change much above this
    Re = 900000;
end
if Re == 0
        Re = 100000;
end
CL_alfa = 57.30*airfoil_data{Re/100000,5}; % Convertion to 1/rad

%fun_alfa = @(a_interp) interp1(airfoil_data{Re/100000,1},airfoil_data{Re/100000,2}, a_interp);
cc = polyfit(airfoil_data{Re/100000,1},airfoil_data{Re/100000,2},1);
fun_alfa = @(a) cc(1)*a + cc(2);
alfa_L0 = deg2rad(fzero(fun_alfa,1));

alfat = deg2rad(alfa).*ones(1,n);
CL_alfa = CL_alfa.*ones(1,n);
alfa_L0 = alfa_L0.*ones(1,n);

%% COEFFICIENTS An

B = (alfat(2:n-1)-alfa_L0(2:n-1))';
for i=n-1:-1:2 %theta0
    for j=n-1:-1:2
        A(i-1,j-1)=4*b*sin((j-1)*theta0(i))/CL_alfa(i)/c(i)+(j-1)*sin((j-1)*theta0(i))/sin(theta0(i));
    end
end
% Linear system
An = A\B;

%% OUTPUT LLT

AR = b^2/S;
CL = AR*pi*An(1);
CDi=pi*AR*sum((1:n-2)'.*(An.^2));
L = 0.5*1.225*Vinf^2*CL*S;
Di = 0.5*1.225*Vinf^2*CDi*S;

%% ADD VISCOUS EFFECTS
Dv = visc_drag(c,y0,Vinf,alfa,airfoil_data);
Mfoc = moment_f(c,y0,Vinf,alfa,airfoil_data);
Dtot = (Di + Dv)*1.1;
CD = Dtot/(0.5*1.225*Vinf^2*S);
E = L/Dtot;

end
