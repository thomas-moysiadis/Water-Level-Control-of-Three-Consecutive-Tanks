clear
clc

%Define Constants
H1sp = 0.12; %m
H2sp = 0.1; %m
H3sp = 0.1; %m
H10 = 0.1; %m
H20 = 0.08; %m
H30 = 0.12; %m
H1max = 0.35; %m
a = 0.25; %m
w = 0.035; %m
H2max = 0.35; %m
b = 0.345; %m
c = 0.1; %m
H3max = 0.35; %m
R = 0.364; %m
wl1 = 0;
wl2 = 0;
wl3 = 0;
wl4 = 0;
qmax = 1.4*10^(-4); %m^3/s
Q1max = 5.5*10^(-5); %m^3/s
Q2max = 5.5*10^(-5); %m^3/s
Q3max = 5.5*10^(-5); %m^3/s
T = 0.01; %s
Np = 50;
Nc = 30;
Qw = eye(3*Np);
Rw = 70000*eye(4*Nc);
tf = 50; %s

%Equilibrium Points
H1s = H1sp;
H2s = H2sp;
H3s = H3sp;
qs = 2.4*10^(-5);
C1s = qs/H1s^(1/2) - wl1;
C2s = qs/H2s^(1/2) - wl4;
C3s = (qs - wl4*H2s^(1/2) + wl2*H2s^(1/2) - wl3*H3s^(1/2)) / H3s^(1/2);

%Linear Model
Ac = [-(C1s+wl1)/(2*a*w*H1s^(1/2)) 0 0; (C1s+wl1)/(2*H1s^(1/2)*(c*w+H2s/H2max*b*w)) -((C1s+wl1)*H1s^(1/2)*b*w)/(H2max*(c*w+H2s/H2max*b*w)^2)-((C2s+wl4)*(c*w+H2s/H2max*b*w)/(2*H2s^(1/2))-(C2s+wl4)*H2s^(1/2)*b*w/H2max)/(c*w+H2s/H2max*b*w)^2 0; 0 (C2s+wl2)/(2*w*H2s^(1/2)*(R^2-(R-H3s)^2)^(1/2)) (C2s+wl2)*H2s^(1/2)*(R-H3s)/(w*(R^2-(R-H3s)^2)^(3/2))-((C3s+wl3)*w*(R^2-(R-H3s)^2)^(1/2)/(2*(H3s)^(1/2))+(C3s+wl3)*H3s^(1/2)*w*(R-H3s)/(R^2-(R-H3s)^2)^(1/2))/(w^2*((R^2-(R-H3s)^2)^(1/2))^2)];
Bc = [1/(a*w) -(H1s)^(1/2)/(a*w) 0 0; 0 H1s^(1/2)/(c*w+H2s/H2max*b*w) -(H2s)^(1/2)/(c*w+H2s/H2max*b*w) 0; 0 0 H2s^(1/2)/(w*(R^2-(R-H3s)^2)^(1/2)) -(H3s)^(1/2)/(w*(R^2-(R-H3s)^2)^(1/2))];
Cc = eye(3);
Dc = zeros(3,4);

%Discretization
sysc = ss(Ac,Bc,Cc,Dc);
sysd = c2d(sysc,T);
[Ad, Bd, Cd, Dd] = ssdata(sysd);

%Total Model
A = [Ad zeros(3,3); Cd*Ad eye(3,3)];
B = [Bd; Cd*Bd];
C = [zeros(3,3) eye(3,3)];
D = zeros(3,4);

%Total Matrices
F = zeros(3*Np,6);
for i = 1:Np
    F(3*i-2:3*i,:) = C*A^i;
end

Phi = zeros(3*Np,4*Nc);
for j = 1:Nc
    for i = j:Np
        Phi(3*i-2:3*i,4*j-3:4*j) = C*A^(i-j)*B;
    end
end

%Compute Optimal Control Action
Y_dotsp = zeros(3*Np,1);

Qmax = [Q1max; Q2max; Q3max];
Usp = [qs; C1s; C2s; C3s];
xsp = [H1sp; H2sp; H3sp];
x0 = [H10; H20; H30];
U_h = zeros(4*Nc,1);
x_dot = x0 - xsp;
dx_dot = zeros(3,1);
xtot_dot = [dx_dot; x_dot];
y0 = x0;
t1 = 0:T:tf;
Utot = zeros(length(t1)-1,4);
DUtot = zeros(length(t1)-1,4);
xtot = zeros(length(t1),3);
xtot(1,:) = x0';

for i = 1:length(t1)-1
    %Optimal Control Action
    DU_dot = inv(Phi'*Phi + Rw)*Phi'*(Y_dotsp - F*xtot_dot);
    for k = 1:Nc
        if k == 1
            U_h(1:4,1) = U_h(1:4,1) + DU_dot(1:4,1);
        else
            U_h(4*k-3:4*k,1) = U_h(4*(k-1)-3:4*(k-1),1) + DU_dot(4*k-3:4*k,1);
        end
    end
    U = U_h(1:4,1);
    
    %System Constraints
    if U(1,1) > qmax
        U(1,1) = qmax;
    elseif U(1,1) < 0
        U(1,1) = 0;
    end
    
    for j = 1:3
        if U(j+1,1)*sqrt(xtot_dot(j+3,1)+xsp(j,1)) > Qmax(j,1)
            U(j+1,1) = Qmax(j,1) / sqrt(xtot_dot(j+3,1)+xsp(j,1));
        elseif U(j+1,1) < 0
            U(j+1,1) = 0;
        end
    end
    
    Utot(i,1) = U(1,1);
    Utot(i,2:4) = (U(2:4,1) .* sqrt((xtot_dot(4:6,1) + xsp)))';
    DUtot(i,1) = DU_dot(1,1);
    if i == 1
        DUtot(i,2:4) = Utot(i,2:4);
    else
        DUtot(i,2:4) = Utot(i,2:4) - Utot(i-1,2:4);
    end
    
    %Disturbances
    if i >= 300 && i <= 800
        wl1 = 5*10^(-5);
        wl2 = 0;
        wl3 = 0;
        wl4 = 0;
    elseif i >= 1500 && i <= 2000
        wl1 = 0;
        wl2 = 10^(-4);
        wl3 = 0;
        wl4 = 10^(-4);
    elseif i >= 3000 && i<= 3500
        wl1 = 0;
        wl2 = 0;
        wl3 = 10^(-5);
        wl4 = 0;
    else
        wl1 = 0;
        wl2 = 0;
        wl3 = 0;
        wl4 = 0;
    end
    %System Feedback
    [t,y] = ode45(@(t,y) ode_fun(t,y,U,wl1,wl2,wl3,wl4), t1(i):T:t1(i)+T, y0);
    y0 = y(end,:)';
    
    x = y0;
    x_dot = x - xsp;
    dx_dot = x_dot - xtot_dot(4:6,1);
    xtot_dot = [dx_dot; x_dot];
    xtot(i+1,:) = x';
    
end

%Plots
f1 = figure;
hold on
stairs(t1(1:i), Utot(:,1), "r");
stairs(t1(1:i), Utot(:,2), "b");
stairs(t1(1:i), Utot(:,3), "m");
stairs(t1(1:i), Utot(:,4), "k");
legend(["q" "Q1" "Q2" "Q3"]);
xlabel("Time (s)");
ylabel("Input (m^3/s)");
title("MPC inputs");
grid on
hold off

f2 = figure;
hold on
stairs(t1(1:i), DUtot(:,1), "r");
stairs(t1(1:i), DUtot(:,2), "b");
stairs(t1(1:i), DUtot(:,3), "m");
stairs(t1(1:i), DUtot(:,4), "k");
legend(["Dq" "DQ1" "DQ2" "DQ3"]);
xlabel("Time (s)");
ylabel("Input increment (m^3/s)");
title("MPC input increments");
grid on
hold off

f3 = figure;
hold on
plot(t1, xtot(:,1), "r");
plot(t1, xtot(:,2), "b");
plot(t1, xtot(:,3), "m");
legend(["H1" "H2" "H3"]);
xlabel("Time (s)");
ylabel("Output (m)");
title("MPC outputs");
grid on
hold off
