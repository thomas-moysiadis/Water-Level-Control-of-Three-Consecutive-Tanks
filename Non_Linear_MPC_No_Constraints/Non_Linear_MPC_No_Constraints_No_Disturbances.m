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
Np = 15;
Nc = 8;
Qw = 1;
Rw = 5;
tf = 10; %s

%Equilibrium Points
H1s = H1sp;
H2s = H2sp;
H3s = H3sp;
qs = 2.4*10^(-5);
C1s = qs/H1s^(1/2) - wl1;
C2s = qs/H2s^(1/2) - wl4;
C3s = (qs - wl4*H2s^(1/2) + wl2*H2s^(1/2) - wl3*H3s^(1/2)) / H3s^(1/2);
Qsp = [qs; C1s; C2s; C3s]';

%Non Linear MPC Construction
Qmax = [Q1max; Q2max; Q3max];
xsp = [H1sp; H2sp; H3sp]';
x0 = [H10; H20; H30];
xk = x0;
U0 = zeros(4,1);
y0 = x0;
t1 = 0:T:tf;
Utot = zeros(length(t1)-1,4);
DUtot = zeros(length(t1)-1,4);
xtot = zeros(length(t1),3);
xtot(1,:) = x0';

for i = 1:length(t1)-1
    nlobj = nlmpc(3,3,4);
    
    nlobj.Ts = T;
    nlobj.PredictionHorizon = Np;
    nlobj.ControlHorizon = Nc;
    
    nlobj.Model.StateFcn = @myStateFunction;
    nlobj.Jacobian.StateFcn = @myStateJacobian;
    nlobj.Model.IsContinuousTime = true;

    nlobj.Weights.ManipulatedVariables = [Rw Rw Rw Rw];
    nlobj.Weights.ManipulatedVariablesRate = [Qw Qw Qw Qw];
    nlobj.Optimization.UseSuboptimalSolution = true;
    
    options = nlmpcmoveopt;
    options.MVTarget = Qsp;
    %validateFcns(nlobj,xk,U0,[],[],xsp,Qsp);
    %Optimal Control Action
    U = nlmpcmove(nlobj,xk,U0,xsp,[],options);
    
    %System Constraints
    if U(1) > qmax
        U(1) = qmax;
    elseif U(1) < 0
        U(1) = 0;
    end
    
    for j = 1:3
        if U(j+1) > Qmax(j,1) / sqrt(x0(j))
            U(j+1) = Qmax(j,1) / sqrt(x0(j));
        elseif U(j+1) < 0
            U(j+1) = 0;
        end
    end
    
    Utot(i,1) = U(1);
    Utot(i,2:4) = U(2:4,1)'.*sqrt(x0');
    if i == 1
        DUtot(i,:) = Utot(i,:);
    else
        DUtot(i,:) = Utot(i,:) - Utot(i-1,:);
    end
    
    %System Feedback
    [t,y] = ode45(@(t,y) ode_fun(t,y,U,wl1,wl2,wl3,wl4), t1(i):T:t1(i)+T, y0);
    y0 = y(end,:)';
    
    xk = y0;
    xtot(i+1,:) = xk';
    U0 = U;
    
    disp(i)
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
