function [dHdt] = ode_fun(t,y,U,wl1,wl2,wl3,wl4)

H1max = 0.35; %m
a = 0.25; %m
w = 0.035; %m
H2max = 0.35; %m
b = 0.345; %m
c = 0.1; %m
H3max = 0.35; %m
R = 0.364; %m

dH1dt = U(1)/(a*w) - U(2)*y(1)^(1/2)/(a*w) - wl1*y(1)^(1/2)/(a*w);
dH2dt = U(2)*y(1)^(1/2)/(c*w+y(2)/H2max*b*w) - U(3)*y(2)^(1/2)/(c*w+y(2)/H2max*b*w) + wl1*y(1)^(1/2)/(c*w+y(2)/H2max*b*w) - wl4*y(2)^(1/2)/(c*w+y(2)/H2max*b*w);
dH3dt = U(3)*y(2)^(1/2)/(w*(R^2-(R-y(3))^2)^(1/2)) - U(4)*y(3)^(1/2)/(w*(R^2-(R-y(3))^2)^(1/2)) + wl2*y(2)^(1/2)/(w*(R^2-(R-y(3))^2)^(1/2)) - wl3*y(3)^(1/2)/(w*(R^2-(R-y(3))^2)^(1/2));
dHdt = [dH1dt; dH2dt; dH3dt];

end
