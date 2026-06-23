function [dHdt] = myStateFunction(X,U)

H1max = 0.35; %m
a = 0.25; %m
w = 0.035; %m
H2max = 0.35; %m
b = 0.345; %m
c = 0.1; %m
H3max = 0.35; %m
R = 0.364; %m
params = zeros(4,1);

dH1dt = U(1)/(a*w) - U(2)*X(1)^(1/2)/(a*w) - params(1)*X(1)^(1/2)/(a*w);
dH2dt = U(2)*X(1)^(1/2)/(c*w+X(2)/H2max*b*w) - U(3)*X(2)^(1/2)/(c*w+X(2)/H2max*b*w) + params(1)*X(1)^(1/2)/(c*w+X(2)/H2max*b*w) - params(4)*X(2)^(1/2)/(c*w+X(2)/H2max*b*w);
dH3dt = U(3)*X(2)^(1/2)/(w*(R^2-(R-X(3))^2)^(1/2)) - U(4)*X(3)^(1/2)/(w*(R^2-(R-X(3))^2)^(1/2)) + params(2)*X(2)^(1/2)/(w*(R^2-(R-X(3))^2)^(1/2)) - params(3)*X(3)^(1/2)/(w*(R^2-(R-X(3))^2)^(1/2));
dHdt = [dH1dt; dH2dt; dH3dt];

end
