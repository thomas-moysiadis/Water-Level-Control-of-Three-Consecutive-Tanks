function [J] = obj_fun(DU,Phi,F,Rw,Y_dotsp,xtot_dot)

J = (Y_dotsp-F*xtot_dot)'*(Y_dotsp-F*xtot_dot) - 2*DU'*Phi'*(Y_dotsp-F*xtot_dot) + DU'*(Phi'*Phi+Rw)*DU;

end
