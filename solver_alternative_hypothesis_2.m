function F=solver_alternative_hypothesis_2(x,psi_xylem,rh,a1,a2,r3,L_ad,L_ab,gs_ad,gs_ab,p_sat,v_bar,R,T,k_vap,k_cw,p_atm,A_c_sym,A_c_vap,A_c_cw)

F(1) = x(1) + x(3) + x(4) - psi_xylem;
F(2) = x(2) + x(3)*r3*exp(r3*L_ad) - x(4)*r3*exp(-r3*L_ad);
F(3) = x(1) - x(5) - a2/a1*(x(3) + x(4) -x(7) - x(8));
F(4) = (k_cw*A_c_cw + k_vap*A_c_vap)*(1 + a2/a1)*r3*(x(4)*exp(-r3*L_ad) - x(3)*exp(r3*L_ad)) + (1e6)*gs_ad*(A_c_sym + A_c_vap + A_c_cw)*p_sat*(1 - rh +v_bar*(psi_xylem - x(3) - x(4) + r3*L_ad*(x(4)*exp(-r3*L_ad) - x(3)*exp(r3*L_ad)) - a2*x(3)*exp(r3*L_ad)/a1 - a2*x(4)*exp(-r3*L_ad)/a1)/(R*T))/p_atm;
F(5) = x(5) + x(7) + x(8) - psi_xylem;
F(6) = x(6) + x(7)*r3*exp(-r3*L_ab) - x(8)*r3*exp(r3*L_ab);
F(7) = x(6) - x(2) - a2/a1*r3*(x(7) - x(8) - x(3) + x(4));
F(8) = (k_cw*A_c_cw + k_vap*A_c_vap)*(1 + a2/a1)*r3*(x(8)*exp(r3*L_ab) - x(7)*exp(-r3*L_ab)) - (1e6)*gs_ab*(A_c_sym + A_c_vap + A_c_cw)*p_sat*(1 - rh +v_bar*(psi_xylem - x(7) - x(8) - r3*L_ab*(x(8)*exp(r3*L_ab) - x(7)*exp(-r3*L_ab)) - a2*x(7)*exp(-r3*L_ab)/a1 - a2*x(8)*exp(r3*L_ab)/a1)/(R*T))/p_atm;

end