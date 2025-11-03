function F=solver_separate_bs_and_mesophyll(x,psi_xylem,rh,a1,a2,b1,b2,r3,s3,L_ad,L_ab,L_bs,gs_ad,gs_ab,p_sat,v_bar,R,T,k_mes_vap,k_bs_vap,k_mes_cw,k_bs_cw,p_atm,A_sym,A_vap,A_cw)

F(1) = x(1) + x(3) + x(4) - psi_xylem; % psi_cc = psi_xyl at z=0+

F(2) = x(5) + x(7) + x(8) - psi_xylem; % psi_cc = psi_xyl at z=0-

F(3) = x(9) + x(10)*L_bs + x(11)*exp(s3*L_bs) + x(12)*exp(-s3*L_bs) - x(1) - x(2)*L_bs - x(3)*exp(r3*L_bs) - x(4)*exp(-r3*L_bs); % equality of psi_cc at adaxial bundle sheath-mesophyll interface

F(4) = x(13) - x(14)*L_bs + x(15)*exp(-s3*L_bs) + x(16)*exp(s3*L_bs) - x(5) + x(6)*L_bs - x(7)*exp(-r3*L_bs) - x(8)*exp(r3*L_bs); % equality of psi_cc at abaxial bundle sheath-mesophyll interface

F(5) = x(1) - x(5) - a2/a1*(x(3) + x(4) -x(7) - x(8)); % equality of psi_apo at adaxial-abaxial interface

F(6) = x(6) - x(2) - a2/a1*r3*(x(7) - x(8) - x(3) + x(4)); % equality of apoplasmic flux at adaxial-abaxial interface

F(7) = x(2) + x(3)*r3*exp(r3*L_bs) - x(4)*r3*exp(-r3*L_bs) - x(10) - x(11)*s3*exp(s3*L_bs) + x(12)*s3*exp(-s3*L_bs); % equality of cell-to-cell flux at adaxial bundle sheath-mesophyll interface

F(8) = x(6) + x(7)*r3*exp(-r3*L_bs) - x(8)*r3*exp(r3*L_bs) - x(14) - x(15)*s3*exp(-s3*L_bs) + x(16)*s3*exp(s3*L_bs); % equality of cell-to-cell flux at abaxial bundle sheath-mesophyll interface

F(9) = x(10) + x(11)*s3*exp(s3*L_ad) - x(12)*s3*exp(-s3*L_ad); % no flux b.c. on cell-to-cell path at z = L_ad

F(10) = x(14) + x(15)*s3*exp(-s3*L_ab) - x(16)*s3*exp(s3*L_ab); % no flux b.c. on cell-to-cell path at z = -L_ab

F(11) = (k_bs_cw*A_cw + k_bs_vap*A_vap)*(x(2) - a2/a1*x(3)*r3*exp(r3*L_bs) + a2/a1*x(4)*r3*exp(-r3*L_bs)) - (k_mes_cw*A_cw + k_mes_vap*A_vap)*(x(10) - b2/b1*x(11)*s3*exp(s3*L_bs) + b2/b1*x(12)*s3*exp(-s3*L_bs)); % equality of apoplasmic flux at adaxial bundle sheath-mesophyll interface

F(12) = x(9) + x(10)*L_bs - b2/b1*x(11)*exp(s3*L_bs) - b2/b1*x(12)*exp(-s3*L_bs) - x(1) - x(2)*L_bs + a2/a1*x(3)*exp(r3*L_bs) + a2/a1*x(4)*exp(-r3*L_bs); % equality of psi_apo at adaxial bundle sheath-mesophyll interface

F(13) = x(5) - x(6)*L_bs - a2/a1*x(7)*exp(-r3*L_bs) - a2/a1*x(8)*exp(r3*L_bs) - x(13) + x(14)*L_bs + b2/b1*x(15)*exp(-s3*L_bs) + b2/b1*x(16)*exp(s3*L_bs); % equality of psi_cc at abaxial bundle sheath-mesophyll interface

F(14) = (k_mes_cw*A_cw + k_mes_vap*A_vap)*(x(14)  - b2/b1*x(15)*s3*exp(-s3*L_bs) + b2/b1*x(16)*s3*exp(s3*L_bs)) - (k_bs_cw*A_cw + k_bs_vap*A_vap)*(x(6) - a2/a1*x(7)*r3*exp(-r3*L_bs) + a2/a1*x(8)*r3*exp(r3*L_bs)); % equality of apoplasmic flux at abaxial bundle sheath-mesophyll interface

F(15) = (k_mes_cw*A_cw + k_mes_vap*A_vap)*(1 + b2/b1)*s3*(x(12)*exp(-s3*L_ad) - x(11)*exp(s3*L_ad)) + (1e6)*gs_ad*(A_sym + A_vap + A_cw)*p_sat*(1 - rh +v_bar*(x(9) + x(10)*L_ad - b2*x(11)*exp(s3*L_ad)/b1 - b2*x(12)*exp(-s3*L_ad)/b1)/(R*T))/p_atm; % flux b.c. at adaxial stomates

F(16) = (k_mes_cw*A_cw + k_mes_vap*A_vap)*(1 + b2/b1)*s3*(x(16)*exp(s3*L_ab) - x(15)*exp(-s3*L_ab)) - (1e6)*gs_ab*(A_sym + A_vap + A_cw)*p_sat*(1 - rh +v_bar*(x(13) - x(14)*L_ab - b2*x(15)*exp(-s3*L_ab)/b1 - b2*x(16)*exp(s3*L_ab)/b1)/(R*T))/p_atm; % flux b.c. at abaxial stomates

end