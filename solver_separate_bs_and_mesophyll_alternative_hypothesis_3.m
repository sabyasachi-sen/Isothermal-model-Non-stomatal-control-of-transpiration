function F=solver_separate_bs_and_mesophyll_alternative_hypothesis_3(x,psi_xylem,rh,b1,b2,r1,r2,r3,r4,c1,c2,c3,c4,s3,L_ad,L_ab,L_bs,gs_ad,gs_ab,p_sat,v_bar,R_gas,T,k_mes_vap,k_bs_vap,k_mes_cw,k_bs_cw,p_atm,A_sym,A_vap,A_cw)

F(1) = x(1) + x(2) + x(3) + x(4) - (x(5) + x(6) + x(7) + x(8)); % equality of psi_cc at z=0+

F(2) = r1*x(1) + r2*x(2) + r3*x(3) + r4*x(4) - (r1*x(5) + r2*x(6) + r3*x(7) + r4*x(8)); % equality of psi_cc at z=0-

F(3) = x(9) + x(10)*L_bs + x(11)*exp(s3*L_bs) + x(12)*exp(-s3*L_bs) - (psi_xylem + x(1)*exp(r1*L_bs) + x(2)*exp(r2*L_bs) + x(3)*exp(r3*L_bs) + x(4)*exp(r4*L_bs)); % equality of psi_cc at adaxial bundle sheath-mesophyll interface

F(4) = x(13) - x(14)*L_bs + x(15)*exp(-s3*L_bs) + x(16)*exp(s3*L_bs) - (psi_xylem + x(5)*exp(-r1*L_bs) + x(6)*exp(-r2*L_bs) + x(7)*exp(-r3*L_bs) + x(8)*exp(-r4*L_bs)); % equality of psi_cc at abaxial bundle sheath-mesophyll interface

F(5) = c1*x(1) + c2*x(2) + c3*x(3) + c4*x(4) - (c1*x(5) + c2*x(6) + c3*x(7) + c4*x(8)); % equality of psi_apo at adaxial-abaxial interface

F(6) = c1*x(5)*r1 + c2*x(6)*r2 + c3*x(7)*r3 + c4*x(8)*r4 - (c1*x(1)*r1 + c2*x(2)*r2 + c3*x(3)*r3 + c4*x(4)*r4); % equality of apoplasmic flux at adaxial-abaxial interface

F(7) = x(1)*r1*exp(r1*L_bs) + x(2)*r2*exp(r2*L_bs) + x(3)*r3*exp(r3*L_bs) + x(4)*r4*exp(r4*L_bs) - (x(10) + x(11)*s3*exp(s3*L_bs) - x(12)*s3*exp(-s3*L_bs)); % equality of cell-to-cell flux at adaxial bundle sheath-mesophyll interface

F(8) = x(5)*r1*exp(-r1*L_bs) + x(6)*r2*exp(-r2*L_bs) + x(7)*r3*exp(-r3*L_bs) + x(8)*r4*exp(-r4*L_bs) - (x(14) + x(15)*s3*exp(-s3*L_bs) - x(16)*s3*exp(s3*L_bs)); % equality of cell-to-cell flux at abaxial bundle sheath-mesophyll interface

F(9) = x(10) + x(11)*s3*exp(s3*L_ad) - x(12)*s3*exp(-s3*L_ad); % no flux b.c. on cell-to-cell path at z = L_ad

F(10) = x(14) + x(15)*s3*exp(-s3*L_ab) - x(16)*s3*exp(s3*L_ab); % no flux b.c. on cell-to-cell path at z = -L_ab

F(11) = (k_bs_cw*A_cw + k_bs_vap*A_vap)*(c1*x(1)*r1*exp(r1*L_bs) + c2*x(2)*r2*exp(r2*L_bs) + c3*x(3)*r3*exp(r3*L_bs) + c4*x(4)*r4*exp(r4*L_bs)) - (k_mes_cw*A_cw + k_mes_vap*A_vap)*(x(10) - b2/b1*x(11)*s3*exp(s3*L_bs) + b2/b1*x(12)*s3*exp(-s3*L_bs)); % equality of apoplasmic flux at adaxial bundle sheath-mesophyll interface

F(12) = x(9) + x(10)*L_bs - b2/b1*x(11)*exp(s3*L_bs) - b2/b1*x(12)*exp(-s3*L_bs) - (psi_xylem + c1*x(1)*exp(r1*L_bs) + c2*x(2)*exp(r2*L_bs) + c3*x(3)*exp(r3*L_bs) + c4*x(4)*exp(r4*L_bs)); % equality of psi_apo at adaxial bundle sheath-mesophyll interface

F(13) = psi_xylem + c1*x(5)*exp(-r1*L_bs) + c2*x(6)*exp(-r2*L_bs) + c3*x(7)*exp(-r3*L_bs) + c4*x(8)*exp(-r4*L_bs) - (x(13) - x(14)*L_bs - b2/b1*x(15)*exp(-s3*L_bs) - b2/b1*x(16)*exp(s3*L_bs)); % equality of psi_cc at abaxial bundle sheath-mesophyll interface

F(14) = (k_mes_cw*A_cw + k_mes_vap*A_vap)*(x(14)  - b2/b1*x(15)*s3*exp(-s3*L_bs) + b2/b1*x(16)*s3*exp(s3*L_bs)) - (k_bs_cw*A_cw + k_bs_vap*A_vap)*(c1*x(5)*r1*exp(-r1*L_bs) + c2*x(6)*r2*exp(-r2*L_bs) + c3*x(7)*r3*exp(-r3*L_bs) + c4*x(8)*r4*exp(-r4*L_bs)); % equality of apoplasmic flux at abaxial bundle sheath-mesophyll interface

F(15) = (k_mes_cw*A_cw + k_mes_vap*A_vap)*(1 + b2/b1)*s3*(x(12)*exp(-s3*L_ad) - x(11)*exp(s3*L_ad)) + (1e6)*gs_ad*(A_sym + A_vap + A_cw)*p_sat*(1 - rh +v_bar*(x(9) + x(10)*L_ad - b2*x(11)*exp(s3*L_ad)/b1 - b2*x(12)*exp(-s3*L_ad)/b1)/(R_gas*T))/p_atm; % flux b.c. at adaxial stomates

F(16) = (k_mes_cw*A_cw + k_mes_vap*A_vap)*(1 + b2/b1)*s3*(x(16)*exp(s3*L_ab) - x(15)*exp(-s3*L_ab)) - (1e6)*gs_ab*(A_sym + A_vap + A_cw)*p_sat*(1 - rh +v_bar*(x(13) - x(14)*L_ab - b2*x(15)*exp(-s3*L_ab)/b1 - b2*x(16)*exp(s3*L_ab)/b1)/(R_gas*T))/p_atm; % flux b.c. at abaxial stomates

end