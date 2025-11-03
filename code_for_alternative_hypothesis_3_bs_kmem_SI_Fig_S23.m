% code_for_alternative_architecture_bs_kmem_SI_Fig_S22.m
% Description: This code executes the mathematical model developed in SI
% Section S10 in the paper titled, "Loss of plasma membrane conductance of 
% outside-xylem zone explains non-stomatal control of transpiration" for 
% alternative architecture 3: bs membranes define the xylem-bs and 
% bs-mesophyll interface and undergoes loss of conductance under drought.
% We generate 3 predictions (mean and bounds of uncertainty interval) with 
% npoints (see below) being the number of points in the discretized domain 
% for each prediction. npoints are distributed across four ranges in 
% psi_xyl: -0.4 to 0; -0.8 to -0.4; -1.2 to -0.8 and -1.6 to -1.2 (MPa). 
% The prediction mean and bounds of uncertainty interval is then compared to 
% the data mean +/- S.E. 
% This code uses the following MATLAB toolboxes:
% 1. Optimization Toolbox
% 2. Symbolic Math Toolbox
% Code developed by: Sabyasachi Sen
% Email: ss3945@cornell.edu
% Affiliation: Cornell University
% Created: February 8, 2025
% Last Modified: February 10, 2025
% Copyright (c) 2025 Sabyasachi Sen
% If you use this code in your research, please cite as:
% "Sabyasachi Sen. (2025). 
% code_for_alternative_architecture_bs_kmem_SI_Fig_S22.m
% Retrieved from manuscript titled: Loss of plasma membrane conductance of 
% outside-xylem zone explains non-stomatal control of transpiration"
% This code is provided under the MIT License.
% For full license terms, see: https://opensource.org/licenses/MIT

clear;
close all;
clc;

syms var;

%% Thermodynamic parameters (see Table S3 in SI Section S10)

leaf_temp=309; % leaf temperature (=36C) in Kelvin
p_sat=5946.178; % saturation pressure in Pa at 36C
R=8.314; % universal gas constant in SI units Pa.m3/K/mol
v_bar=18e-6; % molar volume of water in m3/mol
c=39.45; % molar concentration of air in mol/m3
D=2.5e-5; % diffusivity of water vapor in air in m2/s
p_atm=101325; % atmospheric pressure in Pa
kai_0=p_sat/p_atm; % reference mole fraction of water vapor at leaf temperature

%% Introduce uncertainty in parameters to generate uncertainty interval in predictions (see Tables S4 and S5 in SI Section S10)

% Notation: variable name: high, mid, low values
A_cc_arr=[5e-9, 6e-9, 7e-9]; % x-sectional area for flow in cell-to-cell path (in m2)
A_vap_arr=[5e-9, 5.5e-9, 6e-9]; % x-sectional area for flow in apoplasmic path (in m2)
bs_tortuosity_arr=[1, 2, 4]; % bundle sheath tortuosity (dimensionless)
mes_tortuosity_arr=[1, 1.5, 4]; % mesophyll tortuosity (dimensionless)
kappa_bs_mem_max_arr=[10e-6, 9e-6, 8e-6]; % maximum membrane conductivity in mol/m/s/MPa
kappa_bs_mem_min_arr=[3e-7, 2e-7, 1e-7]; % minimum membrane conductivity in mol/m/s/MPa
kappa_bs_mem_trigger_arr=[-0.24e6, -0.23e6, -0.22e6]; % xylem stress at which 'characteristic' loss of membrane conductivity occurs (in Pa)
kappa_bs_mem_decline_sensitivity_arr=[1e-5, 1.1e-5, 1.2e-5]; % determines sensitivity of membrane conductivity to xylem stress in MPa^(-1)

%% Constant leaf anatomical parameters (see Table S4 in SI Section S10)

L_total=2e-4; % total thickness of leaf in m
L_ad=L_total/2; % Assume half-thickness of leaf 
L_ab=L_total/2; % Assume half-thickness of leaf 
L_bs=L_ad/3; % Assume bundle sheath occupies 1/3rd of leaf thickness
A_cw=0.08e-9; % x-sectional area for flow in cell wall (in m2)

%% Constant transport parameters (see Table S5 in SI Section S10)

kappa_mes_mem=8e-6; % constant, high conductivity mesophyll membranes
gs_ad=90e-3; % adaxial stomatal conductance in mol/m2/s
gs_ab=90e-3; % abaxial stomatal conductance in mol/m2/s
gs=(gs_ad+gs_ab)*1e3; % total stomatal conductance in mmol/m2/s: will be used to calculate transpiration rate

%% Driving forces for transpiration

% psi_xylem sweep
npoints=30; % number of nodes in the discretized domain
psi_xylem=flipud(linspace(-1.6,-0.1,npoints)')*1e6; % xylem water potential in Pa
xylem_means=[-0.2;-0.6;-1;-1.4]; % for plotting in bins

% VPD sweep 
rel_hum = flipud(linspace(0.36,0.7,npoints)'); % relative humidity of bulk air outside the leaf
vpd = p_sat*(1-rel_hum)/1e3; % vpd in kPa: the mean vpd in our measurements is ~3.8 kPa which is the last element in the vpd array here

%% Finite difference discretization parameters:

N=1000; % Total # of nodes including boundary nodes
h=L_total/(N-1); % interval between nodes

%% Coordinates along leaf thickness direction (in metres)

coords_ad=linspace(0,L_ad,N)'; % adaxial coordinates
coords_ab=linspace(-L_ab,0,N)'; % abaxial coordinates
coords_bs_ad=coords_ad(1:round(N/3)); % adaxial bundle sheath coordinates
coords_bs_ab=coords_ab(round(2*N/3)+1:end); % abaxial bundle sheath coordinates
coords_mes_ad=coords_ad(round(N/3):end); % adaxial mesophyll coordinates
coords_mes_ab=coords_ab(1:round(2*N/3)+1); % abaxial mesophyll coordinates
coordinates_ad=coords_ad*1e6; % coordinates for adaxial section: microns for plotting
coordinates_ab=coords_ab*1e6; % coordinates for abaxial section: microns for plotting

%% Loop through 3 sets of parameters

for ii=1:3

    % Variable leaf anatomical parameters and transport parameters (see Tables S4 and S5 in SI Section S10)

    bs_tortuosity=bs_tortuosity_arr(ii); % bundle sheath tortuosity
    mes_tortuosity=mes_tortuosity_arr(ii); % mesophyll tortuosity

    k_bs_vap=c*D*v_bar*kai_0/(R*leaf_temp)/bs_tortuosity/1e-6; % hydraulic conductivity of vapor in bundle sheath in mol/m/s/MPa; this assumes driving force for transport is water potential. Ref: Rockwell 2014
    k_mes_vap=c*D*v_bar*kai_0/(R*leaf_temp)/mes_tortuosity/1e-6; % hydraulic conductivity of vapor in mesophyll in mol/m/s/MPa; this assumes driving force for transport is water potential. Ref: Rockwell 2014
    
    k_bs_cw=2.7e-6/bs_tortuosity; % hydraulic conductivity of bundle sheath cell wall in mol/m/s/MPa; this assumes driving force for transport is water potential gradient. Refs: Michael potato parenchyma 1997; Rockwell et al. 2022
    k_mes_cw=2.7e-6/mes_tortuosity; % hydraulic conductivity of mesophyll cell wall in mol/m/s/MPa; this assumes driving force for transport is water potential gradient. Refs: Michael potato parenchyma 1997; Rockwell et al. 2022

    A_cc=A_cc_arr(ii); % see SI Section S10 for choice of A_cc
    A_vap=A_vap_arr(ii); % see SI Section S10 for choice of A_vap
                
    %% Tissue hydraulic conductivity parameters 
    
    % Membrane hydraulic parameters for current iteration
    kappa_bs_mem_max=kappa_bs_mem_max_arr(ii); 
    kappa_bs_mem_decline_sensitivity=kappa_bs_mem_decline_sensitivity_arr(ii); 
    kappa_bs_mem_trigger=kappa_bs_mem_trigger_arr(ii); 
    kappa_bs_mem_min=kappa_bs_mem_min_arr(ii); 
    % sigmoid function
    kappa_bs_mem_arr=kappa_bs_mem_max*1./(1+exp(-kappa_bs_mem_decline_sensitivity*(psi_xylem-kappa_bs_mem_trigger)))+kappa_bs_mem_min; 
    
    % plot kappa_mem vs psi_xyl
    figure(1);hold on;
    plot(psi_xylem/1e6,kappa_bs_mem_arr,'k-','linewidth',2);hold on;
    pbaspect([1 1 1]);
    ylabel('\kappa_{mem}^{bs} (mol/m/s/MPa)');ylim([1e-8 2e-5]);
    xlabel('\psi_{xyl} (MPa)');xlim([-1.7 0]);xticks(flipud(xylem_means));
    set(gca,'fontsize',24,'ticklabelinterpreter','tex','fontname','Arial','fontweight','bold','yscale','log');
    box on;ax=gca;ax.LineWidth = 2;
    
    %% Solution
    
    % Initialize arrays that will contain the solutions: psi_cc and psi_apo (in Pa)
    psi_bs_cc_ad_arr=zeros(length(psi_xylem),length(rel_hum)); % psi_cc in bundle sheath adaxial domain
    psi_bs_cc_ab_arr=zeros(length(psi_xylem),length(rel_hum)); % psi_cc in bundle sheath abaxial domain
    psi_bs_apo_ad_arr=zeros(length(psi_xylem),length(rel_hum)); % psi_apo in bundle sheath adaxial domain
    psi_bs_apo_ab_arr=zeros(length(psi_xylem),length(rel_hum)); % psi_apo in bundle sheath abaxial domain

    psi_mes_cc_ad_N=zeros(length(psi_xylem),length(rel_hum)); % psi_cc in mesophyll adaxial domain
    psi_mes_cc_ab_N=zeros(length(psi_xylem),length(rel_hum)); % psi_cc in mesophyll abaxial domain
    psi_mes_apo_ad_N=zeros(length(psi_xylem),length(rel_hum)); % psi_apo in mesophyll adaxial domain
    psi_mes_apo_ab_N=zeros(length(psi_xylem),length(rel_hum)); % psi_apo in mesophyll abaxial domain
    
    for i=1:length(psi_xylem)
    
        kappa_bs_mem=kappa_bs_mem_arr(i); % mesophyll membrane conductivity 
        k_cc=kappa_bs_mem; % cell-to-cell conductivity 
        
        for j=1:length(rel_hum)
    
            % see SI Section S10 for details of these constants
            % a1, a2, b1, and b2 multiply the derivates in SI equations S30 and S31
    
            % bundle sheath hydraulic properties
            a1=kappa_bs_mem/(k_cc*A_cc); 
            a2=kappa_bs_mem/(k_bs_cw*A_cw + k_bs_vap*A_vap);
            r=double(solve(var^4 - (2*a1+a2)*var^2 + a1*a2 == 0, var));
            r1=r(1);
            r2=r(2);
            r3=r(3);
            r4=r(4);
    
            rat=(2*a1-r.^2)./a1;
            c1=rat(1);
            c2=rat(2);        
            c3=rat(3);
            c4=rat(4);
            
            % mesophyll hydraulic properties        
            b1=kappa_mes_mem/(k_cc*A_cc); 
            b2=kappa_mes_mem/(k_mes_cw*A_cw + k_mes_vap*A_vap);
            s3=+sqrt(b1+b2);

            init=[-1e6,-1e8,-10,-1e6,-2e6,-2e8,-20,-2e6,-1.5e6,-1.5e8,-25,-1.5e6,-2.5e6,-2.5e8,-30,-2.5e6];
            x=fsolve(@(x)solver_separate_bs_and_mesophyll_alternative_hypothesis_3(x,psi_xylem(i),rel_hum(j),b1,b2,r1,r2,r3,r4,c1,c2,c3,c4,s3,L_ad,L_ab,L_bs,gs_ad,gs_ab,p_sat,v_bar,R,leaf_temp,k_mes_vap,k_bs_vap,k_mes_cw,k_bs_cw,p_atm,A_cc,A_vap,A_cw),init);
            
            %  psi_bs
            psi_bs_cc_ad_arr = psi_xylem(i) + x(1)*exp(r1*coords_bs_ad) + x(2)*exp(r2*coords_bs_ad) + x(3)*exp(r3*coords_bs_ad) + x(4)*exp(r4*coords_bs_ad);
            psi_bs_apo_ad_arr = psi_xylem(i) + c1*x(1)*exp(r1*coords_bs_ad) + c2*x(2)*exp(r2*coords_bs_ad) + c3*x(3)*exp(r3*coords_bs_ad) + c4*x(4)*exp(r4*coords_bs_ad);
            psi_bs_cc_ab_arr = psi_xylem(i) + x(5)*exp(r1*coords_bs_ab) + x(6)*exp(r2*coords_bs_ab) + x(7)*exp(r3*coords_bs_ab) + x(8)*exp(r4*coords_bs_ab);
            psi_bs_apo_ab_arr = psi_xylem(i) + c1*x(5)*exp(r1*coords_bs_ab) + c2*x(6)*exp(r2*coords_bs_ab) + c3*x(7)*exp(r3*coords_bs_ab) + c4*x(8)*exp(r4*coords_bs_ab);
            
            % psi_mes
            psi_mes_cc_ad_arr = x(9) + x(10)*coords_mes_ad + x(11)*exp(s3*coords_mes_ad) + x(12)*exp(-s3*coords_mes_ad);
            psi_mes_apo_ad_arr = x(9) + x(10)*coords_mes_ad - b2*x(11)*exp(s3*coords_mes_ad)/b1 - b2*x(12)*exp(-s3*coords_mes_ad)/b1;
            psi_mes_cc_ab_arr = x(13) + x(14)*coords_mes_ab + x(15)*exp(s3*coords_mes_ab) + x(16)*exp(-s3*coords_mes_ab);
            psi_mes_apo_ab_arr = x(13) + x(14)*coords_mes_ab - b2*x(15)*exp(s3*coords_mes_ab)/b1 - b2*x(16)*exp(-s3*coords_mes_ab)/b1;        

            % Combining the two cell-to-cell domains
            psi_cc_ab_arr=[psi_mes_cc_ab_arr(1:end-1);psi_bs_cc_ab_arr];
            psi_cc_ad_arr=[psi_bs_cc_ad_arr(1:end-1);psi_mes_cc_ad_arr];
            
            % Combining the two apoplasm domains
            psi_apo_ab_arr=[psi_mes_apo_ab_arr(1:end-1);psi_bs_apo_ab_arr];
            psi_apo_ad_arr=[psi_bs_apo_ad_arr(1:end-1);psi_mes_apo_ad_arr];

            % Terminal apoplasm and symplasm potentials on adaxial side 
            psi_mes_cc_ad_N(i,j)=psi_mes_cc_ad_arr(end)/1e6; % divide by 1e6 to convert from Pa to MPa
            psi_mes_apo_ad_N(i,j)=psi_mes_apo_ad_arr(end)/1e6; % divide by 1e6 to convert from Pa to MPa 
            
            % Terminal apoplasm and symplasm potentials on abaxial side 
            psi_mes_cc_ab_N(i,j)=psi_mes_cc_ab_arr(end)/1e6; % divide by 1e6 to convert from Pa to MPa
            psi_mes_apo_ab_N(i,j)=psi_mes_apo_ab_arr(end)/1e6; % divide by 1e6 to convert from Pa to MPa
                    
        end
        
        kappa_bs_mem_exp=floor(log10(kappa_bs_mem));
        kappa_bs_mem_a=kappa_bs_mem/(10^kappa_bs_mem_exp);

        % Plot psi_apoplast vs z
        figure(2);
        plot([psi_cc_ab_arr; psi_cc_ad_arr]/1e6,[coordinates_ab; coordinates_ad],[psi_apo_ab_arr; psi_apo_ad_arr]/1e6,[coordinates_ab; coordinates_ad],'linewidth',4);
        pbaspect([1 1 1]);
        xlabel('\psi(MPa)','interpreter','tex');xlim([-20 0]);
        ylabel('z (\mu m)','interpreter','tex');ylim([-110 110]);
        set(gca,'fontsize',24,'ticklabelinterpreter','tex','FontName','Arial','FontWeight','Bold');
        title(sprintf('\\kappa_{mem}^{bs} = %.2f \\times 10^{%d} mol/m/s/MPa',kappa_bs_mem_a,kappa_bs_mem_exp),'FontSize',20,'Interpreter','tex');
        legend1=legend(sprintf('\\psi_{c-c}'),sprintf('\\psi_{apo}'),'interpreter','tex');
        set(legend1,'Position',[0.157149113714695 0.748901596398512 0.159461545944213 0.147857143765404],'Interpreter','tex','FontSize',20,'FontName','Arial');
        box on;ax=gca;ax.LineWidth = 2;
        drawnow();
        thisFrame=getframe(gcf);
        size(thisFrame.cdata)
        % Write this frame out to a new video file.
        myMovie1(i) = thisFrame;

    end
    
    %% Calculate goxz
    % Note: here, e_xyl, e_i - all "e"s are vapor pressures
    
    e_xyl=p_sat.*exp(v_bar*psi_xylem./(R*leaf_temp)); % equivalent vapor pressure at xylem water potential 
    e_air=p_sat-vpd(end)*1e3; % vapor pressure in the air outside the leaf
    e_i_high_vpd=p_sat.*exp(v_bar*psi_mes_apo_ad_N(:,end)*1e6./(R*leaf_temp)); % equivalent vapor pressure at AquaDust water potential in the ssc
    rh=e_i_high_vpd./p_sat; % relative humidity in the ssc
    e=gs*(e_i_high_vpd-e_air)/p_atm; % calculate transpiration rate e
    goxz_all=e./(e_xyl-e_i_high_vpd)*p_atm; % calculate g_oxz
    goxz_all=goxz_all.*ones(length(psi_xylem),length(psi_xylem));
    goxz=goxz_all(:,1); % extract g_oxz (here, there is no variation of g_oxz with VPD so any column will do)
    
    %% Calculate model prediction mean and max/min to generate uncertainty intervals
    
    % binned rh predictions
    rh_bin1=rh(psi_xylem>-0.4e6);
    rh_bin2=rh(psi_xylem<-0.4e6 & psi_xylem>-0.8e6);
    rh_bin3=rh(psi_xylem<-0.8e6 & psi_xylem>-1.2e6);
    rh_bin4=rh(psi_xylem<-1.2e6 & psi_xylem>=-1.6e6);
    
    % binned goxz predictions
    goxz_bin1=goxz(psi_xylem>-0.4e6);
    goxz_bin2=goxz(psi_xylem<-0.4e6 & psi_xylem>-0.8e6);
    goxz_bin3=goxz(psi_xylem<-0.8e6 & psi_xylem>-1.2e6);
    goxz_bin4=goxz(psi_xylem<-1.2e6 & psi_xylem>=-1.6e6);

    if ii==1 % upper bound on uncertainty interval

        rh_bin1_max=mean(rh_bin1);
        rh_bin2_max=mean(rh_bin2);
        rh_bin3_max=mean(rh_bin3);
        rh_bin4_max=mean(rh_bin4);

        goxz_bin1_max=mean(goxz_bin1);
        goxz_bin2_max=mean(goxz_bin2);
        goxz_bin3_max=mean(goxz_bin3);
        goxz_bin4_max=mean(goxz_bin4);
    
    elseif ii==2 % mean prediction

        rh_bin1_mean=mean(rh_bin1);rh_bin1_pos_err=rh_bin1_max-rh_bin1_mean; % compute difference between upper bound on uncertainty interval and mean
        rh_bin2_mean=mean(rh_bin2);rh_bin2_pos_err=rh_bin2_max-rh_bin2_mean;
        rh_bin3_mean=mean(rh_bin3);rh_bin3_pos_err=rh_bin3_max-rh_bin3_mean;
        rh_bin4_mean=mean(rh_bin4);rh_bin4_pos_err=rh_bin4_max-rh_bin4_mean;

        goxz_bin1_mean=mean(goxz_bin1);goxz_bin1_pos_err=goxz_bin1_max-goxz_bin1_mean; % compute difference between upper bound on uncertainty interval and mean
        goxz_bin2_mean=mean(goxz_bin2);goxz_bin2_pos_err=goxz_bin2_max-goxz_bin2_mean;
        goxz_bin3_mean=mean(goxz_bin3);goxz_bin3_pos_err=goxz_bin3_max-goxz_bin3_mean;
        goxz_bin4_mean=mean(goxz_bin4);goxz_bin4_pos_err=goxz_bin4_max-goxz_bin4_mean;
    
    else % lower bound on uncertainty interval

        rh_bin1_min=mean(rh_bin1);rh_bin1_neg_err=rh_bin1_mean-rh_bin1_min; % compute difference between mean and lower bound on uncertainty interval
        rh_bin2_min=mean(rh_bin2);rh_bin2_neg_err=rh_bin2_mean-rh_bin2_min;
        rh_bin3_min=mean(rh_bin3);rh_bin3_neg_err=rh_bin3_mean-rh_bin3_min;
        rh_bin4_min=mean(rh_bin4);rh_bin4_neg_err=rh_bin4_mean-rh_bin4_min;

        goxz_bin1_min=mean(goxz_bin1);goxz_bin1_neg_err=goxz_bin1_mean-goxz_bin1_min; % compute difference between mean and lower bound on uncertainty interval
        goxz_bin2_min=mean(goxz_bin2);goxz_bin2_neg_err=goxz_bin2_mean-goxz_bin2_min;
        goxz_bin3_min=mean(goxz_bin3);goxz_bin3_neg_err=goxz_bin3_mean-goxz_bin3_min;
        goxz_bin4_min=mean(goxz_bin4);goxz_bin4_neg_err=goxz_bin4_mean-goxz_bin4_min;
    end

end

%% Plotting prediction mean and uncertainty interval

rh_neg_err=[rh_bin1_neg_err;rh_bin2_neg_err;rh_bin3_neg_err;rh_bin4_neg_err];
rh_pos_err=[rh_bin1_pos_err;rh_bin2_pos_err;rh_bin3_pos_err;rh_bin4_pos_err];
rh_means_model=[rh_bin1_mean;rh_bin2_mean;rh_bin3_mean;rh_bin4_mean];

figure(2);
bar(xylem_means,rh_means_model,0.5,'linewidth',2);hold on;
pbaspect([1 1 1]);
er=errorbar(xylem_means,rh_means_model,rh_neg_err,rh_pos_err,'linewidth',2);
er.Color=[0 0 0];                            
er.LineStyle='none';  
ylabel('\it{h}_{ssc}^{apo}');ylim([0.88 1.01]);
xlabel('\psi_{xyl} (MPa)');xlim([-1.7 0]);xticks(flipud(xylem_means));
set(gca,'fontsize',24,'ticklabelinterpreter','tex','FontName','Arial','fontweight','bold');
box on;ax=gca;ax.LineWidth = 2;

goxz_neg_err=[goxz_bin1_neg_err;goxz_bin2_neg_err;goxz_bin3_neg_err;goxz_bin4_neg_err];
goxz_pos_err=[goxz_bin1_pos_err;goxz_bin2_pos_err;goxz_bin3_pos_err;goxz_bin4_pos_err];
goxz_means_model=[goxz_bin1_mean;goxz_bin2_mean;goxz_bin3_mean;goxz_bin4_mean];

figure(3);
bar(xylem_means,goxz_means_model,0.5,'linewidth',2);hold on;
pbaspect([1 1 1]);
er=errorbar(xylem_means,goxz_means_model,goxz_neg_err,goxz_pos_err,'linewidth',2);
er.Color=[0 0 0];                            
er.LineStyle='none';
ylabel('{\itg}_{oxz} (mmol/m^2/s)');ylim([4e2 6e4]);
xlabel('\psi_{xyl} (MPa)');xlim([-1.7 0]);xticks(flipud(xylem_means));
set(gca,'fontsize',24,'ticklabelinterpreter','tex','FontName','Arial','fontweight','bold','yscale','log');
box on;ax=gca;ax.LineWidth = 2;

%% plot h_ssc^apo expt and model together

% rh expt data from Fig. 2c
rh_means_expt=[0.992847001815885;0.946963460686416;0.929521991233521;0.919537568810321];
rh_se_expt=[0.00180002577421479;0.0201784472288282;0.0131077713211281;0.0122904973618615];

figure(4);
b1=bar(xylem_means,[rh_means_expt zeros(length(rh_means_expt),1)],'grouped','linewidth',2);hold on;
x=NaN(2,4);
for i = 1:2
    x(i,:) = b1(i).XEndPoints;
end
errorbar(x',[rh_means_expt zeros(length(rh_means_expt),1)],[rh_se_expt zeros(length(rh_se_expt),1)],[rh_se_expt zeros(length(rh_se_expt),1)],'k','linestyle','none','linewidth',2);

b2=bar(xylem_means,[zeros(length(rh_means_model),1) rh_means_model],'grouped','linewidth',2);hold on;
for i = 1:2
    x(i,:) = b2(i).XEndPoints;
end
errorbar(x',[zeros(length(rh_means_model),1) rh_means_model],[zeros(length(rh_neg_err),1) rh_neg_err],[zeros(length(rh_pos_err),1) rh_pos_err],'k','linestyle','none','linewidth',2);

xlabel('\psi_{xyl} (MPa)');xlim([-1.7 0]);xticks(flipud(xylem_means));
ylabel('\it{h}_{ssc}^{apo}');ylim([0.88 1.01]);
set(gca,'fontsize',24,'ticklabelinterpreter','tex','FontName','Arial','fontweight','bold');
pbaspect([1 1 1]);
legend([b1(1) b2(2)],'Experiment','Model','Position',[0.265625 0.757142857142857 0.283928571428571 0.129761904761905]);
b2(2).FaceColor(1,:)=[0.8500 0.3250 0.0980];
box on;ax=gca;ax.LineWidth = 2;

%% plot goxz expt and model together

% goxz expt data from Fig. 2e
goxz_means_expt=[38556.2258047020;5667.33817184295;2422.37978608609;1541.23262663729];
goxz_se_expt=[13600.6249811751;1367.98145675526;303.839465584318;290.498308255051];

figure(5);
b1=bar(xylem_means,[goxz_means_expt zeros(length(goxz_means_expt),1)],'grouped','linewidth',2);hold on;
x=NaN(2,4);
for i = 1:2
    x(i,:) = b1(i).XEndPoints;
end
errorbar(x',[goxz_means_expt zeros(length(goxz_means_expt),1)],[goxz_se_expt zeros(length(goxz_se_expt),1)],[goxz_se_expt zeros(length(goxz_se_expt),1)],'k','linestyle','none','linewidth',2);

b2=bar(xylem_means,[zeros(length(goxz_means_model),1) goxz_means_model],'grouped','linewidth',2);hold on;
for i = 1:2
    x(i,:) = b2(i).XEndPoints;
end
errorbar(x',[zeros(length(goxz_means_model),1) goxz_means_model],[zeros(length(goxz_neg_err),1) goxz_neg_err],[zeros(length(goxz_pos_err),1) goxz_pos_err],'k','linestyle','none','linewidth',2);

xlabel('\psi_{xyl} (MPa)');xlim([-1.7 0]);xticks(flipud(xylem_means));
ylabel('{\itg}_{oxz} (mmol/m^2/s)');ylim([6e2 8e4]);
set(gca,'fontsize',24,'ticklabelinterpreter','tex','FontName','Arial','fontweight','bold','yscale','log');
pbaspect([1 1 1]);
legend([b1(1) b2(2)],'Experiment','Model','Position',[0.265625 0.757142857142857 0.283928571428571 0.129761904761905]);
b2(2).FaceColor(1,:)=[0.8500 0.3250 0.0980];
box on;ax=gca;ax.LineWidth = 2;

%% save workspace variables

filename='workspace_variables.mat';
save(filename);