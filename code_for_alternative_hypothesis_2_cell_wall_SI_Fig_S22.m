% code_for_alternative_architecture_cell_wall_SI_Fig_S21.m
% Description: This code executes the mathematical model developed in SI
% Section S10 in the paper titled, "Loss of plasma membrane conductance of 
% outside-xylem zone explains non-stomatal control of transpiration" for 
% alternative architecture 2: cell wall is the dominant hydraulic 
% conductor and undergoes loss of hydraulic conductance under drought stress. 
% We generate 3 predictions (mean and bounds of uncertainty 
% interval) with npoints (see below) being the number of points in the 
% discretized domain for each prediction. npoints are distributed across 
% four ranges in psi_xyl: -0.4 to 0; -0.8 to -0.4; -1.2 to -0.8 and -1.6 to
% -1.2 (MPa). The prediction mean and bounds of uncertainty interval is 
% then compared to the data mean +/- S.E. 
% This code uses the following MATLAB toolboxes:
% 1. Optimization Toolbox
% Code developed by: Sabyasachi Sen
% Email: ss3945@cornell.edu
% Affiliation: Cornell University
% Created: February 8, 2025
% Last Modified: February 10, 2025
% Copyright (c) 2025 Sabyasachi Sen
% If you use this code in your research, please cite as:
% "Sabyasachi Sen (2025). 
% code_for_alternative_architecture_cell_wall_SI_Fig_S21.m
% Retrieved from manuscript titled: Loss of plasma membrane conductance of 
% outside-xylem zone explains non-stomatal control of transpiration"
% This code is provided under the MIT License.
% For full license terms, see: https://opensource.org/licenses/MIT

clear;
close all;
clc;

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
tortuosity_arr=[1, 2, 3]; % uniform tortuosity (dimensionless)
k_cw_max_arr=[2.9e-6, 2.7e-6, 2.5e-6]; % maximum cell wall conductivity in mol/m/s/MPa
k_cw_trigger_arr=[-0.27e6, -0.25e6, -0.23e6]; % xylem stress at which 'characteristic' loss of cell wall conductivity occurs (in Pa)
k_cw_decline_sensitivity_arr=[0.8e-5, 1e-5, 1.2e-5]; % determines sensitivity of cell wall conductivity to xylem stress in MPa^(-1)
% see vapor conductivity for minimum cell wall conductivity

%% Constant leaf anatomical parameters (see Table S4 in SI Section S10)

L_total=2e-4; % total thickness of leaf in m
L_ad=L_total/2; % Assume half-thickness of leaf 
L_ab=L_total/2; % Assume half-thickness of leaf 
A_cw=0.08e-9; % x-sectional area for flow in cell wall (in m2)

%% Constant transport parameters (see Table S5 in SI Section S10)

kappa_mem=1e-5; % hydraulic conductivity of membranes at cell-air space interface in mol/m/s/MPa: Rockwell review undersaturation - use higher end of L_p 
k_cc=5e-8; % hydraulic conductivity of cell-to-cell path in mol/m/s/MPa: Rockwell review undersaturation - use higher end of L_p 
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
coordinates_ad=coords_ad*1e6; % coordinates for adaxial section: microns for plotting
coordinates_ab=coords_ab*1e6; % coordinates for abaxial section: microns for plotting

%% Loop through 3 sets of parameters

for ii=1:3

    % Variable leaf anatomical parameters and transport parameters (see Tables S4 and S5 in SI Section S10)

    tortuosity=tortuosity_arr(ii); % tortuosity
    k_vap=c*D*v_bar*kai_0/(R*leaf_temp)/tortuosity/1e-6; % hydraulic conductivity of vapor in mol/m/s/MPa; this assumes driving force for transport is water potential. Ref: Rockwell 2014
    
    A_cc=A_cc_arr(ii); % see SI Section S10 for choice of A_cc
    A_vap=A_vap_arr(ii); % see SI Section S10 for choice of A_vap
                
    %% Tissue hydraulic conductivity parameters 
    
    % Cell-to-cell hydraulic parameters for current iteration
    k_cw_max=k_cw_max_arr(ii); 
    k_cw_decline_sensitivity=k_cw_decline_sensitivity_arr(ii); 
    k_cw_trigger=k_cw_trigger_arr(ii); 
    k_cw_min=k_vap; 
    % sigmoid function
    k_cw_arr=k_cw_max*1./(1+exp(-k_cw_decline_sensitivity*(psi_xylem-k_cw_trigger)))+k_cw_min; 
    
    % plot kappa_mem vs psi_xyl
    figure(1);hold on;
    plot(psi_xylem/1e6,k_cw_arr,'k-','linewidth',2);hold on;
    pbaspect([1 1 1]);
    ylabel('k_{c-c} (mol/m/s/MPa)');ylim([8e-8 5e-6]);
    xlabel('\psi_{xyl} (MPa)');xlim([-1.7 0]);xticks(flipud(xylem_means));
    set(gca,'fontsize',24,'ticklabelinterpreter','tex','fontname','Arial','fontweight','bold','yscale','log');
    box on;ax=gca;ax.LineWidth = 2;
    
    %% Solution
    
    % Initialize arrays that will contain the solutions: psi_cc and psi_apo (in Pa)

    psi_cc_ad_N=zeros(length(psi_xylem),length(rel_hum)); % psi_cc in adaxial domain
    psi_cc_ab_N=zeros(length(psi_xylem),length(rel_hum)); % psi_cc in abaxial domain
    psi_apo_ad_N=zeros(length(psi_xylem),length(rel_hum)); % psi_apo in adaxial domain
    psi_apo_ab_N=zeros(length(psi_xylem),length(rel_hum)); % psi_apo in abaxial domain
    
    for i=1:length(psi_xylem)
    
        k_cw=k_cw_arr(i); % cell wall conductivity 
        
        for j=1:length(rel_hum)
    
            % see SI Section S10 for details of these constants
            % a1, a2, b1, and b2 multiply the derivates in SI equations S30 and S31

            % Composite parameters for ease of computation
            c_cc=kappa_mem/(k_cc*A_cc); 
            c_apo=kappa_mem/(k_cw*A_cw + k_vap*A_vap);
            
            a1=c_cc;
            a2=c_apo;
            r3=+sqrt(a1+a2);
            
            % solve for the unknown x's ('alpha's in SI equations S44 and S45)            
            init=[-1e6,-1e8,-10,-1e6,-2e6,-2e8,-20,-2e6];
            x=fsolve(@(x)solver_alternative_hypothesis_2(x,psi_xylem(i),rel_hum(j),a1,a2,r3,L_ad,L_ab,gs_ad,gs_ab,p_sat,v_bar,R,leaf_temp,k_vap,k_cw,p_atm,A_cc,A_vap,A_cw),init);

            psi_cc_ad_arr=x(1) + x(2)*coords_ad + x(3)*exp(r3*coords_ad) + x(4)*exp(-r3*coords_ad);
            psi_apo_ad_arr=x(1) + x(2)*coords_ad - a2*x(3)*exp(r3*coords_ad)/a1 - a2*x(4)*exp(-r3*coords_ad)/a1;
            psi_cc_ab_arr=x(5) + x(6)*coords_ab + x(7)*exp(r3*coords_ab) + x(8)*exp(-r3*coords_ab);
            psi_apo_ab_arr=x(5) + x(6)*coords_ab - a2*x(7)*exp(r3*coords_ab)/a1 - a2*x(8)*exp(-r3*coords_ab)/a1;        
            
            % Terminal apoplasm and symplasm potentials on adaxial side
            psi_cc_ad_N(i,j)=psi_cc_ad_arr(end)/1e6; % divide by 1e6 to convert from Pa to MPa
            psi_apo_ad_N(i,j)=psi_apo_ad_arr(end)/1e6; % divide by 1e6 to convert from Pa to MPa
            
            % Terminal apoplasm and symplasm potentials on abaxial side 
            psi_cc_ab_N(i,j)=psi_cc_ab_arr(end)/1e6; % divide by 1e6 to convert from Pa to MPa
            psi_apo_ab_N(i,j)=psi_apo_ab_arr(end)/1e6; % divide by 1e6 to convert from Pa to MPa
            
        end
        
        k_cw_exp=floor(log10(k_cw));
        k_cw_a=k_cw/(10^k_cw_exp);

        % Plot psi_apoplast vs z
        figure(2);
        plot([psi_cc_ab_arr; psi_cc_ad_arr]/1e6,[coordinates_ab; coordinates_ad],[psi_apo_ab_arr; psi_apo_ad_arr]/1e6,[coordinates_ab; coordinates_ad],'linewidth',4);
        pbaspect([1 1 1]);
        xlabel('\psi(MPa)','interpreter','tex');xlim([-20 0]);
        ylabel('z (\mu m)','interpreter','tex');ylim([-110 110]);
        set(gca,'fontsize',24,'ticklabelinterpreter','tex','FontName','Arial','FontWeight','Bold');
        title(sprintf('k_{cw} = %.2f \\times 10^{%d} mol/m/s/MPa',k_cw_a,k_cw_exp),'FontSize',20,'Interpreter','tex');
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
    e_i_high_vpd=p_sat.*exp(v_bar*psi_apo_ad_N(:,end)*1e6./(R*leaf_temp)); % equivalent vapor pressure at AquaDust water potential in the ssc
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
ylabel('\it{h}_{ssc}^{apo}');ylim([0.9 1.01]);
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
ylabel('\it{h}_{ssc}^{apo}');ylim([0.9 1.01]);
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