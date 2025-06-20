%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% MAIN_FRAME SPATIAL TETHYS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Zbio,rsd,aR,Zs,...
    EvL_Zs,Inf_Zs,Bio_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Ks_Zs,Dz,...
    ms,Kbot,Krock,zatm,...
    Ccrown,Cbare,Crock,Curb,Cwat,...
    Color_Class,OM_H,OM_L,PFT_opt_H,PFT_opt_L,d_leaf_H,d_leaf_L,...
    SPAR,Phy,Soil_Param,Interc_Param,SnowIce_Param,VegH_Param,VegL_Param,fpr,...
    VegH_Param_Dyn,VegL_Param_Dyn,...
    Stoich_H,aSE_H,Stoich_L,aSE_L,fab_H,fbe_H,fab_L,fbe_L,...
    ZR95_H,ZR95_L,In_max_urb,In_max_rock,K_usle,...
    Deb_Par,Zs_deb,... 
    Sllit,Kct,ExEM,ParEx_H,Mpar_H,ParEx_L,Mpar_L]=PARAMETERS_ALL_Plynlimon(ANSWER,Psan,Pcla,Porg,md_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cur_dir=cd;
% cd(code_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PARAMETERS SOILS AND VEGETATION  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ANSWER; 
%%%%
%%% Sitka Spruce  1
%%% Grassland   2
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
ms=10; %%% Soil Layer
zatm =[30.0 3.0 ] ; %% Reference Height
zatm = zatm(ANSWER);
fpr = 1;
%%%%
aR =50;
%Kh=Ks*aR;
Kbot = [0.01 0.01 ]; %% [mm/h] Conductivity at the bedrock layer
Krock = [NaN NaN ]; %% [mm/h] Conductivity of Fractured Rock
Kbot= Kbot(ANSWER);
Krock=Krock(ANSWER);
%%%%%%%%%% SOIL INPUT
Color_Class = 0;
%%%%
switch ANSWER
     %case {1,2,3,8,9}
    case {1}
        %%%% LAND COVER PARTITION  100  Sitka Spruce   
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [ 1 0 ]>0;
     %case 6
    case {2}
        %%%% LAND COVER PARTITION 100 GRassland 
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; 
        % Ccrown = [0 1];
        Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [0 1]>0;
    otherwise
        disp('INDEX FOR SOIL VEGETATION PARAMETER INCONSISTENT')
        return
end
%%%%%%%%%%%%%%%%%%%
SPAR=2; %%% SOIL PARAMETER TYPE
%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s =lan_s*ones(1,ms);
cv_s = cv_s*ones(1,ms);
%%%
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%%%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
%[s_SVG,bVG]=Soil_parameters_VG(Phy,Osat,Ohy,nVG,alpVG,0);
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
Zs = [ 0 10 20 50 100 150 200 300 400 500 600 ]; %%% [ms+1]
if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end
Zdes = 10; %%% Evaporation depth
Zinf = 10; %%% Infiltration depth
Zbio = 250; 
[EvL_Zs]=Evaporation_layers(Zs,Zdes); %%% Evaporation Layer fraction 
[Inf_Zs]=Evaporation_layers(Zs,Zinf); %%% Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio); 
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
Dz=zeros(1,ms); 
for ii = 1:ms
    if ii>1
        Dz(ii)= (dz(ii)+ dz(ii-1))/2; %%% Delta Depth Between Middle Layer  [mm]
    else
        Dz(ii)=dz(1)/2; %%% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end
%%%%%%%%%%%%%%%%% OTHER PARAMETER
In_max_urb=5;
In_max_rock=2; %% [mm]
%%%%%%%%%%%%% SNOW PARAMETER
TminS=-0.8;%% Threshold temperature snow
TmaxS= 2.8;%% Threshold temperature snow
ros_max1=580; %600; %%% [kg/m^3]
ros_max2=300; %450; %%% [kg/m^3]
Th_Pr_sno = 8.0; %%% [mm/day] Threshold Intensity of snow to consider a New SnowFall
%%%% ICE Parameter
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice 
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice 
Aice = 0.35; %% [-] Ice albedo 
WatFreez_Th = -8; %% [°C] Threshold for freezing lake water
dz_ice = 0.54; %% [mm / h] Water Freezing Layer progression without snow-layer  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExEM = [0.0 0.0 ];
ExEM = ExEM(II); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS VEGETATION
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile 
%%%
ZR95_H = [500 0 ]; %% [mm]
ZR95_L = [0 250 ]; %% [mm]
ZR50_H = [NaN NaN];
ZR50_L = [NaN NaN ];
ZRmax_H = [NaN NaN ];
ZRmax_L = [NaN NaN ];
ZR95_H =ZR95_H(II); ZR50_H =ZR50_H(II); ZRmax_H =ZRmax_H(II);
ZR95_L =ZR95_L(II); ZR50_L =ZR50_L(II); ZRmax_L =ZRmax_L(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
%5 Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_H_In= [0.1 0.1 ]; %%[mm/LAI]
Sp_LAI_L_In= [0.2 0.2 ]; %%[mm/LAI]
Sp_LAI_H_In =Sp_LAI_H_In(II);
Sp_LAI_L_In =Sp_LAI_L_In(II);
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [0.25 0.25 ]; %%[cm]
d_leaf_L= [0.8 0.8 ];  %% [cm]
d_leaf_H =d_leaf_H(II);
d_leaf_L =d_leaf_L(II);
%%%%%%%% Biochemical parameter
KnitH=[0.35 0.35]; %%% Canopy Nitrogen Decay
KnitL=[0.15 0.15 ];
mSl_H = [0 0 ];%% [m2 PFT /gC]  Linear increase in Sla with LAI 
mSl_L = [0.0 0.0];  % 
KnitH =KnitH(II);  mSl_H =mSl_H(II); 
KnitL =KnitL(II);  mSl_L =mSl_L(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Photosynthesis Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FI_H=[0.081 0.081  ];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[800 800 ] ; %%[Pa]
a1_H=[7 7 ];
go_H=[0.01 0.01 ];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3 3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649 0.649  ];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[72 72 ]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf Inf Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[2.3 2.3 ]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%
FI_H=FI_H(II); Do_H=Do_H(II); a1_H=a1_H(II); go_H=go_H(II);
CT_H=CT_H(II); DSE_H=DSE_H(II); Ha_H=Ha_H(II); gmes_H=gmes_H(II);
rjv_H=rjv_H(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_L=[0.081 0.081 ];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[1000 1000 ]; %%[Pa]
a1_L=[6 6 ];
go_L=[0.01 0.01  ];% % [mol / s m^2] minimum Stomatal Conductance
CT_L=[3 3 3];  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[0.649 0.649 ];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[55 55]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[Inf Inf Inf ]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L=[2.2 2.2 ]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%
FI_L=FI_L(II); Do_L=Do_L(II); a1_L=a1_L(II); go_L=go_L(II);
CT_L=CT_L(II); DSE_L=DSE_L(II); Ha_L=Ha_L(II); gmes_L=gmes_L(II);
rjv_L=rjv_L(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi_sto_50_H = [-2.5 -2.5] ;%% [MPa]  Water Potential at 50% loss conductivity
Psi_sto_00_H = [-0.5 -0.5 ]; %% [MPa]  Water Potential at 2% loss conductivity
%%% Leaf
PsiL00_H = [ -1.0 -1.0 ] ;%%[MPa]  Water Potential at 50% loss conductivity
PsiL50_H =  [-3.2 -3.2 ]; %% [MPa]  Water Potential at 2% loss conductivity
Kleaf_max_H = [10 10 ] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200 1200 ];  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = [ 15 15 ] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000 80000 ];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-5.0 -5.0 ]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [150 150]; %%%% [kg / m^3 sapwood MPa]
%%------------------------
%
Psi_sto_00_L = [-0.5 -0.5];%  %% [MPa]  Water Potential at PLCs loss conductivity
Psi_sto_50_L = [ -3.0 -3.0];%  %% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L =  [-0.9 -0.9]; %% [MPa]  Water Potential at PLCs% loss conductivity
PsiL50_L = [-4.0 -4.0] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [5 5] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [1200 1200];  %%%  [500 - 3000]%  [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L = [0.0 0.0] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [ 80000 80000];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [ -4.5 -4.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [150 150]; %%% [kg / m^3 sapwood MPa]
%%%
Psi_sto_50_H =Psi_sto_50_H(II);  Psi_sto_00_H =Psi_sto_00_H(II); 
PsiL00_H = PsiL00_H(II); PsiL50_H=PsiL50_H(II);  Kleaf_max_H=Kleaf_max_H(II); 
Cl_H=Cl_H(II); Axyl_H=Axyl_H(II); Kx_max_H=Kx_max_H(II); PsiX50_H=PsiX50_H(II); Cx_H=Cx_H(II); 
Psi_sto_50_L =Psi_sto_50_L(II);  Psi_sto_00_L =Psi_sto_00_L(II); 
PsiL00_L = PsiL00_L(II); PsiL50_L=PsiL50_L(II);  Kleaf_max_L=Kleaf_max_L(II); 
Cl_L=Cl_L(II); Axyl_L=Axyl_L(II); Kx_max_L=Kx_max_L(II); PsiX50_L=PsiX50_L(II); Cx_L=Cx_L(II); 

%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L); 

%%%% Growth Parameters 
PsiG50_H= [-0.5 -0.5 ];  %%[MPa]
PsiG99_H= [-2.5 -2.5 ];  %%[MPa]
gcoef_H = [3.5 3.5 ]; % [gC/m2 day]
%%------  
PsiG50_L= [-3.0 -3.0]; 
PsiG99_L= [-4.0 -4.0]; 
gcoef_L =[3.5 3.5]; % [gC/m2 day]
%%%%%
PsiG50_H=PsiG50_H(II); PsiG99_H=PsiG99_H(II); gcoef_H=gcoef_H(II); 
PsiG50_L=PsiG50_L(II); PsiG99_L=PsiG99_L(II); gcoef_L=gcoef_L(II); 

OPT_PROP_H =[2 0 ];
OPT_PROP_L =[0 13 ];
OPT_PROP_H=OPT_PROP_H(II);
OPT_PROP_L=OPT_PROP_L(II);
for i=1:cc
    %%%%%%%% Vegetation Optical Parameter
    [PFT_opt_H(i)]=Veg_Optical_Parameter(OPT_PROP_H(i));
    [PFT_opt_L(i)]=Veg_Optical_Parameter(OPT_PROP_L(i));
end

OM_H=1;
OM_L=1;
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.011 0.011 ]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [42 42 ]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
r_H = [0.058 0.058 ];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_H= [0.25 0.25  ]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
aSE_H= [0 0 ]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_H= [1/250 1/250]; %%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  [78/365 78/365  ]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H = [-20 -20 ]; %% [°C] Cold Leaf Shed
drn_H=  [1/700 1/700 ]; %% turnover root  [1/d]
dsn_H= [1/900 1/900]; % normal transfer rate sapwood [1/d]
age_cr_H= [950 950]; %% [day] Critical Leaf Age
Bfac_lo_H= [0.99 0.99 ]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN NaN  ]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [5.5 5.5]; %% Mean Temperature for Leaf onset
Tls_H = [NaN NaN  ]; %% Mean Temperature for Leaf Shed
PAR_th_H= [NaN NaN ];
dmg_H= [30 30  ]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.001 0.001];
Trr_H = [0.25 0.25]; %% Translocation rate [gC /m^2 d]
mjDay_H = [180 180 ]; %% Maximum Julian day for leaf onset
LDay_min_H =[12.8 12.8 ]; %% Minimum Day duration for leaf onset
LtR_H = [0.7 0.7]; %%% Leaf to Root ratio maximum
Mf_H= [1/80 1/80  ]; %% fruit maturation turnover [1/d]
Wm_H= [1/(40.1*365) 1/(40.1*365)  ] ; % wood turnover coefficient [1/d]. % result from BG has to be rewtitten here!!
eps_ac_H = [0.1 0.1]; %% Allocation to reserve parameter [0-1]
LDay_cr_H = [9.45 9.45 ]; %%%  Threshold for senescence day light [h]
Klf_H =[1/40 1/40 ]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74 0.74  ]; %% fraction above-ground sapwood and reserve
fbe_H = [0.26 0.26 ]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1 0.1 ]; %% 
%%%%

Sl_H =Sl_H(II); Nl_H=Nl_H(II); 
r_H=r_H(II); gR_H=gR_H(II); aSE_H=aSE_H(II); dd_max_H=dd_max_H(II);
dc_C_H=dc_C_H(II); Tcold_H=Tcold_H(II); drn_H=drn_H(II);
dsn_H=dsn_H(II);  age_cr_H=age_cr_H(II);
Bfac_lo_H=Bfac_lo_H(II); Bfac_ls_H=Bfac_ls_H(II);
Tlo_H = Tlo_H(II);  Tls_H=Tls_H(II);
dmg_H = dmg_H(II); LAI_min_H=LAI_min_H(II);
Trr_H = Trr_H(II);  mjDay_H=mjDay_H(II);
LDay_min_H= LDay_min_H(II); LtR_H =LtR_H(II);
Mf_H= Mf_H(II);  Wm_H= Wm_H(II);  eps_ac_H = eps_ac_H(II);
LDay_cr_H = LDay_cr_H(II);  Klf_H = Klf_H(II);
fab_H = fab_H(II); fbe_H = fbe_H(II); ff_r_H = ff_r_H(II);

for i=1:cc
    [Stoich_H(i)]=Veg_Stoichiometric_Parameter(Nl_H(i));
    [ParEx_H(i)]=Exudation_Parameter(0);
    [Mpar_H(i)]=Vegetation_Management_Parameter;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L = [ 0.022 0.022 ]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L= [23 23  ]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
r_L = [ 0.062 0.062 ];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_L= [0.25 0.25   ]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
aSE_L= [2 2  ]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_L= [ 1/45 1/45 ];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L =  [7/365 7/365  ]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_L = [-2.0 -2.0 ]; %% [°C] Cold Leaf SLed
drn_L=  [1/450 1/450  ]; %% turnover root  [1/d]
dsn_L= [ 1/365 1/365 ]; % normal transfer rate sapwood [1/d]
age_cr_L= [180 180  ]; %% [day] Critical Leaf Age
Bfac_lo_L= [0.99 0.99  ]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN NaN  ]; %% Leaf Shed Water Stress [0-1]
Tlo_L = [3.5 3.5 ]; %% Mean Temperature for Leaf onset
Tls_L = [NaN NaN  ]; %% Mean Temperature for Leaf Shed
PAR_th_L= [NaN NaN ];
dmg_L= [20 20  ]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_L = [0.1 0.1  ];
Trr_L = [3.0 3.0 ]; %% Translocation rate [gC /m^2 d]
mjDay_L = [250 250 ]; %% Maximum Julian day for leaf onset
LDay_min_L =[10.0 10.0 ]; %% Minimum Day duration for leaf onset
LtR_L = [0.4 0.4 ]; %%% Leaf to Root ratio maximum
Mf_L= [1/50 1/50 ]; %% fruit maturation turnover [1/d]
Wm_L= [0 0 ] ; % wood turnover coefficient [1/d]
eps_ac_L = [0.2 0.2  ]; %% Allocation to reserve parameter [0-1]
LDay_cr_L = [11.4 11.4 ]; %%%  Threshold for senescence day light [h]
Klf_L =[1/50 1/50  ]; %% Dead Leaves fall turnover [1/d]
fab_L = [ 0 0  ]; %% fraction above-ground sapwood and reserve
fbe_L = [ 1 1  ]; %% fraction below-ground sapwood and reserve
ff_r_L= [0.1 0.1 ]; %% 
%%%
Sl_L =Sl_L(II); Nl_L=Nl_L(II);
r_L=r_L(II); gR_L=gR_L(II); aSE_L=aSE_L(II); dd_max_L=dd_max_L(II);
dc_C_L=dc_C_L(II); Tcold_L=Tcold_L(II); drn_L=drn_L(II);
dsn_L=dsn_L(II);  age_cr_L=age_cr_L(II);
Bfac_lo_L=Bfac_lo_L(II); Bfac_ls_L=Bfac_ls_L(II);
Tlo_L = Tlo_L(II);  Tls_L=Tls_L(II);
dmg_L = dmg_L(II); LAI_min_L=LAI_min_L(II);
Trr_L = Trr_L(II);  mjDay_L=mjDay_L(II);
LDay_min_L= LDay_min_L(II); LtR_L =LtR_L(II);
Mf_L= Mf_L(II);  Wm_L= Wm_L(II);  eps_ac_L = eps_ac_L(II);
LDay_cr_L = LDay_cr_L(II);  Klf_L = Klf_L(II);
fab_L = fab_L(II); fbe_L = fbe_L(II); ff_r_L = ff_r_L(II);

for i=1:cc
    [Stoich_L(i)]=Veg_Stoichiometric_Parameter(Nl_L(i));
    [ParEx_L(i)]=Exudation_Parameter(0);
    [Mpar_L(i)]=Vegetation_Management_Parameter;
end

%%%%
Vmax_H = [62 62 ]; %
Vmax_L = [58 58 ]; %
Vmax_H =Vmax_H(II); Vmax_L =Vmax_L(II);
%%%%%%%%%%%%%%%%%
Restating_parameters;
if dbThick>5
    nst=md_max;
    k=(dbThick/5)^(1/(nst-1));
    Zs_deb = [0 5*k.^(1:nst-1)]; %% [mm]
    clear nst k
else
end
%cd(cur_dir)
end