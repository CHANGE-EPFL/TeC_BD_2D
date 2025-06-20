%%%%%%%%% this is the file to hotstart the simulation
clear Datam
OutDir=TITLE_SAVE;
Fstep=strcat('Final_step_',TITLE_SAVE);
%%%%%%%%%%%%%%%%%%%%%% to get the correct meteology data
%%%%%%%%%%%% Data loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TR_I=0;
load('Hafren_rain_deposition_Rain_correct_m3_limited.mat')
TR_I(:,:)=0;
%%%%%%% Stations -->  1 Jan 1976 - 31 Dec 2010
load([pwd filesep 'Data_Plynlimon_Tanllwyth_run.mat']);
%%% Pr = 2580 mm/yr
Pre_Rif = Pre'; Ta_Rif = Ta'; Ws_Rif = Ws';
ea_Rif = ea'; Tdew_Rif=Tdew'; esat_Rif = esat';
N_Rif = N;  SAD1_Rif = SAD1'; SAD2_Rif = SAD2'; SAB1_Rif = SAB1'; SAB2_Rif = SAB2';
PARD_Rif = PARD'; PARB_Rif = PARB';
Pr_dist = Pr;
qta_Sta=[Zbas]; %%[m a.s.l.]
%%%%
load([pwd filesep 'Data_Plynlimon_Carreg_Wen_run.mat']);
%%% Pr = 2724 mm/yr
Pre_Rif =[Pre_Rif Pre']; Ta_Rif = [Ta_Rif Ta']; Ws_Rif = [Ws_Rif Ws'];
ea_Rif = [ea_Rif ea']; Tdew_Rif=[Tdew_Rif Tdew']; esat_Rif = [esat_Rif esat'];
N_Rif = [N_Rif N];  SAD1_Rif = [SAD1_Rif SAD1']; SAD2_Rif = [SAD2_Rif SAD2'];
SAB1_Rif = [SAB1_Rif SAB1']; SAB2_Rif = [ SAB2_Rif SAB2'];
PARD_Rif = [PARD_Rif PARD'];  PARB_Rif = [PARB_Rif PARB'];
qta_Sta = [qta_Sta Zbas];
Pr_dist=[Pr_dist Pr];
%%%%
clear Pr Pre Ta Ws ea Tdew esat N SAD1 SAD2 SAB1 SAB2 Zbas PARB PARD Rdif Rsw U
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Discharge --  2093 mm/yr ;; 3217 upper
%%%% I would like to set the start and end day

            index1 = find(Date >= dateNum1, 1, 'first');
            index2 = find(Date >= dateNum2, 1, 'first');
            x1=index1; %%%% as part2
            x2=index2;
%%%%
% x1= 87649;% 1;%% 2month warm up the model
% x2= 87649+8760*2;%8760;%175296;% 8760;%306816;%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%
Date=Date(x1:x2);
Pr_dist=Pr_dist(x1:x2,:);
TR_I=TR_I(x1:x2,:);
%%%%%%%%%%%%% strong rain fo test
%Pr_dist(1:100,:)=50;
%%%%%%%%%%%%%
Pre_Rif=Pre_Rif(x1:x2,:);
Ta_Rif=Ta_Rif(x1:x2,:);
Ws_Rif=Ws_Rif(x1:x2,:); ea_Rif=ea_Rif(x1:x2,:);
Tdew_Rif=Tdew_Rif(x1:x2,:); esat_Rif=esat_Rif(x1:x2,:);
%%%%%%%
SAD1_Rif=SAD1_Rif(x1:x2,:); SAB2_Rif=SAB2_Rif(x1:x2,:);
SAD2_Rif=SAD2_Rif(x1:x2,:); SAB1_Rif=SAB1_Rif(x1:x2,:);
N_Rif=N_Rif(x1:x2,:);
PARB_Rif =PARB_Rif(x1:x2,:); PARD_Rif = PARD_Rif(x1:x2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_bef= 0; t_aft= 1;
%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
Ds=esat_Rif-ea_Rif; %% [Pa] Vapor Pressure Deficit
Ds(Ds<0)=0;
%%%%%%%%%%%%%%%%%%%%%%%%% 330-380
load([pwd filesep 'Inputs' filesep 'Ca_Data.mat']);
d1 = find(abs(Date_CO2-Date(1))<1/36);d2 = find(abs(Date_CO2-Date(end))<1/36);
Ca=Ca(d1:d2);
clear d1 d2 Date_CO2


%%%%
N_time_step=x2-x1;%8760;%87648;% 8760;
Nd_time_step = ceil(N_time_step/24)+1;

%%%% general
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE
%%%%%%%%
L_day=zeros(length(Datam),1);
for j=2:24:length(Datam)
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')