%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Datam
OutDir=TITLE_SAVE;
Fstep=strcat('Final_step_',TITLE_SAVE);
%%%%%%%%%%%% Data loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% 
load([pwd filesep 'Data_Erlenbach_run.mat']);
%%% put all variables in one colume

qta_Sta=[Zbas]; %%[m a.s.l.]
%%%%

%%%%
%clear Pr Pre Ta Ws ea Tdew esat N SAD1 SAD2 SAB1 SAB2 Zbas PARB PARD Rdif Rsw U
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
% Pr_dist=Pr_dist(x1:x2,:);

%%%%%%%%%%%%% strong rain fo test
Pr=Pr(x1:x2);
%Pr(1:20,:)=10;
%%%%%%%%%%%%%
Pre=Pre(x1:x2);
Ta=Ta(x1:x2);
Ws=Ws(x1:x2); ea=ea(x1:x2);
Tdew=Tdew(x1:x2); esat=esat(x1:x2);
%%%%%%%
SAD1=SAD1(x1:x2); SAB2=SAB2(x1:x2);
SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2);
N=N(x1:x2);
PARB =PARB(x1:x2); PARD = PARD(x1:x2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_bef= 0.25; t_aft= 0.75;
%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
Ds=esat-ea; %% [Pa] Vapor Pressure Deficit
Ds(Ds<0)=0;
%%%%%%%%%%%%%%%%%%%%%%%%% 
load([pwd filesep 'Inputs' filesep 'Ca_Data.mat']);
d1 = find(abs(Date_CO2-Date(1))<1/36);d2 = find(abs(Date_CO2-Date(end))<1/36);
Ca=Ca(d1:d2);
clear d1 d2 Date_CO2
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ws(Ws<=0)=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_time_step=x2-x1;%8760;%87648;% 8760;
Nd_time_step = ceil(N_time_step/24)+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% GENERAL PARAMETER %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%