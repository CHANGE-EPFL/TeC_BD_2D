%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  ROUTING_MODULE_BG             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[TR_sur,P,LEAK_nt,N_check,q_runon,q_channel_out,Qi_in,Slo_pot,Q_exit,Qsub_exit,...
    TR_exit,...
    TR_exitS,...
    TRQexit,...
    T_pot,QpointH,QpointC,UpointH,UpointC,Utot_H,Utot_C]= ROUTING_MODULE_BG_rain(N_nt,TR_I,TR_sur,QsurR,Asur,Ared,EG,If,Pr,P,Vtm1,LK,dt,dth,Rd,Rh,Qi_out,q_channel_in,...
    cellsize,Area,DTM,NMAN_H,NMAN_C,MRough,WC,SN,T_flow,T_potI,Slo_top,ms_max,POT,ZWT,OPT_HEAD,Xout,Yout)


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% OUTPUT
%%% TR_sur [g/mm]: nutrient concentration in surface flow
%%% P [g/m2]: 8 dissolved P pool to rout, the last two are all DOC. 
%%% LEAK_nt [g/m2]: nutrient leakage to bedrock
%%% N_check: nutrient mass balance check for each nutrient

%%%%%%%%%%%%
%%% q_runon  [mm]  %% Runon
%%% Qi_in, [mm] %% Subsurface Lateral flow
%%% Slo_pot [fraction] %% Slope of Hydraulic head
%%% Q_exit [mm] %% discharge from domain surface
%%% Qsub_exit [mm] %% discharge from domain subsurface
%%% q_channel_out [mm] %%% Water in channels
%%% T_pot -- new flow direction cell of matrixs for subsurface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% INPUT
%%% N_nt: number of dissolved nutrient
%%% Asur,Ared: related to slope, used for water balance
%%% EG [mm]: soil evaporation
% % % If [mm]: infiltration
% % % Pr [mm]: precipitation
% % % Vtm1 [mm]: soil water volume before the hydrology module
% % % LK [mm]: water leakage
% % % TR_I [g/mm] : chemical concentration in the rain

%%%%%%%%%%%%
%%% dt [s] time step
%%% Rd Dunne Runoff  [mm]
%%% Rh Horton Runoff [mm]
%%% Qi_out  Lateral Subsurface [mm/h]
%%% q_channel_in [mm] Water in channels in
%%%% cellsize [m]
%%% Area [m^2] watershed area
%%% DTM
%%% NMAN_H  [s/(m^1/3)] Manning Coefficient hillslope
%%% NMAN_C  [s/(m^1/3)] Manning Coefficient channels
%%% WC [m] Channel width
%%% SN [-] stream network identifier
%%% T_flow [ Sparse mn x mn]  Flow matrix surface
%%% T_potI  ms cells [ Sparse mn x mn]  Flow matrixs for subsurface
%%% Slo_top [ fraction] Topographic Slope
%%% ms_max Soil layers
%%% POT Head in a cell [mm]
%%% Zwt [mm] water table depth %%%
%%% OPT_HEAD options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m_cell,n_cell]=size(DTM);
Qsur=Rh+Rd; %% [mm] %% Hillslope surface water
Qsur(Qsur<0)=0; %%%% Numerical Instability Issue
Qi_out = Qi_out*dth;  %% Subsurface flow [mm]
opt_cons_CUE=1;
[BiogeoPar]=Biogeochemistry_Parameter(opt_cons_CUE);
%
TR_I=repmat(TR_I,n_cell*m_cell,1);
TR_I=reshape(TR_I,m_cell,n_cell,7);
%%%%%%%%% ice is considerd as pure water
%%%%%V=Vtm1+sum((Oicetm1).*dz - Vice)
P_old=P; %% this is the P before the BG part

% TR_I=0; % deposition
ft_up = 1;% 0.9 ;%coefficient of evapoconcentration. 1=no evapoconcentration (=0.9-1 in Benettin et al., 2015)
Pr_old=Pr; % [mm/h]  Precipitation
Pr_old=reshape(Pr_old,m_cell,n_cell);

Pr= Pr + QsurR; %%% water on the surface
Pr=reshape(Pr,m_cell,n_cell);
QsurR=reshape(QsurR,m_cell,n_cell);

TR_If  = (Pr_old.*TR_I + (QsurR).*TR_sur)./(Pr) ; % nutrient concentration in infiltration [g/mm]
logic_index = repmat(Pr==0, 1, 1, N_nt);
TR_If(logic_index) = 0;


V_total=sum(Vtm1,2)+If; % volume after infiltration

aNH4 =BiogeoPar.aNH4;
aNO3 = BiogeoPar.aNO3;
aP = BiogeoPar.aP;
aK = BiogeoPar.aK;
aDON= BiogeoPar.aDON;
aDOP= BiogeoPar.aDOP;
aDOC = BiogeoPar.aDOC;
solu_coeff=[aNH4;aNO3;aP;aK;aDON;aDOP;aDOC];

%%%%%%%% nutrient concentration in the soil, assume the concentration is
%%%%%%%% same in each soil layer, before infiltration
TR_avg=zeros(m_cell*n_cell,N_nt);
for i=1:min(N_nt,6)
    TR_avg(:,i)=max(solu_coeff(i)*P(:,i)*cellsize*cellsize,0)./sum(Vtm1,2); %[g/mm] average for each layer,before infiltration
end
if N_nt==7
    TR_avg(:,7)=max(solu_coeff(7)*(P(:,7)+P(:,8))*cellsize*cellsize,0)./sum(Vtm1,2);
end


% average nutrient concentration with infiltration
TR_avg = (sum(Vtm1,2).*TR_avg + If.*reshape(TR_If,m_cell*n_cell,N_nt))./V_total; %% [g/mm] , 
TR_avg(V_total==0,:)=0;



%%%%% mass in from infiltration, update P
DOC_fr=P_old(:,7)./(P_old(:,7)+P_old(:,8));
DOC_fr(isnan(DOC_fr))=0;

for i =1:min(N_nt,6)
    P(:,i)=P(:,i)+If.*reshape(TR_If(:,:,i)/cellsize^2,m_cell*n_cell,1);% g to g/m2
end
if N_nt==7
    P(:,7)=P(:,7)+If.*reshape(TR_If(:,:,7)/cellsize^2,m_cell*n_cell,1).*DOC_fr;
    P(:,8)=P(:,8)+If.*reshape(TR_If(:,:,7)/cellsize^2,m_cell*n_cell,1).*(1-DOC_fr);
end

%%%%%%%%%%
%%% Volume second update, consider E [mm]
V_total=V_total-EG*dth./(Asur.*Ared); %Asur Ared here to keep water balance (same in hydrology module)
V_total(V_total<=0) = 0;
V_total(isnan(V_total))=0;


TR_up = EG.*dth./(Asur.*Ared)*ft_up.*TR_avg;           % [C] Tracer uptaken by ET
TR_up(V_total==0,:) = 0;              %
TR1_ETup = sum(TR_up); % [C] total tracer out for ET, but I did not set it as output

%%%%%%%%%%% mass out for ET
for i =1:min(N_nt,6)
    P(:,i)=P(:,i)-TR_up(:,i)./cellsize^2;
end
if N_nt==7
    P(:,7)=P(:,7)-TR_up(:,7)./cellsize^2.*DOC_fr;
    P(:,8)=P(:,8)-TR_up(:,7)./cellsize^2.*(1-DOC_fr);
end
%%%%%%%%%

%%%% soil concentration update due to ET
TR_avg = ((V_total+(1-ft_up)*EG*dth./(Asur.*Ared)).*(TR_avg))./V_total; % Tracer-1 [C/mm]
TR_avg(V_total==0,:) = 0;


%%%%%%%%%%%%%%
Rd=reshape(Rd,m_cell*n_cell,1); %[mm]
LK=reshape(LK,m_cell*n_cell,1);

[LEAK_nt]= Biogeo_mass_leak(P,LK+Rd,V_total,BiogeoPar); % g/m2
% LEAK_X contains leakage of Rd and leakage to the bottom

rRd=Rd./(LK+Rd); rRd(isnan(rRd))=0; %[-]
Rd_nt=rRd.*LEAK_nt;%[g/m2]
LEAK_nt=LEAK_nt-Rd_nt;%% this is the leadkage to the bottom%[g/m2]

% mass out due to leak and exceed to surface
for i =1:min(N_nt,6)
    P(:,i)=P(:,i)-Rd_nt(:,i)-LEAK_nt(:,i); %[g/m2]
end
if N_nt==7
    P(:,7)=P(:,7)-(Rd_nt(:,7)+LEAK_nt(:,7)).*DOC_fr;
    P(:,8)=P(:,8)-(Rd_nt(:,7)+LEAK_nt(:,7)).*(1-DOC_fr);
end
%%%%%%%%%%%%%%%%%

%%%% calculate the nurtient concentration in the surface
TR_sur= (reshape(Rd.*TR_avg,m_cell,n_cell,N_nt) + Rh.*TR_If)./Qsur ; %%[C/mm]
logic_index = repmat(Qsur == 0, 1, 1, N_nt);
TR_sur(logic_index) = 0;

% TR2_sur=TR_sur(:,:,2); only for test

%%%%%%%% for mass balance check, save the variables before nutrient
%%%%%%%% routing
TR_sur_old=TR_sur;
Qsur_old=Qsur;
TR_If(isnan(TR_If))=0;
TR_If_old= TR_If;

%%%%%%%%% nutrient concentration in the soil, a different shape of TR_avg
TR=zeros(m_cell,n_cell,ms_max,N_nt); % [g/mm]
for i=1:N_nt
    TR(:,:,:,i)=reshape(repmat(TR_avg(:,i),1,ms_max),m_cell,n_cell,ms_max); % assume same concentration in each layer
end
% TR1=TR(:,:,:,1); only for test
% TR2=TR(:,:,:,2);
% TR3=TR(:,:,:,3);




%% concentration for this time step rountin
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
Q_exit = 0; %% Surface Flow exits the domain



Qsub_exit = 0; % Subsurface Flow exits the domain
TR_QoutS=zeros(N_nt,1);
% TR2_QoutS = 0;
% TR3_QoutS = 0;
npoint = length(Xout);
QpointH = zeros(npoint,1);
QpointC = zeros(npoint,1);
UpointH = zeros(npoint,1);
UpointC = zeros(npoint,1);
Qi_outR=zeros(m_cell,n_cell,ms_max);
TR_Qseep=zeros(m_cell,n_cell,N_nt); %% seep to channel [g/mm]


Qi_seep=zeros(m_cell,n_cell,ms_max);
Qi_fall=zeros(m_cell,n_cell);
tH_store = zeros(m_cell,n_cell);
tC_store = zeros(m_cell,n_cell);
%%%%%%%
TR_Qout= zeros(N_nt,1); %% nutrient surface out [g]
TR_ch=zeros(m_cell,n_cell,N_nt);%nutrient channal concentration [g/mm]
TR_QoutC= zeros(N_nt,1); % channel out [g]


TRQexit=zeros(npoint,N_nt); % point exit [g]
TR_QiM=zeros(m_cell,n_cell,ms_max,N_nt); % nutrient mass after subsurface routing [g]
TR_QsurM=zeros(m_cell,n_cell,N_nt); % nutrient mass after surface routing [g]
TR_QchM=zeros(m_cell,n_cell,N_nt); % nutrient mass after channel routing [g]

% TR2_Qout= 0;TR2_ch=zeros(m_cell,n_cell);
% TR2_QoutC=0;TR2Qexit=zeros(npoint,1);


Qi_out_old=Qi_out; % save for later mass balance check

%%% mass out due to previous step subsurface flow, but concentration unchange.

for i =1:min(N_nt,6)
    P(:,i)=max(P(:,i)-reshape(sum(Qi_out_old(:,:,:).*TR(:,:,:,i)/(cellsize^2),3),m_cell*n_cell,1),0);

end
if N_nt==7
    P(:,7)=max(P(:,7)-reshape(sum(Qi_out_old(:,:,:).*TR(:,:,:,7)/(cellsize^2),3),m_cell*n_cell,1).*DOC_fr,0);
    P(:,8)=max(P(:,8)-reshape(sum(Qi_out_old(:,:,:).*TR(:,:,:,7)/(cellsize^2),3),m_cell*n_cell,1).*(1-DOC_fr),0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBSURFACE ROUTING
for jk=1:ms_max
    %%%% No Seepage
    Qi_seep(:,:,jk)=Qi_out(:,:,jk).*(SN); %%[mm]
    Qi_out(:,:,jk)= Qi_out(:,:,jk).*(1-SN); %% [mm]
    %%% SUBSURFACE ROUTING
    Qi_outR(:,:,jk)=Flow_Routing_Step2(DTM,T_potI{jk},Qi_out(:,:,jk)); %% [mm]
    %%%%%% routing for nurtition, only for test
    % [TR2_QiM(:,:,jk)]=Flow_Routing_Step2(DTM,T_potI{jk},Qi_out(:,:,jk).*TR2(:,:,jk)); %%[C] % unit of TR1 is g/mm

    for i = 1: N_nt
        [TR_QiM(:,:,jk,i)]=Flow_Routing_Step2(DTM,T_potI{jk},Qi_out(:,:,jk).*TR(:,:,jk,i)); %%[C] % unit of TR1 is g/mm
        TR_QoutS(i,1) = TR_QoutS(i,1) + (sum(sum(Qi_out(:,:,jk).*TR(:,:,jk,i)))- sum(sum(TR_QiM(:,:,jk,i))));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    Qsub_exit = Qsub_exit + (sum(sum(Qi_out(:,:,jk)))- sum(sum(Qi_outR(:,:,jk)))); %%%  %%[mm]
    % TR2_QoutS = TR2_QoutS + (sum(sum(Qi_out(:,:,jk).*TR2(:,:,jk)))- sum(sum(TR2_QiM(:,:,jk)))); %%%for test 
    Qi_out(:,:,jk)=Qi_outR(:,:,jk);   %% [mm]
end
Qi_in = Qi_out ; %%% [mm] Lateral subsurface for next step

% these codes only rout the mass, it's the same
% for i=1:N_nt
% for jk=1:ms_max
%     Qi_out(:,:,jk)= Qi_out_old(:,:,jk).*(1-SN); %% [mm]
%     %%%%%% routing for nurtition
%     [TR_QiM(:,:,jk,i)]=Flow_Routing_Step2(DTM,T_potI{jk},Qi_out(:,:,jk).*TR(:,:,jk,i)); %%[C] % unit of TR1 is g/mm
%     TR_QoutS(i,1) = TR_QoutS(i,1) + (sum(sum(Qi_out(:,:,jk).*TR(:,:,jk,i)))- sum(sum(TR_QiM(:,:,jk,i)))); %%%  %%[g]
% end
% end
% clear i


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i= 1:N_nt
    TR_Qseep(:,:,i)=sum(Qi_seep.*TR(:,:,:,i),3);
end
% TR2_Qseep= sum(Qi_seep.*TR2,3);


Qi_seep=sum(Qi_seep,3); %% [mm] Seepage flow from soils to channels


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%update P due to subsurface flow, get the back nutrient
%%%%%%%%%
for i =1:min(N_nt,6)
    P(:,i)=P(:,i)+sum(reshape(TR_QiM(:,:,:,i)./(cellsize^2),m_cell*n_cell,ms_max),2);
end
if N_nt==7
    P(:,7)=P(:,7)+sum(reshape(TR_QiM(:,:,:,7)./(cellsize^2),m_cell*n_cell,ms_max),2).*DOC_fr;
    P(:,8)=P(:,8)+sum(reshape(TR_QiM(:,:,:,7)./(cellsize^2),m_cell*n_cell,ms_max),2).*(1-DOC_fr);
end

%%%%%%%%%
MASK=ones(m_cell,n_cell); MASK(isnan(DTM))=0;

TR_Qi_fall=zeros(m_cell,n_cell,N_nt);
for i=1:N_nt
    TR_Qi_fall(:,:,i)=0*MASK; % surface flow fall in channel [g]
end

% 
% TR2_Qi_fall = 0*MASK;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OVERLAND FLOW ROUTING
dti= 60; %%[s] Internal Time step for Surface Routing
cdti = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(sum(Qsur))>0
    while cdti < dt
        % for jj=1:1:dt/dti;
        %%%%%% SURFACE VELOCITY %%%%%%%%%%%
        Y = (Qsur/1000); %% [m]
        %%%%% Ponding - Microroughness
        Y = Y - MRough ; Y(Y<0)=0;
        %%%%
        t= (cellsize.*NMAN_H)./(Y.^(2/3).*sin(atan(Slo_top)).^0.5); %%% [s]
        t(isnan(t))=0;
        tH_store = tH_store + t*dti/dt;
        kdt= dti./t; %[]
        kdt(kdt>1)=1;
        kdt(kdt<0)=1; %% For numerical instabilities
        kdt(isnan(kdt))=0;
        %%%% Surface Routing

        [QsurM]=Flow_Routing_Step2(DTM,T_flow,kdt.*Qsur); %%[mm]
        for i =1:N_nt
            [TR_QsurM(:,:,i)]=Flow_Routing_Step2(DTM,T_flow,kdt.*Qsur.*TR_sur(:,:,i)); %%[mm]
        end
        % [TR2_QsurM]=Flow_Routing_Step2(DTM,T_flow,kdt.*Qsur.*TR2_sur); %%[mm]

        QsurM2 = QsurM.*(1-SN); %%[mm]
        Qi_fall = Qi_fall + QsurM.*(SN); %% [mm]
        I_fall = sum(sum(QsurM.*(SN)));

        TR_QsurM2 = TR_QsurM.*(1-SN); %%[C]
        % TR2_QsurM2 = TR2_QsurM.*(1-SN); %%[C]


        TR_Qi_fall = TR_Qi_fall + TR_QsurM.*(SN); %% [C]
        TR_fall = sum(sum(TR_QsurM.*(SN))); %%[C]
        % TR2_Qi_fall = TR2_Qi_fall + TR2_QsurM.*(SN); %% [C]
        % TR2_fall = sum(sum(TR2_QsurM.*(SN))); %%[C]

        TR_QsurM = TR_QsurM2; %%[C]
        TR_QsurR = TR_QsurM + (TR_sur.*Qsur - (kdt.*Qsur).*TR_sur) ; %%[C]

        % TR2_QsurM = TR2_QsurM2; %%[C]
        % TR2_QsurR = TR2_QsurM + (TR2_sur.*Qsur - (kdt.*Qsur).*TR2_sur) ; %%[C]

        %%%
        QsurM = QsurM2;
        QsurR = QsurM + (Qsur - kdt.*Qsur) ; %% [mm]
        %%%%%%%%%%%%%%
        Q_exit= Q_exit + (sum(sum(Qsur))- sum(sum(QsurR))- I_fall); %% [mm]
        TR_Qout= TR_Qout +  squeeze(  sum(sum(TR_sur.*Qsur)) - sum(sum(TR_QsurR)) - TR_fall  ); %%[C]
        % TR2_Qout= TR2_Qout +  (  sum(sum(TR2_sur.*Qsur)) - sum(sum(TR2_QsurR)) - TR2_fall  ); %%[C]

        for ipo=1:npoint
            QpointH(ipo)= QpointH(ipo) + kdt(Yout(ipo),Xout(ipo))*Qsur(Yout(ipo),Xout(ipo)); %%%%%[mm]
        end
        %%%
        Qsur=QsurR; %%[mm]
        TR_sur = TR_QsurR./Qsur;  %%[C/mm]
        logic_index = repmat(Qsur == 0, 1, 1, N_nt);
        TR_sur(logic_index) = 0;
        % TR_sur(Qsur==0)=0;

        % TR2_sur = TR2_QsurR./Qsur;  %%[C/mm]
        % TR2_sur(Qsur==0)=0;
        cdti = cdti +dti;
        dti =  min(min(t(t>0)));
        if cdti+dti>dt; dti = dt-cdti; end
    end
else
    Qsur = zeros(m_cell,n_cell);

    TR_sur =zeros(m_cell,n_cell,N_nt); %%[C/mm]
    TR_QsurR=zeros(m_cell,n_cell,N_nt); %%[C]
    TR_Qi_fall =zeros(m_cell,n_cell,N_nt);  %% [C]
    % 
    % TR2_sur =zeros(m_cell,n_cell); %%[C/mm]
    % TR2_Qi_fall =0;  %% [C]
end
%%%%%%%%%%
q_runon = Qsur;
%%% TR_sur is already updated within the routing..


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHANNEL FLOW ROUTING, TR1_Qseep [g] already multiply TR1
q_channel = q_channel_in + Qi_seep + Qi_fall; %%% [mm]
TR_ch = ( q_channel_in.*(TR_ch) + TR_Qseep +  TR_Qi_fall)./q_channel ; %[C/mm]

logic_index = repmat(q_channel == 0, 1, 1, N_nt);
TR_ch(logic_index) = 0;

% TR_ch(q_channel==0)=0;

% TR2_ch = ( q_channel_in.*(TR2_ch) + TR2_Qseep +  TR2_Qi_fall)./q_channel ; %[C/mm]
% TR2_ch(q_channel==0)=0;


dti= 2; %%[s] Internal Time step for Surface Routing
cdti = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(sum(q_channel))>0
    % for jj=1:1:dt/dti;
    while cdti < dt
        %%%%%% SURFACE VELOCITY %%%%%%%%%%%
        Y = (q_channel/1000).*cellsize./WC ; %% [m]
        t= (cellsize.*NMAN_C)./(Y.^(2/3).*sin(atan(Slo_top)).^0.5); %%% [s]
        t(isnan(t))=0;
        tC_store = tC_store + t*dti/dt;
        kdt= dti./t; %[]
        kdt(kdt>1)=1;
        kdt(kdt<0)=1; %% For numerical instabilities
        kdt(isnan(kdt))=0;
        %%%% Surface Routing
        [QchM]=Flow_Routing_Step2(DTM,T_flow,kdt.*q_channel); %%[mm]
        for i =1:N_nt
            [TR_QchM(:,:,i)]=Flow_Routing_Step2(DTM,T_flow,(kdt.*q_channel).*TR_ch(:,:,i)); %%[C]
        end
        % [TR2_QchM]=Flow_Routing_Step2(DTM,T_flow,(kdt.*q_channel).*TR2_ch); %%[C]

        QchR = QchM + (q_channel - kdt.*q_channel) ; %% [mm]
        TR_QchR = TR_QchM + (TR_ch.*q_channel - (kdt.*q_channel).*TR_ch) ; %%[C]
        % TR2_QchR = TR2_QchM + (TR2_ch.*q_channel - (kdt.*q_channel).*TR2_ch) ; %%[C]

        %%%%%%%%%%%%%%
        Q_exit= Q_exit + (sum(sum(q_channel))- sum(sum(QchR))); %% [mm]
        TR_QoutC= TR_QoutC +  squeeze((  sum(sum(TR_ch.*q_channel)) - sum(sum(TR_QchR)) )); %%[C]
        % TR2_QoutC= TR2_QoutC +  (  sum(sum(TR2_ch.*q_channel)) - sum(sum(TR2_QchR)) ); %%[C]


        for ipo=1:npoint
            QpointC(ipo)= QpointC(ipo) + kdt(Yout(ipo),Xout(ipo))*q_channel(Yout(ipo),Xout(ipo)); %%%%%[mm]
            for i=1:N_nt
                TRQexit(ipo,i) = TRQexit(ipo,i) + kdt(Yout(ipo),Xout(ipo))*TR_ch(Yout(ipo),Xout(ipo),i).*q_channel(Yout(ipo),Xout(ipo)); %[C]
            end
            % TR2Qexit(ipo) = TR2Qexit(ipo) + kdt(Yout(ipo),Xout(ipo))*TR2_ch(Yout(ipo),Xout(ipo)).*q_channel(Yout(ipo),Xout(ipo)); %[C]

        end
        %%%%%%%%%%%%%%%%
        q_channel=QchR; %%[mm]
        TR_ch = TR_QchR./q_channel;  %%[C/mm]
        logic_index = repmat(q_channel == 0, 1, 1, N_nt);
        TR_ch(logic_index) = 0;
        %TR_ch(q_channel==0)=0;
        % TR2_ch = TR2_QchR./q_channel;  %%[C/mm]
        % TR2_ch(q_channel==0)=0;
        cdti = cdti +dti;
        dti =  min(min(t(t>0)));
        if cdti+dti>dt; dti = dt-cdti; end
    end
else
    q_channel = zeros(m_cell,n_cell);
    TR_ch =zeros(m_cell,n_cell,N_nt); %%[C/mm]
    TR_QoutC =zeros(N_nt,1);  %% [C] % I think out means out of the domain, can go out in many cells, but exit means we only track the outlet point?
    TRQexit = 0;
    % TR2_ch =zeros(m_cell,n_cell); %%[C/mm]
    % TR2_QoutC =0;  %% [C]
    % TR2Qexit = 0;
end
q_channel_out = q_channel;%[mm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% mass balance for subsurface routing
% mass out by subsurface routing from last time step - left mass(back to P) this time step - out domain mass
CK_Sub=(squeeze(sum(sum(sum(Qi_out_old.*TR,3))))-squeeze(sum(sum(sum(TR_QiM))))-squeeze(sum(sum(TR_Qseep)))-TR_QoutS)./squeeze(sum(sum(sum(Qi_out_old.*TR,3))));
%%%%%%%%%% mass balance for surface routing. original, out, left
CK_S=(squeeze(sum(sum(Qsur_old.*TR_sur_old)))-TR_Qout-squeeze(sum(sum(TR_QsurR)))-squeeze(sum(sum(TR_Qi_fall))))./squeeze(sum(sum(Qsur.*TR_sur)));
%%%% mass balance for each nutrient

%%%% 
% old mass- leak - evaporation - exceed to surface + infiltration - subsurface out - seep to channel= new mass 
%%%%%%
if N_nt <7
    N_check=(sum(P_old(:,1:N_nt))-sum(P(:,1:N_nt))-sum(LEAK_nt(:,1:N_nt))-sum(TR_up./cellsize^2)-sum(Rd_nt(:,1:N_nt))...
        +sum(If.*reshape(TR_If_old/cellsize^2,m_cell*n_cell,N_nt))...
        -reshape(TR_QoutS/cellsize^2,1,[])...
        -sum(reshape(TR_Qseep./(cellsize^2),m_cell*n_cell,N_nt)))...
        ./sum(P_old(:,1:N_nt))
else
    N_check(1:6)=(sum(P_old(:,1:6))-sum(P(:,1:6))-sum(LEAK_nt(:,1:6))-sum(TR_up(:,1:6)./cellsize^2)-sum(Rd_nt(:,1:6))...
        +sum(If.*reshape(TR_If_old(:,:,1:6)/cellsize^2,m_cell*n_cell,6))...
        -reshape(TR_QoutS(1:6)/cellsize^2,1,[])...
        -sum(reshape(TR_Qseep(:,:,1:6)./(cellsize^2),m_cell*n_cell,6)))...
        ./sum(P_old(:,1:6));
    N_check(7)=(sum(sum(P_old(:,7:8)))-sum(sum(P(:,7:8)))-sum(LEAK_nt(:,7))-sum(TR_up(:,7)./cellsize^2)-sum(Rd_nt(:,7))...
        +sum(If.*reshape(TR_If_old(:,:,7)/cellsize^2,m_cell*n_cell,1))...
        -reshape(TR_QoutS(7)/cellsize^2,1,[])...
        -sum(reshape(TR_Qseep(:,:,7)./(cellsize^2),m_cell*n_cell,1)))...
        ./sum(sum(P_old(:,7:8)));
end

%%%%%%%%%%%%%%
for ipo=1:npoint
    UpointH(ipo)= cellsize./tH_store(Yout(ipo),Xout(ipo)); %%% Surface Velocity [m/s]
    UpointC(ipo)= cellsize./tC_store(Yout(ipo),Xout(ipo)); %%% Surface Velocity [m/s]
end
Utot_H= cellsize./tH_store; %%% Surface Velocity [m/s]
Utot_C= cellsize./tC_store; %%% Surface Velocity [m/s]
%%%%%%%%%%%%%%%%%%%%%%%%
Q_exit = Q_exit*(cellsize^2)/Area; %%[mm] includes surface and channel
Qsub_exit = Qsub_exit*(cellsize^2)/Area; %%[mm] subsurface

TR_exit = (TR_Qout+TR_QoutC); % [g] surface + channel
TR_exitS=TR_QoutS; %[g] subsurface



%%%%%%%%%%%%%%%%%%%%%%%%
%QpointH = (QpointH/dth)*(cellsize^2)/(3600000); %% [m^3/s]
%QpointC = (QpointC/dth)*(cellsize^2)/(3600000); %% [m^3/s]
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%%%%%%% HYDRAULIC HEAD PART %%%
T_pot=cell(1,ms_max);
Slo_pot = zeros(m_cell,n_cell,ms_max);
met_fdir=1;
if OPT_HEAD == 1
    %%%% NEW FLOW DIRECTIONS %%%%
    %%%% Estimation of Energy Slopes
    for jk=1:ms_max
        H = DTM + 0.001*reshape(POT(:,jk),m_cell,n_cell).*cos(atan(Slo_top)); %%% Hydraulic Head [m]
        %H = 0.001*reshape(POT(:,jk),m_cell,n_cell).*cos(atan(Slo_top)); %%% Water Potential Head [m]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if met_fdir == 1
            %%% D-Inf flow direction method
            [R_pot, Slo_pot(:,:,jk)] = dem_flow(H,cellsize,cellsize);
            T_pot{jk} = flow_matrix(H,R_pot,cellsize,cellsize); %% Flow Matrix %%%%
        else
            Slo_pot(:,:,jk)=Slope_Aspect_indexes(H,cellsize,'mste');
            x=0:cellsize:(0+cellsize*(n_cell-1));
            y=0:cellsize:(0+cellsize*(m_cell-1));
            [x,y]=meshgrid(x,y);
            [Mpot] = flowdir(x,y,H,'type','multi'); %% Multiple D-Inf Quinn et al., 1993
            %[Mpot] = flowdir(x,y,H,'type','single'); %% D8  O'Callaghan & Mark, 1984
            T_pot{jk}=speye(m_cell*n_cell,m_cell*n_cell)-Mpot';
        end
    end
    Slo_pot(Slo_pot<0)=0;
else
    for jk=1:ms_max;
        T_pot{jk}= T_flow;
        Slo_pot(:,:,jk)=Slo_top;  %%%
    end
end
return


