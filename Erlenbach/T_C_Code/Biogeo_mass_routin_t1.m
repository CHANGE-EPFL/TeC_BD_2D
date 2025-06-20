%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction BIOGEOCHEMISTRY_DYNAMIC   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[LEAK_NH4,LEAK_NO3,LEAK_P,LEAK_K,LEAK_DOC,LEAK_DON,LEAK_DOP]= Biogeo_mass_routin(B,Lk,V,Routing_out_volume,BiogeoPar)  
%%%%%%%%%%%%%% Leakage of Nitrogen Compound 
%%%% Proportional to leakage from the nitrogen zone with fraction of
%%%% dissolved ammonium and nitrate 
%%% For simplicity leakage only at the base of soil column 
%%%
VSUM=sum(V,2)+Routing_out_volume+Lk;


%%%% Porporato et al., 2003 
aNH4 =BiogeoPar.aNH4;  
aNO3 = BiogeoPar.aNO3;  
aP = BiogeoPar.aP;
aK = BiogeoPar.aK; 
aDOC = BiogeoPar.aDOC; 
aDON= BiogeoPar.aDON;
aDOP= BiogeoPar.aDOP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% we are calculating the mass
% here we may have a unit problem?? this is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actually [gN/ gsoil], should *(ZBIOG*rsd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to get g/m2
LEAK_NH4 = max(0,min(B(:,31),aNH4*B(:,31).*(Lk)./(VSUM))); %% [gN/m2]
LEAK_NO3 = max(0,min(B(:,32),aNO3*B(:,32).*(Lk)./(VSUM))); %% [gN/m2]
LEAK_P = max(0,min(B(:,43),aP*B(:,43).*(Lk)./(VSUM))); %% [gP/m2]
LEAK_K = max(0,min(B(:,52),aK*B(:,52).*(Lk)./(VSUM))); %% [gK/m2]
LEAK_DOC = max(0,min(B(:,12)+B(:,13),aDOC*(B(:,12)+B(:,13)).*(Lk)./(VSUM))); %% [gC/m2]
LEAK_DON = max(0,min(B(:,33),aDON*B(:,33).*(Lk)./(VSUM))); %% [gN/m2]
LEAK_DOP = max(0,min(B(:,47),aDOP*B(:,47).*(Lk)./(VSUM))); %% [gP/m2]
% %%%%%%%%%%%%%%%%%%


return 
