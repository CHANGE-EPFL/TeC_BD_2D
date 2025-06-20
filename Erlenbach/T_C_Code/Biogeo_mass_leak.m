%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction BIOGEOCHEMISTRY_DYNAMIC   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[LEAK_nt]= Biogeo_mass_leak(B,Lk,Vtotal,BiogeoPar)  
%%%%%%%%%%%%%% Leakage of Nitrogen Compound 
%%%% Proportional to leakage from the nitrogen zone with fraction of
%%%% dissolved ammonium and nitrate 
%%% For simplicity leakage only at the base of soil column 
%%%
VSUM=Vtotal;


%%%% Porporato et al., 2003 
aNH4 =BiogeoPar.aNH4;  
aNO3 = BiogeoPar.aNO3;  
aP = BiogeoPar.aP;
aK = BiogeoPar.aK; 
aDON= BiogeoPar.aDON;
aDOP= BiogeoPar.aDOP;
aDOC = BiogeoPar.aDOC; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% we are calculating the mass
% here we may have a unit problem?? this is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actually [gN/ gsoil], should *(ZBIOG*rsd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to get g/m2
LEAK_NH4 = max(0,min(B(:,1),aNH4*B(:,1).*(Lk)./(VSUM))); %% [gN/m2]
LEAK_NO3 = max(0,min(B(:,2),aNO3*B(:,2).*(Lk)./(VSUM))); %% [gN/m2]
LEAK_P = max(0,min(B(:,3),aP*B(:,3).*(Lk)./(VSUM))); %% [gP/m2]
LEAK_K = max(0,min(B(:,4),aK*B(:,4).*(Lk)./(VSUM))); %% [gK/m2]
LEAK_DON = max(0,min(B(:,5),aDON*B(:,5).*(Lk)./(VSUM))); %% [gN/m2]
LEAK_DOP = max(0,min(B(:,6),aDOP*B(:,6).*(Lk)./(VSUM))); %% [gP/m2]
LEAK_DOC = max(0,min(B(:,7)+B(:,8),aDOC*(B(:,7)+B(:,8)).*(Lk)./(VSUM))); %% [gC/m2]
% %%%%%%%%%%%%%%%%%%
LEAK_nt=[LEAK_NH4,LEAK_NO3,LEAK_P,LEAK_K,LEAK_DON,LEAK_DOP,LEAK_DOC];


return 
