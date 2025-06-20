function[dB_flux,P,BLit,Bam,Bem,Nuptake_H,Puptake_H,Kuptake_H,Nuptake_L,Puptake_L,Kuptake_L,RexmyI,...
    R_litter,R_microbe,R_litter_sur,R_ew,VOL,N2flx,Min_N,Min_P,R_bacteria,RmycAM,RmycEM,Prod_B,Prod_F,BfixN,NavlI,LitFirEmi]= Biogeo_MODULE_PAR(ISOIL_L,ISOIL_H,Rexmy_L,Rexmy_H,ManIH,ManIL,Ptm1,ZBIOG,rsd,PH,Ts,Ta,Psi_s,Se,Se_fc,V,VT,Ccrown,Bio_Zs,RfH_Zs,RfL_Zs,...
    Lk,T_H,T_L,Broot_H,Broot_L,LAI_H,LAI_L,...
    SupN_H,SupP_H,SupK_H,SupN_L,SupP_L,SupK_L,RexmyI,ExEM,NavlI,Pcla,Psan,...
    B_IO,jDay,AAET,cc_max)
dB_flux=zeros(55,8);
P=zeros(1,55);
Nuptake_H= zeros(1,cc_max);
Puptake_H= zeros(1,cc_max);
Kuptake_H= zeros(1,cc_max);
Nuptake_L= zeros(1,cc_max);
Puptake_L= zeros(1,cc_max);
Kuptake_L= zeros(1,cc_max);
Bam=zeros(1,1);
Bem=zeros(1,1);
R_litter=zeros(1,1);
R_microbe=zeros(1,1);
R_litter_sur=zeros(1,1);
R_ew=zeros(1,1);
VOL=zeros(1,1);
N2flx=zeros(1,1);
Min_N=zeros(1,1);
Min_P=zeros(1,1);
R_bacteria=zeros(1,1);
RmycAM=zeros(1,1);
RmycEM=zeros(1,1);
Prod_B=zeros(1,1);
Prod_F=zeros(1,1);
BfixN=zeros(1,1);
LitFirEmi=zeros(1,2);



IS= Ccrown*squeeze(ISOIL_L) + Ccrown*squeeze(ISOIL_H); % Plant exports to Litter and soil Low/High Vegetation Composed of 18 different components of Export Types/Element
Rexmy= Ccrown*squeeze(Rexmy_L) + Ccrown*squeeze(Rexmy_H); % Root exudates (1),export of carbon toward mycorrhiza (2) and export to root-noduli (3), High vegetation
FireA = max(max(ManIH),max(ManIL)); % management indicator, also FireA later on is related to it


[dB_flux,P,Nuptake_H,Puptake_H,Kuptake_H,Nuptake_L,Puptake_L,Kuptake_L,RexmyI,...
    R_litter,R_microbe,R_litter_sur,R_ew,VOL,N2flx,Min_N,Min_P,R_bacteria,RmycAM,RmycEM,Prod_B,Prod_F,BfixN,NavlI,LitFirEmi]= BIOGEO_UNIT_0Leak(Ptm1,IS,ZBIOG,rsd,PH,Ts,Ta,Psi_s,Se,Se_fc,V,VT,Ccrown,Bio_Zs,RfH_Zs,RfL_Zs,...
    Lk,T_H,T_L,Broot_H,Broot_L,LAI_H,LAI_L,...
    SupN_H,SupP_H,SupK_H,SupN_L,SupP_L,SupK_L,Rexmy,RexmyI,ExEM,NavlI,Pcla,Psan,...
    B_IO,jDay,FireA,AAET);

BLit=0.002*sum(P(1:5))*Ccrown; %% %%[kg DM / m2] Total litter content on the surface, dry matter
Bam =  P(20); %%[gC/m2] % %%% B20 AM-Mycorrhizal - C
Bem =  P(21); %%[gC/m2] %%% B21 EM-Mycorrhizal - C


return