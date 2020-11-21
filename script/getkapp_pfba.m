function [Kapp_matrix,V_matrix,count]=getkapp_pfba(abundance_mapp,g_vect,sol_flux,conv,thre)
R_index=abundance_mapp.reacind;
R_ab_matrix=abundance_mapp.abun;
R_names=abundance_mapp.genes;
n=size(R_ab_matrix,2);
Kapp_matrix=zeros(length(R_index),n);
V_matrix=zeros(size(sol_flux,1),n);

%number of reactions with flux larger than the threshold over the whole
%model
count.nonzero=zeros(n,1); 

%number of reactions with flux larger than the threshold with mapped
%enzymes
count.nonzero2=zeros(n,1); 

%number of reactions with flux larger than the threshold and with non-nan
%abundance
count.kapp=zeros(n,1);


%number of reactions with non-nan abundance inside mapped one
count.withabun=zeros(n,1);

%number of unique enzymes catalyzing single reactions
count.homoenzyme=zeros(n,1);

%number of active reactions with active gprRules versus total reactions
%with active rule
count.ratio=zeros(n,1);

for cond=1:n
    V=sol_flux(:,cond);
    V(V<=thre)=0;
    V(V<0)=0;
    V_matrix(:,cond)=V*conv;
    const=find(g_vect(:,cond)==1);
    Kapp_matrix(:,cond)=V(R_index)./R_ab_matrix(:,cond)*conv/3600;
    count.nonzero(cond)=sum(V>thre);
    count.nonzero2(cond)=sum(V(abundance_mapp.reacind)>thre);
    count.kapp(cond)=length(intersect(find(Kapp_matrix(:,cond)~=0),find(~isnan(Kapp_matrix(:,cond)))));
    count.withabun(cond)=sum(~isnan(R_ab_matrix(:,cond)));
    count.homoenzyme(cond)=length(unique(R_names(~isnan(R_ab_matrix(:,cond)))));
    count.ratio(cond)=sum(V(const)>thre)/sum(g_vect(:,cond)~=0);
end
end