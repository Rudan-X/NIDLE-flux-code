%%
%Solver setup
changeCobraSolver('gurobi', 'MILP');

%%
%reconstruction model
load('../data/iJO1366.mat');

%31 growth conditions and 18 KO strainsc
load('../data/conditions31');
load('../data/data18_file');
load('../data/abundance_file');
%%
%Parse the experimental data
[model_irrev,model_rev,lb,ub]=parse_bounds(iJO1366,conditions31,'core');
%model_irrev: uptake reactions as reversible reactions and then converted
%into irreversible model
%model_rev: model with reactions in both directions
%lb and ub matrix: rows(reactions) x columns(conditions)

[KOdata_B1,KO_lb,KO_ub]=parse_KO_data(model_irrev,model_rev,'B1',conditions31,'core',KOdata_file);
[KOdata_B2,~,~]=parse_KO_data(model_irrev,model_rev,'B2',conditions31,'core',KOdata_file);

%Protein abundance mapped to genes
abundance_31=from58to31(model_irrev,abundance_file,conditions31);
toremove=find(all(isnan(abundance_31.abun),2));
abundance_31.abun(toremove,:)=[];
abundance_31.genes(toremove,:)=[]; 

[genes_u,ind_u]=unique(KOdata_B1.genes);
abundance_18A.genes=string(genes_u);
abundance_18A.abun=KOdata_B1.abun(ind_u,:);
[genes_u,ind_u]=unique(KOdata_B2.genes);
abundance_18B.genes=string(genes_u);
abundance_18B.abun=KOdata_B2.abun(ind_u,:);


%Protein abundance mapped to reactions
abundance=parse_abundance(model_irrev, conditions31,abundance_file);

%Get GPR states
[abun_mapped,g_vect_31_comp]=parse_rules(model_irrev,abundance_31);
[abun_mapped_18A,g_vect_18A]=parse_rules(model_irrev,abundance_18A);
[abun_mapped_18B,g_vect_18B]=parse_rules(model_irrev,abundance_18B);

%%
%Running pFBA and NIDLE-flux
scale=1e6; %the resulting flux will be 1e-6 times the original scale
conv=1e-6;
%threshold for pFBA to consider non-zero: 1e-10

res_31_p=struct();
res_31_p.flux=pFBA(model_irrev,lb,ub,scale,'core'); 
res_31_e10_comp=NIDLE(model_irrev,g_vect_31_comp,lb,ub,1e-4,scale,'core');


res_18_p=struct();
res_18_p.flux=pFBA(model_irrev,KO_lb,KO_ub,scale,'core');
res_18A_e10=NIDLE(model_irrev,g_vect_18A,KO_lb,KO_ub,1e-4,scale,'core');
res_18B_e10=NIDLE(model_irrev,g_vect_18B,KO_lb,KO_ub,1e-4,scale,'core');

%%
%Computing Kapp
%_p31 refers to pFBA at 31 conditions
%_n18A referst to NIDLE-flux at replicate 1 on 18 strains
%count vector stores statistics of the results

[Kapp_p31,V_p31,count_p31]=getkapp_pfba(abundance,g_vect_31_comp, res_31_p.flux,conv,1e-4+1e-5);
[Kapp_n31,V_n31,count_n31]=getkapp(abundance,g_vect_31_comp,res_31_e10_comp.flux,conv,1e-4+1e-5);

[Kapp_p18A,V_p18,count_p18A]=getkapp_pfba(KOdata_B1,g_vect_18A,res_18_p.flux,conv,1e-4+1e-5);
[Kapp_p18B,~,count_p18B]=getkapp_pfba(KOdata_B2,g_vect_18B,res_18_p.flux,conv,1e-4+1e-5);

[Kapp_n18A,V_n18A,count_n18A]=getkapp(KOdata_B1,g_vect_18A,res_18A_e10.flux,conv,1e-4+1e-5);
[Kapp_n18B,V_n18B,count_n18B]=getkapp(KOdata_B2,g_vect_18B,res_18B_e10.flux,conv,1e-4+1e-5);

%Compute k_max
[kmax_p31,kmax_p31_vect]=getkmax(Kapp_p31,V_p31(abundance.reacind,:),abundance);
[kmax_n31,kmax_n31_vect]=getkmax(Kapp_n31,V_n31(abundance.reacind,:),abundance);

[kmax_n18]=getkmax_av(Kapp_n18A,Kapp_n18B,KOdata_B1);
[kmax_p18]=getkmax_av(Kapp_p18A,Kapp_p18B,KOdata_B1);

%Find idle enzymes
idle_p31=idle_enzyme(model_irrev,g_vect_31_comp,V_p31,abundance_31,1e-10+1e-11);
idle_n31=idle_enzyme(model_irrev,g_vect_31_comp,V_n31,abundance_31,1e-10+1e-11);

idle_p18A=idle_enzyme(model_irrev,g_vect_18A,V_p18,abundance_18A,1e-10+1e-11);
idle_p18B=idle_enzyme(model_irrev,g_vect_18B,V_p18,abundance_18B,1e-10+1e-11);

idle_n18A=idle_enzyme(model_irrev,g_vect_18A,V_n18A,abundance_18A,1e-10+1e-11);
idle_n18B=idle_enzyme(model_irrev,g_vect_18B,V_n18B,abundance_18A,1e-10+1e-11);




