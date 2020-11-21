function [KOdata,KO_lb,KO_ub]=parse_KO_data(model_irrev,model_rev,replica,media,biomass_choose,KOdata_file)
    if strcmp(biomass_choose,'wt')
        biom=find(contains(model_irrev.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M'));
    elseif strcmp(biomass_choose,'core')
        biom=find(contains(model_irrev.rxns,'BIOMASS_Ec_iJO1366_core_53p95M'));
    end
    
    name=append(media.carbonreac,'_f');
    [~,cs_ind]=ismember(media.carbonreac,model_rev.rxns);
    %Part one: create a abudance matrix (rows=reactions, columns=strains)
    KO_reactions_u=unique(KOdata_file.reac(:,1));
    ind_b=find(contains(KO_reactions_u,'_b'));
    reac_names=KO_reactions_u;
    reac_names(ind_b)=erase(KO_reactions_u(ind_b),'_b');
    [KO_reactions_u2,~,~]=unique(reac_names); 
    %%These reactions are the ones we need to create the protein abundance
    
    [~,ind_model]=ismember(KO_reactions_u2,model_rev.rxns);
    ind_model=ind_model(ind_model~=0);
    [~,back_model]=ismember(KO_reactions_u(ind_b),model_irrev.rxns);
    back_model=back_model(back_model~=0);
    reacind_model=[ind_model;back_model];
  
    % '_b' appended for backward reactions, used to map the reactions to
    % model
    reac_model=[model_rev.rxns(ind_model);model_irrev.rxns(back_model)];
    
    n=length(reacind_model);
    KOdata.abun=NaN(n,18);
    KOdata.reac=model_irrev.rxns(reacind_model);
    KOdata.reac_orig=reac_model;
    KOdata.reacind=reacind_model;
    KOdata.genes=model_irrev.grRules(reacind_model);
	KOdata.flux_av=NaN(n,18);
    KOdata.flux_ub=NaN(n,18);
    KOdata.flux_lb=NaN(n,18);
    KOdata.kapp=NaN(n,18);


    for strain=1:18
        %index of reactions of the strain within the excel file
        ind=intersect(find(contains(KOdata_file.sample_id,KOdata_file.strains(strain))),find(contains(KOdata_file.sample_rep,replica)));

        [~,ind_inabun]=ismember(KOdata_file.reac(ind,1),reac_model);
        %[~,ind_inabun]=ismember(KOdata_file.bnumd(ind),model_irrev.genes);
        ind=ind(ind_inabun~=0);
        %index of reactions of the strain within the abundance matrix
        ind_inabun=ind_inabun(ind_inabun~=0);
        
        %abundance given in fmol/gDW. Only considered if it is >=50pmol
        abun=KOdata_file.abun;
        abun(abun<50000)=NaN;
        KOdata.abun(ind_inabun,strain)=abun(ind)*1e-12;
        
        
        KOdata.flux_av(ind_inabun,strain)=KOdata_file.flux_ave(ind);
        KOdata.flux_ub(ind_inabun,strain)=KOdata_file.flux_ub(ind);
        KOdata.flux_lb(ind_inabun,strain)=KOdata_file.flux_lb(ind);
        KOdata.kapp(ind_inabun,strain)=KOdata_file.kapp(ind);
        %ind=intersect(find(~isnan(KOdata.flux_av(:,strain))),find(isnan(KOdata.kapp(:,strain))));
        %KOdata.flux_av(ind,strain)=0;
    end
    
    
    
    %Part2: create lb and ub matrices

    n_r=size(model_irrev.S,2);
    KO_lb=zeros(n_r,18);
    KO_ub=zeros(n_r,18);
    for i=1:18
        KO_lb(:,i)=model_irrev.lb;
        KO_ub(:,i)=model_irrev.ub;
    
        KO_genes=KOdata_file.strain_mutations{i};
        [~,gene_ind]=ismember(KO_genes,model_irrev.proteins);
        gene_ind=gene_ind(gene_ind~=0);

        reac_ind=find(contains(model_irrev.grRules,model_irrev.genes(gene_ind)));

        %for every reaction containing at least one of theses genes
        for j=1:length(reac_ind)
            check=split(model_irrev.grRules(reac_ind(j)),' or ');
            check2=split(model_irrev.grRules(reac_ind(j)),' and ');
            if length(check)>1 && length(check2)==1 %isoenzyme

                [~,ind2]=ismember(check,model_irrev.genes(gene_ind)); %check if the genes belong to the KO set
                if sum(ind2~=0)==length(check) %all possible isoenzymes are KO
                    KO_ub(reac_ind(j),i)=0;
                end

            elseif length(check2)>1 && length(check)==1 %only having 'and'
                KO_ub(reac_ind(j),i)=0;

            elseif length(check)==1 && length(check2)==1 %single genes
                KO_ub(reac_ind(j),i)=0;

            elseif length(check)>1 && length(check2)>1 %complex (&)or(&)
                count=0;
                for m=1:length(check)
                    check3=split(check(m),' and ');
                    check3=erase(check3,'(');
                    check3=erase(check3,')');

                    [~,ind]=ismember(check3,model_irrev.genes(gene_ind));
                    if sum(ind~=0)>0 
                        count=count+1;
                    end
                end
                if count==length(check)
                    KO_ub(reac_ind(j),i)=0;
                end
            end
        end
        
        

        %model_irrev=deleteModelGenes(model_rev,model_rev.genes(dele_genes),0);
        

        if KOdata_file.gro_rate(i,2)~=0
            KO_lb(biom,i)=KOdata_file.gro_rate(i,1);
            KO_ub(biom,i)=KOdata_file.gro_rate(i,2);
        end

        ind_glu=find(contains(model_irrev.rxns,'EX_glc__D_e_f'));
        if KOdata_file.glc_up(i,2)~=0
            KO_lb(model_irrev.match(ind_glu),i)=KOdata_file.glc_up(i,1);
            KO_ub(model_irrev.match(ind_glu),i)=KOdata_file.glc_up(i,2);
        end

        limited_f=setdiff(cs_ind,ind_glu); %limit other carbon sources
        limited_b=model_irrev.match(limited_f);

        KO_ub(limited_b,i)=0;
        ind_ace=find(contains(model_irrev.rxns,'EX_ac_e_f'));

        if KOdata_file.ace_se(i,2)~=0 % there were extra info in labeling experiment
            KO_lb(ind_ace,i)=KOdata_file.ace_se(i,1);
            KO_ub(ind_ace,i)=KOdata_file.ace_se(i,2); 
        end
        
        ind_suc=find(contains(model_irrev.rxns,'EX_succ_e_f'));
        if KOdata_file.suc_se(i,1)~=0 
            KO_lb(ind_suc,i)=KOdata_file.suc_se(i,1);
            KO_ub(ind_suc,i)=KOdata_file.suc_se(i,2);
        end

%         if KOdata_file.lac_se(i,1)~=0 
%             KO_lb(95,i)=KOdata_file.lac_se(i,1);
%             KO_ub(95,i)=KOdata_file.lac_se(i,2);
%         end
    end
end