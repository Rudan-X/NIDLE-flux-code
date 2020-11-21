function [abundance_31]=parse_abundance(model_irrev,conditions31,abundance_file) 
    [~,cond_ind]=ismember(conditions31.cond,abundance_file.cond);

    %one condition not found
    %conditions31(cond_ind==0)
    %check=find(contains(abundance_file.cond,'GLC_CHEM_mu=0.51'));
    %abundance_file.cond(check) %there is a replicate for GLC_CHEM_mu=0.51, take one of them

    cond_ind(cond_ind==0)=18;

    abundance_31.cond=conditions31;
    abundance_31.abun=abundance_file.abun(:,cond_ind); %take only the 31 out of 58
    abundance_31.genes=abundance_file.genes;


    n_r=size(model_irrev.S,2);
    rules_type=zeros(n_r,1);

    for i=1:n_r
        if ~isempty(model_irrev.grRules{i})
            if isempty(find(contains(model_irrev.grRules{i},'or'))) && isempty(find(contains(model_irrev.grRules{i},'and')))
                rules_type(i)=1;
            end
        end
    end

    %homomeric enzymes in model
    homomeric_en_model=model_irrev.grRules(rules_type==1);
    homomeric_ind_model=find(rules_type==1);
    homomeric_reac_model=model_irrev.rxns(rules_type==1);

    %find these enzymes in the abundance file
    [~,enzyme_ind]=ismember(homomeric_en_model,abundance_file.genes);
    
    
    %The reactions we need to extract from the model
    abundance_31.reacind=homomeric_ind_model(enzyme_ind~=0);
    abundance_31.reac=string(homomeric_reac_model(enzyme_ind~=0));

    %The genes we need to extract from the abundance file
    enzyme_ind_found=enzyme_ind(enzyme_ind~=0);
    abundance_31.abun=abundance_31.abun(enzyme_ind_found,:);
    abundance_31.genes=abundance_31.genes(enzyme_ind_found);
    

    %there are some genes with NaN value in all 31 conditions
    tokeep=find(~all(isnan(abundance_31.abun),2));
    abundance_31.abun=abundance_31.abun(tokeep,:);
    abundance_31.genes=abundance_31.genes(tokeep);
    abundance_31.reacind=abundance_31.reacind(tokeep);
    abundance_31.reac=abundance_31.reac(tokeep);

    %save('abundance_31.mat','abundance_31')
end