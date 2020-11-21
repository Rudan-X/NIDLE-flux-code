function [abundance_31]=from58to31(model,abundance_file,conditions31)

[~,cond_ind]=ismember(conditions31.cond,abundance_file.cond);

%one condition not found
%conditions31(cond_ind==0)
%check=find(contains(data31_file.cond,'GLC_CHEM_mu=0.51'));
%abundance_file.cond(check) %there is a replicate for GLC_CHEM_mu=0.51, take one of them
cond_ind(cond_ind==0)=18;
abundance_31.cond=conditions31;

[~,ind]=ismember(abundance_file.genes,model.genes);
ind2=find(ind~=0);
abundance_31.abun=abundance_file.abun(ind2,cond_ind); %take only the 31 out of 58
abundance_31.genes=abundance_file.genes(ind2);
end
