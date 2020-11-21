function [media_model_irrev,media_model,lb,ub]=parse_bounds(iJO1366,media,biomass_choose)


%media_model = removeRxns(iJO1366, {'E. coli biomass objective function (iJO1366) - core - with 53.95 GAM estimate'});
media_model=iJO1366;

%remove the s0001 grRules for the model. There are 3 types
%'s0001', 's0001 or ...' and 'bxxx or s0001'

excep=find(contains(media_model.grRules,'s0001'));
copy=media_model.grRules(excep);
copy=erase(copy,'s0001 or ');
copy=erase(copy,' or s0001');
copy=erase(copy,'s0001');
media_model.grRules(excep)=copy;

%remove the core biomass reaction for the rest of analysis
%14 WT, 19 core
if strcmp(biomass_choose,'wt')
    media_model=removeRxns(media_model,media_model.rxns(19));
elseif strcmp(biomass_choose,'core')
    media_model=removeRxns(media_model,media_model.rxns(14)); 
end


%Since by default the iJO1366 has only glucose exchange as reversible
%reactions, other carbon source exchanges need to be changed to be
%reversible too

[~,cs_ind]=ismember(media.carbonreac,media_model.rxns);

media_model.rev(unique(cs_ind))=1;
media_model.lb(unique(cs_ind))=-1000;
media_model_irrev = convertToIrreversible(media_model);
n_r=size(media_model_irrev.S,2);


if strcmp(biomass_choose,'wt')
    biomass=find(contains(media_model_irrev.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M'));
elseif strcmp(biomass_choose,'core')
    biomass=find(contains(media_model_irrev.rxns,'BIOMASS_Ec_iJO1366_core_53p95M'));
end

%Upper and lower bound from media and growth conditions
lb=zeros(n_r,31);
ub=zeros(n_r,31);

for i=1:31
    lb(:,i)=media_model_irrev.lb;
    ub(:,i)=media_model_irrev.ub;
    
    %minimal medium index 
    min_med_f=cs_ind(i); % the forward reaction==secretion
    min_med_b=media_model_irrev.match(min_med_f); % the backward reaction==uptake
    ub(min_med_f,i)=0; % block the production of secretion of the minimal carbon source
    ub(min_med_b,i)=15;
    
    limited_f=setdiff(cs_ind,min_med_f); 
    limited_b=media_model_irrev.match(limited_f); % reactions which need to be blocked
    ub(limited_b,i)=0; % block the uptake of other carbon sources
    
    lb(biomass,i)=0.95*media.growth(i); %WT biomass index=14, core biomass=19
    ub(biomass,i)=1.05*media.growth(i); %media(:,2)= growth rate

end
end