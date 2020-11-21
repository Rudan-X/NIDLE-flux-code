function [idle]=idle_enzyme(model,g_vect,flux,abundance_n,thre)

idle.count=zeros(size(abundance_n.abun,2),1);
idle.prop=zeros(size(abundance_n.abun,2),1);
idle.unique=zeros(size(abundance_n.abun,2),1);
idle.abun=zeros(size(abundance_n.abun,2),1);
idle.abun_count=zeros(size(abundance_n.abun,2),1);
idle.abun_total=zeros(size(abundance_n.abun,2),1);
for cond=1:size(abundance_n.abun,2)
    notnan=find(~isnan(abundance_n.abun(:,cond)));
    [gene_u,ind_u,~]=unique(abundance_n.genes(notnan));
    
    count=0;
    idle_abun=0;
    active_gpr=find(g_vect(:,cond)==1);
    for gene=1:length(ind_u)
        reac_ind=intersect(find(contains(model.grRules,gene_u(gene))),active_gpr);
        cata_reac=flux(reac_ind,cond);
        ind=find(cata_reac>thre);
        if isempty(ind)
            count=count+1;
            idle_abun=idle_abun+abundance_n.abun(notnan(ind_u(gene)),cond);
        end
    end
    
    idle.count(cond)=count;
    idle.unique(cond)=length(gene_u);
    idle.prop(cond)=count/length(gene_u);
    idle.abun_count(cond)=idle_abun;
    idle.abun_total(cond)=sum(abundance_n.abun(notnan(ind_u),cond));
    idle.abun(cond)=idle_abun/sum(abundance_n.abun(notnan(ind_u),cond));
    idle=orderfields(idle);
end
end
