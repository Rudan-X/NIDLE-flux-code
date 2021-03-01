function [rules_type,vect]=parse_rules(model,abundance_u)
n_r=size(model.S,2);

%1) Find the type of GPR rules_type
ind1=find(contains(model.grRules,' or ('));
ind2=find(contains(model.grRules,') or '));
complex_iso=union(ind1,ind2);

ind1=find(contains(model.grRules,' and ('));
ind2=find(contains(model.grRules,') and '));
complex_and=union(ind1,ind2);

only_and=setdiff(complex_and,complex_iso);

rules_type=zeros(n_r,1);
%rules_type==0: not defined yet
%rules_type==1: no GPR rules_type
%rules_type==2: single gene
%rules_type==3: isoenzyme
%rules_type==4: complex
%rules_type==5: complex isoenzyme
%rules_type==6: complex 'complex'
rules_type(complex_iso)=5;
rules_type(only_and)=6;
l=1:n_r;
rest=setdiff(l,union(complex_iso,only_and));
for i=1:n_r
    if find(rest==i)
        if isempty(model.grRules{i})
            rules_type(i)=1;
        elseif contains(model.grRules{i}, 'or')
            rules_type(i)=3;
        elseif contains(model.grRules{i}, 'and')
            rules_type(i)=4;
        else
            rules_type(i)=2;
        end
    end
end



%Part2: assign vector g according to the conditions

n_cond=size(abundance_u.abun,2);

g_vect=NaN(n_r,n_cond);

single=find(rules_type==2);
iso=find(rules_type==3);
compl=find(rules_type==4);
compl2=find(rules_type==6);
compl_iso=find(rules_type==5);

for cond=1:n_cond
    %1. single genes
    
    for i=1:length(single)
        reac_ind=single(i);
        [~,gene_ind]=ismember(model.grRules(reac_ind),abundance_u.genes);
        if gene_ind~=0     
           g_vect(reac_ind,cond)=abundance_u.abun(gene_ind,cond);
        end
    end 
    
    %2. isoenzymes
    for i=1:length(iso)
        reac_ind=iso(i);
        genes=split(model.grRules{reac_ind},' or ');
        [~,gene_ind]=ismember(genes,abundance_u.genes);
        gene_ind2=gene_ind(gene_ind~=0);
        if ~isempty(gene_ind2)
            g_vect(reac_ind,cond)=nanmax(abundance_u.abun(gene_ind2,cond));
            %we only care if there is at least one enzyme with measured
            %abundance nanmax(NaN,NaN)=NaN
        end
    end

    %2. complex
    for i=1:length(compl)
        reac_ind=compl(i);
        genes=split(model.grRules{reac_ind},' and ');
        [~,gene_ind]=ismember(genes,abundance_u.genes);
        gene_ind2=gene_ind(gene_ind~=0);
        if length(gene_ind2)==length(gene_ind)
            ind=find(isnan(abundance_u.abun(gene_ind2,cond)));
            if ~isempty(ind)
                g_vect(reac_ind,cond)=NaN;
            else
                g_vect(reac_ind,cond)=min(abundance_u.abun(gene_ind2,cond));
            end
        end
    end
    
   %3. complex complex
    for i=1:length(compl2)
        reac_ind=compl2(i);
        genes=strsplit(model.grRules{reac_ind},' and ');
        genes=erase(genes,'(');
        genes=erase(genes,')');
        [~,gene_ind]=ismember(genes,abundance_u.genes);
        gene_ind2=gene_ind(gene_ind~=0);
        if length(gene_ind2)==length(gene_ind)
            ind=find(isnan(abundance_u.abun(gene_ind2,cond)));
            if ~isempty(ind)
                g_vect(reac_ind,cond)=NaN;
            else
                g_vect(reac_ind,cond)=min(abundance_u.abun(gene_ind2,cond));
            end
        end
    end 
    
    % 4. complex isoenzymes
    for i=1:length(compl_iso)
        reac_ind=compl_iso(i);
        complex=split(model.grRules{reac_ind},' or ');
        count2=NaN(length(complex),1);
        for m=1:length(complex)
            genes=split(complex(m),' and ');
            genes=erase(genes,'(');
            genes=erase(genes,')');
            [~,gene_ind]=ismember(genes,abundance_u.genes);
            gene_ind2=gene_ind(gene_ind~=0);
            if length(gene_ind2)==length(gene_ind)
                ind=find(isnan(abundance_u.abun(gene_ind2,cond)));
                if ~isempty(ind)
                    count2(m)=NaN;
                else
                    count2(m)=min(abundance_u.abun(gene_ind2,cond));
                end
            end
        end
        g_vect(reac_ind,cond)=nanmax(count2);
    end
end

vect=zeros(n_r,n_cond);
vect(~isnan(g_vect))=1;

end