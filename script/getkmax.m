function [milp,k_max]=getkmax(Kapp_matrix,V_matrix,homo)
[k_max,ind_max]=max(Kapp_matrix,[],2);
ind=find(k_max~=0);
milp.kmax=k_max(ind);
milp.reac=string(homo.reac(ind));
milp.reacind=homo.reacind(ind);
milp.genes=string(homo.genes(ind));
milp.kapp=Kapp_matrix(ind,:);
milp.conditions=ind_max(ind);
milp.kapp=Kapp_matrix(ind,:);
milp.abun=homo.abun(ind,:);
milp.v=V_matrix(ind,:);
k_max(k_max==0)=NaN;
milp.kmax_comp=k_max;
end