function [milp]=getkmax_av(Kapp_matrix,Kapp_matrix2,homo)
k_max=max(Kapp_matrix,[],2);
k_max2=max(Kapp_matrix2,[],2);
k_max(k_max==0)=NaN;
k_max2(k_max2==0)=NaN;
kmax3=nanmean([k_max,k_max2],2);

ind=find(~isnan(kmax3));
milp.kmax=kmax3(ind);
milp.reac=string(homo.reac(ind));
milp.reacind=homo.reacind(ind);
milp.genes=string(homo.genes(ind));
milp.kapp=Kapp_matrix(ind,:);

milp.kapp=Kapp_matrix(ind,:);
milp.abun=homo.abun(ind,:);

milp.kmax1=k_max;
milp.kmax2=k_max2;
end