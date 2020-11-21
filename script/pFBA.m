function[V_solution_p]=pFBA(model_irrev,lb,ub,scale,biomass_choose)
n_cond=size(lb,2);
[n_m,n_r]=size(model_irrev.S);

f=ones(n_r,1);
V_solution_p=zeros(n_r,n_cond);
if strcmp(biomass_choose,'wt')
    biomass=find(contains(model_irrev.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M'));
elseif strcmp(biomass_choose,'core')
    biomass=find(contains(model_irrev.rxns,'BIOMASS_Ec_iJO1366_core_53p95M'));
end

for cond=1:n_cond
    m=struct();
    m.obj = f;
    Aeq=model_irrev.S;
    beq=zeros(n_m,1);
    Aeq(:,biomass)=Aeq(:,biomass)*scale;
    
    m.A = sparse(Aeq);
    n = size(m.A, 2);
    m.vtype = repmat('C', n, 1);
    m.sense =  repmat('=',size(Aeq,1),1);
    m.rhs = full( beq(:));
    
    m.lb=lb(:,cond)*scale; 
    m.ub=ub(:,cond)*scale;
    m.lb(biomass)=lb(biomass,cond);
    m.ub(biomass)=ub(biomass,cond);
    
    params = struct();
    params.FeasibilityTol=1e-6;
    x = gurobi(m,params);
    
    if ~strcmp(x.status,'INFEASIBLE')
        sol_p=x.x;
        V=sol_p(1:n_r);
        V_solution_p(:,cond)=V;
    end
end
