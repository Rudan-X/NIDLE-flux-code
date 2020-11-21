function [new]=NIDLE(model,g_vect,lb,ub,ep,scale,biomass_choose)

if strcmp(biomass_choose,'wt')
    biomass=find(contains(model.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M'));
elseif strcmp(biomass_choose,'core')
    biomass=find(contains(model.rxns,'BIOMASS_Ec_iJO1366_core_53p95M'));
end


[n_m,n_r]=size(model.S);
n_cond=size(lb,2);
new.y=zeros(n_r,n_cond);
new.flux=zeros(n_r,n_cond);


for cond=1:size(lb,2)
    ind_active=find(g_vect(:,cond)~=0);
    %Assign y variables to reactions with not nan abundance
    ind_f=find(contains(model.rxns(ind_active),'_f'));
    ind_b=find(contains(model.rxns(ind_active),'_b'));

    
    n_y=length(ind_active); 
    y0=n_r+1; %index where the y variable starts


    %%Maximizing sum of y

    %a) Inequality matrix
    %Only Y variable constraints only for reactions associated with genes of known abundance
    %Two inequalities of length(abun_genes) as rows, and n_r+n_y as columns
    %One inequality of length of reversible abundance reactions as rows
    I=eye(n_r);
    I_abun=I(ind_active,:);
    I_y=eye(n_y);
    

    Vmin=lb(:,cond)*scale;
    Vmax=ub(:,cond)*scale;
    Vmin(biomass)=lb(biomass,cond);
    Vmax(biomass)=ub(biomass,cond);

    A1=[ -I_abun, -(Vmin(ind_active)-ep).*I_y;
        I_abun, -(Vmax(ind_active)-ep).*I_y];

    A2=zeros(length(ind_f),n_y);
    for i=1:length(ind_f)
        A2(i,ind_f(i))=1;
        A2(i,ind_b(i))=1;
    end

    A=[A1;[zeros(length(ind_f),n_r),A2]];

    b=[-Vmin(ind_active);ep*ones(n_y,1);ones(length(ind_f),1)];

    Aeq=[model.S,zeros(n_m,n_y)];
    Aeq(:,biomass)=Aeq(:,biomass)*scale;
    beq=zeros(n_m,1);
   
    f=[zeros(1,n_r),-ones(1,n_y)];

    
    int=n_r+1:n_r+n_y;

    lower=[Vmin;zeros(n_y,1)];
    upper=[Vmax;ones(n_y,1)];

    
    m=struct();
    m.obj = f;
    m.A = [sparse(A); sparse(Aeq)]; % A must be sparse
    n = size(m.A, 2);
    m.vtype = repmat('C', n, 1);
    m.vtype(int) = 'I';
    m.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];
    m.rhs = full([b(:); beq(:)]); % rhs must be dense
    m.lb=lower;
    m.ub=upper;
    
    params = struct();
    params.FeasibilityTol=1e-6;
    params.IntFeasTol=1e-3;
    x = gurobi(m,params);

    if ~strcmp(x.status,'INFEASIBLE')
        sol=x.x;
        y=sol(y0:end);

        %Minimizing sum of fluxes
        Aeq2=[Aeq;zeros(1,n_r), ones(1,n_y)];
        beq2=[beq;sum(y)-1];


        f2=zeros(1,n_r+n_y);
        f2(1:n_r)=1; %minimize the total sum of fluxes
        %%sol_new=intlinprog_adap(f2,int,A,b,Aeq2,beq2,lower,upper);
        m.obj = f2;
        m.A = [sparse(A); sparse(Aeq2)]; % A must be sparse

        m.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq2,1),1)];
        m.rhs = full([b(:); beq2(:)]); % rhs must be dense

        x2 = gurobi(m, params);
        if ~strcmp(x2.status,'INFEASIBLE')
            sol_new=x2.x;
            y_new=sol_new(n_r+1:end);

            new.flux(:,cond)=sol_new(1:n_r);
            new.y(ind_active,cond)=y_new;
        end
    end
end
