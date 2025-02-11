function [x,fval]=Join_Flexibility_New(Fn,pz,SF_1,SF_2,U_1,U_2,Gn,Miter, N_Coord, elements,S,T,St,i_v,Sd)
    N=2*length(pz)+6;
    ub=[ones(N-6,1);1.15*ones(6,1)];
    lb=[0.5*ones(N-6,1);0.85*ones(6,1)];
    warning('off', 'all');
    if Fn=="Genetic Algorithm"
        hybridopts = optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'Display', 'off');
        options = optimoptions('ga','Display','none','FunctionTolerance',1e-6,'PopulationSize',Gn, 'MaxGenerations',Miter, 'UseParallel', true, 'UseVectorized', false,'HybridFcn',{@patternsearch,hybridopts});
        fun= @(x)DSM_KFlexible_New(x, SF_1, SF_2, U_1, U_2, S,T, pz, N_Coord,elements,St,i_v,Sd);
        [x,fval] = ga(fun,N,[],[],[],[],lb,ub,[],[],options); % Run the genetic algorithm
    elseif Fn=="Particle Swarm"
        hybridopts = optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'Display', 'off');
        options = optimoptions('particleswarm','Display','none','FunctionTolerance',1e-6,'SwarmSize',Gn, 'MaxIterations',Miter, 'UseParallel', true, 'UseVectorized', false,'HybridFcn',{@patternsearch,hybridopts});
        fun= @(x)DSM_KFlexible_New(x, SF_1, SF_2, U_1, U_2, S,T, pz, N_Coord,elements,St,i_v,Sd);
        [x,fval] = particleswarm(fun,N,lb,ub,options); % Run the particle swarm algorithm
    elseif Fn=="Simulated Annealing"
        hybridopts = optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'Display', 'off');
        options = optimoptions('simulannealbnd','Display','none','FunctionTolerance',1e-6,'ReannealInterval',Gn,'MaxFunctionEvaluations',ceil(Miter*Gn/2),'HybridFcn',{@patternsearch,hybridopts});
        fun= @(x)DSM_KFlexible_New(x, SF_1, SF_2, U_1, U_2, S,T, pz, N_Coord,elements,St,i_v,Sd);
        [x,fval] = simulannealbnd(fun,ub,lb,ub,options); % Run the simulated annealing algorithm
    end
end