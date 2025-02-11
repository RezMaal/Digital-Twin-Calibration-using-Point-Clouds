E=210*10^9; %Young's Modulus
G=81*10^9; %Shear Modulus

load("Element_Properties.mat");
load("FEM_Displacements_New.mat");

St=1; %% support nodes
d1=Dis_Sap_E1;
d2=Dis_Sap_E2;
D_1= reshape(d1',[],1);
D_2= reshape(d2',[],1);
St=sort(St,'descend');
K = assembleGlobalStiffnessMatrix(ones(size(Element_Start,1),2),Node_Coor, [Element_Start,Element_End], A_Iy_Iz_J, E, G);

% Adjusting Stiffness Matrix
for i=1:length(St)
  K(6*St(i)-5:6*St(i),:)=[];       
  K(:,6*St(i)-5:6*St(i))=[];     
  D_1(6*St(i)-5:6*St(i))=[];
  D_2(6*St(i)-5:6*St(i))=[];
end

F_1=K*D_1;
F_2=K*D_2;
[S,T] = assembleGlobalStiffnessMatrix_New(Node_Coor, [Element_Start,Element_End], A_Iy_Iz_J, E, G);
Ax=Node_Coor(Element_End,:)-Node_Coor(Element_Start,:); Ax=Ax./vecnorm(Ax,2,2);
pz=find(Ax(:,3)==1);

load("Real_Deformed.mat");
warning ('off','all');
aa=unique(Node_Coor(:,3));
i_v=[find(Node_Coor(:,3)==aa(1)),find(Node_Coor(:,3)==aa(2)),find(Node_Coor(:,3)==aa(3))];
i_r=setdiff(1:size(d1,1),St);
Zn=[Node_Deformed_1(i_r,3), Node_Deformed_2(i_r,3)];

sim=50; % number of simulations
xp=zeros(sim,2*length(pz));
Par=zeros(sim,2);
Gen=zeros(sim,2);
SA=zeros(sim,2);

% Pre-allocation
Gn=100;
az=Zn(:,2)-Zn(:,1);
Par_Accuracy=zeros(sim,5);
Par_Disp=zeros(sim,length(Nodes));
Par_Z_Score=zeros(sim,length(i_r));
Gen_Accuracy=zeros(sim,5);
Gen_Disp=zeros(sim,length(Nodes));
Gen_Z_Score=zeros(sim,length(i_r));
SA_Accuracy=zeros(sim,5);
SA_Disp=zeros(sim,length(Nodes));
SA_Z_Score=zeros(sim,length(i_r));
ss=3; % coefficient of standard deviation

for l=1:sim
    fprintf('Iteration number: %d\n', l);

    tic;
        [x,fval]=Join_Flexibility_New("Particle Swarm", pz, F_1, F_2, Zn(:,1), Zn(:,2), Gn, 1000, Node_Coor, [Element_Start,Element_End],S,T,St,i_v,ss*Node_RMSE(i_r));
    t=toc;

    Par(l,:)=[fval,t];
    [~,U1,U2]= DSM_KFlexible_New(x, F_1, F_2, Zn(:,1), Zn(:,2), S,T, pz, Node_Coor, [Element_Start,Element_End],St,i_v,ss*Node_RMSE(i_r));
    Acc=abs(az-U2(:,3)+U1(:,3));
    Par_Disp(l,i_r)=Acc;
    Par_Accuracy(l,:)=[prctile(Acc,25),median(Acc),prctile(Acc,95),mean(Acc),sum((Acc./Node_RMSE(i_r))<=3)/length(i_r)];
    Par_Z_Score(l,:)=(Acc./Node_RMSE(i_r))';
    fprintf('Particle Swarm Finished: %d\n',t);

    tic;
        [x,fval]=Join_Flexibility_New("Genetic Algorithm", pz, F_1, F_2, Zn(:,1), Zn(:,2), Gn, 1000, Node_Coor, [Element_Start,Element_End],S,T,St,i_v,ss*Node_RMSE(i_r));
    t=toc;

    Gen(l,:)=[fval,t];  
    [~,U1,U2]= DSM_KFlexible_New(x, F_1, F_2, Zn(:,1), Zn(:,2), S,T, pz, Node_Coor, [Element_Start,Element_End],St,i_v,ss*Node_RMSE(i_r));
    Acc=abs(az-U2(:,3)+U1(:,3));
    Gen_Disp(l,i_r)=Acc;
    Gen_Accuracy(l,:)=[prctile(Acc,25),median(Acc),prctile(Acc,95),mean(Acc),sum((Acc./Node_RMSE(i_r))<=3)/length(i_r)];
    Gen_Z_Score(l,:)=(Acc./Node_RMSE(i_r))';
    fprintf('Genetic Algorithm Finished: %d\n',t);
    
    tic;
        [x,fval]=Join_Flexibility_New("Simulated Annealing", pz, F_1, F_2, Zn(:,1), Zn(:,2), Gn, 1000, Node_Coor, [Element_Start,Element_End],S,T,St,i_v,ss*Node_RMSE(i_r));
    t=toc;

    SA(l,:)=[fval,t];
    [~,U1,U2]= DSM_KFlexible_New(x, F_1, F_2, Zn(:,1), Zn(:,2), S,T, pz, Node_Coor, [Element_Start,Element_End],St,i_v,3*Node_RMSE(i_r));
    Acc=abs(az-U2(:,3)+U1(:,3));
    SA_Disp(l,i_r)=Acc;
    SA_Accuracy(l,:)=[prctile(Acc,25),median(Acc),prctile(Acc,95),mean(Acc),sum((Acc./Node_RMSE(i_r))<=3)/length(i_r)];
    SA_Z_Score(l,:)=(Acc./Node_RMSE(i_r))';
    fprintf('Simulated Annealing Finished: %d\n',t);
    
end

clearvars -except Par_Accuracy Par_Disp Par_Z_Score Gen_Accuracy Gen_Disp Gen_Z_Score SA_Accuracy SA_Disp SA_Z_Score