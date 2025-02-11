function [fun,U1,U2]= DSM_KFlexible_New(x, SF_1, SF_2, U_1, U_2, S,T, pz, N_Coord,elements,St,i_v,Sd)
    n=size(pz,1);    
    w=ones(length(elements),2);
    w(pz,1)=x(1:n);
    w(pz,2)=x(n+1:2*n);
    warning('off', 'all');
    K=Correction_Rotation_Assembly(S,T,w,N_Coord,elements);
    St=sort(St,'descend');  % will need to remove rows/columns from back going forward

    aa=setdiff(1:length(SF_1)/6+length(St),St)';
    S1=zeros(length(SF_1)/6+length(St),6);
    S1(aa,:)=reshape(SF_1,6,[])';
    S2=zeros(length(SF_2)/6+length(St),6);
    S2(aa,:)=reshape(SF_2,6,[])';
    for i=1:size(i_v,2)
        S1(i_v(:,i))=S1(i_v(:,i)).*x(2*n+i);
        S2(i_v(:,i))=S2(i_v(:,i)).*x(2*n+i+3);
    end
    S_1=reshape(S1',[],1);
    S_2=reshape(S2',[],1);

    for i=1:length(St)
      K(6*St(i)-5:6*St(i),:)=[];       % remove row
      K(:,6*St(i)-5:6*St(i))=[];       % remove column
      S_1(6*St(i)-5:6*St(i),:)=[];
      S_2(6*St(i)-5:6*St(i),:)=[];
    end
    
    U_1_e=K\S_1;
    U_2_e=K\S_2;
    U1=reshape(U_1_e,6,[])';
    U1=U1(:,1:3);
    U2=reshape(U_2_e,6,[])';
    U2=U2(:,1:3);
    D1_2=abs((abs(U2(:,3)-U1(:,3))-abs(U_2-U_1)));
    % pcd=ones(length(U2),1);
    % pcd(D1_2-Sd>0)=0.001;
    if ~isempty(D1_2-Sd>0)
        el=0.001/sum(D1_2-Sd>0);
    else
        el=1;
    end
    % fun=mean(D1_2)/prod(pcd);
    fun=mean(D1_2)/el;
end
function C=Correction_Rotation_Assembly(S,T,w,N_Coord,elements)
    numNodes=size(N_Coord,1);
    % Number of elements
    numElements = size(elements, 1);
    startCoord = N_Coord(elements(:,1), :);
    endCoord = N_Coord(elements(:,2), :);
    L = vecnorm(endCoord - startCoord,2,2);
    % Initialize global stiffness matrix
    C = zeros(6*numNodes, 6*numNodes);
    wb=4-w(:,1).*w(:,2);    
    ww1=(4*w(:,2)-2*w(:,1)+w(:,1).*w(:,2))./wb;
    ww2=-2.*L.*w(:,1).*(1-w(:,2))./wb;
    ww3=6./L.*(w(:,1)-w(:,2))./wb;
    ww4=3*w(:,1).*(2-w(:,2))./wb;
    ww5=(4*w(:,1)-2*w(:,2)+w(:,1).*w(:,2))./wb;
    ww6=2.*L.*w(:,2).*(1-w(:,1))./wb;
    ww7=3*w(:,2).*(2-w(:,1))./wb;
    % Loop over each element to compute local and global stiffness matrices
    for i = 1:numElements
        
        R1=eye(6); 
        R1(2,2)=ww1(i); R1(2,6)=ww2(i); R1(6,2)=ww3(i); R1(6,6)=ww4(i);
        R2=eye(6);
        R2(2,2)=ww5(i); R2(2,6)=ww6(i); R2(6,2)=ww3(i); R2(6,6)=ww7(i);
   
        % Transform the local stiffness matrix to global coordinates
        kGlobal = T{i,1}'* ([[S{i,1}(1:6,1:6)* R1,S{i,1}(1:6,7:12)*R2];[S{i,1}(7:12,1:6)* R1,S{i,1}(7:12,7:12)*R2]]) * T{i,1};
        
        % Global DOF indices for start and end nodes
        dof = [6*elements(i,1)-5:6*elements(i,1), 6*elements(i,2)-5:6*elements(i,2)];
        
        % Assemble the global stiffness matrix
        C(dof, dof) = C(dof, dof) + kGlobal;
        
    end
end