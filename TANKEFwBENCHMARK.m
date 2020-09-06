%% Initialisation %%
clc
close all
clearvars
format long g
warning off
%%
cprintf('\n')
cprintf('[10/255,30/255,150/255]','| |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| |\n')
cprintf('[10/255,30/255,150/255]','| |                                                  | |\n')
cprintf('[10/255,30/255,150/255]','| |                                                  | |\n')
cprintf('[10/255,30/255,150/255]','| |             TANKEF vBENCHMARK VERSION            | |\n')
cprintf('[210/255,0,0]          ','| |               (DO NOT CHANGE VALUES)             | |\n')
cprintf('[10/255,30/255,150/255]','| |  MATLAB Codes for FEA of Pressurized Water Tank  | |\n')
cprintf('[10/255,30/255,150/255]','| |                                                  | |\n')
cprintf('[10/255,30/255,150/255]','| |       Mariano Arnaiz and Damian Caballero        | |\n')
cprintf('[10/255,30/255,150/255]','| |           Central University of Venezuela        | |\n')
cprintf('[10/255,30/255,150/255]','| |           University of Hertfordshire            | |\n')
cprintf('[10/255,30/255,150/255]','| |                      (2020)                      | |\n')
cprintf('[10/255,30/255,150/255]','| |--------------------------------------------------| |\n')
cprintf('[10/255,30/255,150/255]','| |--------------------------------------------------| |\n')
cprintf('\n')

%% INPUT DATA:
%% Constants
g_a=9.81; % Gravity acceleration [m/s^2]

%% Mechanical Parameters of the Tank and fluid.
E=30e9;              % Youngs Modulus [Pa]
v=1/6;               % Poissons ratio
Density_fluid=1000; % Density of the tank's fluid [kg/m3]
H_fluid=3; % Height of fluid in the tank [m]

%% Tank Geometry and Mesh Generation
R=2.5;  % Tanks radius [m]
H=3.1;    % Hight of the tank from 0 [m]
wall=0.1; % Tanks wall thickness [m]

% Tank's Perimeter. Make sure that it is closed
pv =[0  0;  R*2  0; R*2  H;   R*2-wall  H;   R*2-wall  wall; wall  wall; wall  H;   0  H;   0  0];

% Generation of mesh
cprintf('black','Generating Mesh with DISTMESH. Please Wait.\n')
fd = { 'l_dpolygon', [], pv }; %Creating polygon function
fh = @(p) ones(size(p,1),1); % Node's Density Function.
[p,TC] = distmesh( fd, fh, 0.02, [0,0; 2*R,H], pv ); % Create Mesh.

%% Rearrage Nodes information
%% Save values to vectors
p0=p; %save original p for plot
p=round(p,4); %make sure zeros are real zeros
vXnodes=p(:,1); %X Coordinates of Nodes
vYnodes=p(:,2); %Z Coordinates of Nodes

% Save the nodes of the boundary nodes
nodes_out_bottom=find((vYnodes==0)==1); % Nodes on the base of the tank
nodes_in_bottom=find((vYnodes==wall)==1); % Nodes on the bottom of the tank
nodes_in_left=find((vXnodes==wall)==1); % Nodes on the inside left wall
nodes_out_left=find((vXnodes==0)==1); % Nodes on the outside left wall
nodes_in_right=find((vXnodes==R*2-wall)==1); % Nodes on the inside right wall
nodes_out_right=find((vXnodes==R*2)==1);  % Nodes on the outside right wall
nodes_top=(find((vYnodes==H)==1));  % Nodes at the top of the tank (end of the walls);


% Number of elements
NUMBER_OF_ELEMENTS=size(TC,1); % Number of elements in the mesh
vNodesNames=1:size(p(:,1),1); % Vector of nodes names

%% FEM MATRIXES
%% Compute D: material property matrix (Rigidiy Matrix)
cprintf('[0/255, 142/255, 158/255]','Computing FE Matrixes.\n')
C1= E*(((1-v))/((1+v)*(1-2*v))); % Constants outside the matriz
C2=v/(1-v); % Constant in the matrix
C3=(1-2*v)/(2*(1-v)); % Constant in the matrix

M1=[1   C2   0   C2 ;  %
    C2  1    0   C2 ;  % D MATRIX FOR ALL ELEMENTS OF THE TANK
    0   0    C3  0  ;  % 
    C2  C2   0   1  ]; %

D= C1*M1; % Rigidity Matrix 

%% Compute det J and A for each element. This information is stored in the TC matrix
for a=1:size(TC,1)
    TC(a,4)=(vXnodes(TC(a,1))-vXnodes(TC(a,3)))*(vYnodes(TC(a,2))-vYnodes(TC(a,3)))-(vXnodes(TC(a,2))-vXnodes(TC(a,3)))*(vYnodes(TC(a,1))-vYnodes(TC(a,3))); %detJ for each element 
    TC(a,5)=abs(TC(a,4))*0.5; % Area of each element    
end

%% Compute radial coordinates and mean radial coordonate for each element
%Find Radial Coordinate for each node
RADIAL=zeros(size(vXnodes,1),1);
for uu=1:size(vXnodes,1)
if vXnodes(uu)<=R
    RADIAL(uu)=R-vXnodes(uu);
else
    RADIAL(uu)=vXnodes(uu)-R;
end    
end
%eliminate zero at the center
indexmin=find(RADIAL==0);
RADIAL(indexmin)=RADIAL(indexmin)+0.0001;

%Find Mean Radial Coordinate for each Element
for a=1:size(TC,1)
mean_x=(vXnodes(TC(a,1))+vXnodes(TC(a,2))+vXnodes(TC(a,3)))/3; % mean x coordinate of each element.

if mean_x<=R
    TC(a,6)=R-mean_x;
else
    TC(a,6)=mean_x-R;
end
end
%eliminate zero at the center
indexmin=find(TC(:,6)==0);
TC(indexmin,6)=TC(indexmin,6)+0.0001;

%% Compute B matrix for each element
B=zeros(4,6,size(TC,1)); %initizialize B as zeros for speed
for f=1:size(TC,1)
Nr=(1/3)/TC(f,6); % N1 = N2 = N3 = 1/3 at the centroid of the element
B(:,:,f)=[(vYnodes(TC(f,2))-vYnodes(TC(f,3)))/TC(f,4) 0 (vYnodes(TC(f,3))-vYnodes(TC(f,1)))/TC(f,4) 0 (vYnodes(TC(f,1))-vYnodes(TC(f,2)))/TC(f,4) 0 ;
          0 (vXnodes(TC(f,3))-vXnodes(TC(f,2)))/TC(f,4) 0 (vXnodes(TC(f,1))-vXnodes(TC(f,3)))/TC(f,4) 0 (vXnodes(TC(f,2))-vXnodes(TC(f,1)))/TC(f,4) ;
          (vXnodes(TC(f,3))-vXnodes(TC(f,2)))/TC(f,4) (vYnodes(TC(f,2))-vYnodes(TC(f,3)))/TC(f,4) (vXnodes(TC(f,1))-vXnodes(TC(f,3)))/TC(f,4) (vYnodes(TC(f,3))-vYnodes(TC(f,1)))/TC(f,4) (vXnodes(TC(f,2))-vXnodes(TC(f,1)))/TC(f,4) (vYnodes(TC(f,1))-vYnodes(TC(f,2)))/TC(f,4);
          Nr 0 Nr 0 Nr 0];
end

%% Compute Ke (element stiffness) for each element
Ke=zeros(6,6,size(TC,1)); %initizialize Ke as zeros
for g=1:size(TC,1)
    %Ke=2*pi*r*Ae*B'*D*B
    Ke(:,:,g)=2*pi*TC(g,6)*TC(g,5)*B(:,:,g)'*D*B(:,:,g);
end

%% Compute K: global stiffenss matric
cprintf('*black','Assembling Global Stiffness Matrix.\n')
K=zeros(max(vNodesNames(:))*2);
for s=1:size(TC,1)
    GDL=[2*TC(s,1)-1 2*TC(s,1) 2*TC(s,2)-1 2*TC(s,2) 2*TC(s,3)-1 2*TC(s,3)]; %Find DOF of all 6 nodes
    for m=1:6
        for n=1:6
            Value=Ke(m,n,s);
            K(GDL(m),GDL(n))=Value+K(GDL(m),GDL(n)); % Add Ke values to K according to DOF
        end
    end
end
K_no_restrictions=K; % Save the original K with no restructions
  
%% Set restriction at the base of the tank by penalty approach
cprintf('black','Setting Restrictions. Bottom of the tank fixed to the ground.\n')

%Find the numbers of the degrees of freedom to be restricted
dof_R_X=(nodes_out_bottom*2)-1; % Degrees of freedom restricted (in the X direction)
dof_R_Y=(nodes_out_bottom*2);% Degrees of freedom restricted (in the Y direction)
dof_R=sort([dof_R_X;dof_R_Y]); % All Degrees of freedom restricted

PENALTY=max(K(:))*10000; % Compute the PENALTY CONSTANT value

%Add a large value to K's diagonal at the coordinates of the restricted degrees of freedom 
for o=1:size(dof_R,1)
        K(dof_R(o),dof_R(o))=K(dof_R(o),dof_R(o))+PENALTY; % Add PENALTY CONSTANT values to K according to DOF
end

% Add restrictions at the top of the tank. If no restriction here comment
% lines from 164 to 171
dof_R_X_top=(nodes_top*2)-1; % Degrees of freedom restricted (in the X direction)
dof_R_Y_top=(nodes_top*2);% Degrees of freedom restricted (in the Y direction)
dof_R_top=sort([dof_R_X_top;dof_R_Y_top]);

for o=1:size(dof_R_top,1)
        K(dof_R_top(o),dof_R_top(o))=K(dof_R_top(o),dof_R_top(o))+PENALTY; % Add PENALTY CONSTANT values to K according to DOF
end

%% Forces Vectors
F=zeros(size(K,1),1); %Empty forces vector
cprintf('*black','Computing Forces Vector.\n')

%% Weight of the fluid (Com\ent lines 166 to 184 so this force is not considered)
% %Vertical force on the base of the tank
% F_y=g_a*Density_fluid*H_fluid;
% %Find the numbers of the degrees of freedom to be loaded
% dof_L_Y=(nodes_in_bottom*2);% Degrees of freedom to be loaded (in the Y direction)
% 
% DBT=[vXnodes(nodes_in_bottom) nodes_in_bottom RADIAL(nodes_in_bottom)];
% sDBT=sortrows(DBT,1); % Organize Nodes in the wall
% 
% %For each pair of nodes that limit the top of the elemnt
% for ii=1:size(sDBT,1)-1
%    Traction1=2*pi*(sDBT(ii+1,3)-sDBT(ii,3))*((2*sDBT(ii,3)+sDBT(ii+1,3))/6)*F_y;
%    Traction2=2*pi*(sDBT(ii+1,3)-sDBT(ii,3))*((sDBT(ii,3)+2*sDBT(ii+1,3))/6)*F_y; 
%     
% F(dof_L_Y(ii))=F(dof_L_Y(ii))-Traction1; %Negative for downwawrd direction
% F(dof_L_Y(ii+1))=F(dof_L_Y(ii+1))-Traction2; %Negative for downwawrd direction
% 
% end
% 
% F=-abs(F); %Make sure force is in negative

%% Horizontal forces on the left wall of the tank
dF_x_L=g_a*Density_fluid*((H_fluid+wall)-p(nodes_in_left,2));
indices = find(p(nodes_in_left,2)>=(H_fluid+wall)); %nodes over the fluid level
dF_x_L(indices) = 0; %remove loads from nodes over the fluid
dof_L_Y_left=(nodes_in_left*2)-1;% Degrees of freedom to be loaded (in the Y direction)

sDBT=sortrows([dF_x_L p(nodes_in_left,2)],2);
Fonwall=pi*((sDBT(1:end-1,1)+sDBT(2:end,1))*0.5).*((-sDBT(1:end-1,2)+sDBT(2:end,2)));
Fonwall(end)=0; %Keep the force in the final node = to 0 to respect equations
for u=1:size(dF_x_L,1)-1
F(dof_L_Y_left(u))=(F(dof_L_Y_left(u))-Fonwall(u));
F(dof_L_Y_left(u+1))=(F(dof_L_Y_left(u+1))-Fonwall(u));%Negative because of direction
end

%% Horizontal forces on the right wall of the tank
dF_x_R=g_a*Density_fluid*((H_fluid+wall)-p(nodes_in_right,2));
indices = find(p(nodes_in_right,2)>=(H_fluid+wall)); %nodes over the fluid level
dF_x_R(indices) = 0; %remove loads from nodes over the fluid
dof_L_Y_right=(nodes_in_right*2)-1;% Degrees of freedom to be loaded (in the Y direction)

sDBT=sortrows([dF_x_R p(nodes_in_right,2)],2); % Organize Nodes in the wall
Fonwall=pi*((sDBT(1:end-1,1)+sDBT(2:end,1))*0.5).*((-sDBT(1:end-1,2)+sDBT(2:end,2)));
Fonwall(end)=0; %Keep the force in the final node = to 0 to respect equations
for u=1:size(dF_x_R,1)-1
F(dof_L_Y_right(u))=(F(dof_L_Y_right(u))+Fonwall(u));
F(dof_L_Y_right(u+1))=(F(dof_L_Y_right(u+1))+Fonwall(u));%Negative because of direction
end

%% COMPUTE STRAIN: THIS REQUIRES INVERSION OF K
cprintf('*[210/255,0,0]','Computing Nodal Displacements. Please Wait.\n')
Q=K\F; % Compute strai
%% GET DISPLACEMENTS
%Store Displacements data in 2 columns
vQ=zeros(size(Q,1)/2,2);
t=1;
for r=1:2:size(Q,1)
  vQ(t,:)=[Q(r),Q(r+1)];
  t=t+1;
end
%Compute new coordinates
vXnodes_Q=vXnodes+vQ(:,1);
vYnodes_Q=vYnodes+vQ(:,2);

%% COMPUTE STRESSES
SIGMA=zeros(4,1,size(TC,1)); %initiate SIGMA
SIGMA_V=zeros(size(TC,1),9); %initiate SIGMA_V
%% Compute the total nodal displacement
node_displacement=[vQ(:,1) vQ(:,2)];   %Nodal displacement vector
%% Compute stresses (sigmas) for each element
for qq=1:size(TC,1)
    qe=[node_displacement(TC(qq,1)) node_displacement(TC(qq,1),2) node_displacement(TC(qq,2)) node_displacement(TC(qq,2),2) node_displacement(TC(qq,3)) node_displacement(TC(qq,3),2)]; %Element nodal displacement
    % SIGMA VECTOR: SigmaX SigmaY TaoXY SigmaO
    SIGMA(:,:,qq)=(D*B(:,:,qq)*qe')/1e6; % Stress calculation for each element
    %C and R from the Mohr's circle
    C=(SIGMA(1,1,qq)+SIGMA(2,1,qq))/2;
    RR=realsqrt((((SIGMA(1,1,qq)-SIGMA(2,1,qq))/2)^2)+(SIGMA(3,1,qq)^2));
    %Vector= [Xcentroid Zcentroid SigmaX SigmaY TaoXY SigmaO Sigma1 Sigma2 TaoMax]
    SIGMA_V(qq,:)=[(vXnodes_Q(TC(qq,1))+vXnodes_Q(TC(qq,2))+vXnodes_Q(TC(qq,3)))/3 (vYnodes_Q(TC(qq,1))+vYnodes_Q(TC(qq,2))+vYnodes_Q(TC(qq,3)))/3 SIGMA(:,:,qq)' C+RR C-RR RR ];
end 

%% FIGURES!
%% Figure 1: Mesh and Nodes
cprintf('black','Plotting Figures. Please Wait.\n')

figure('Name','Figure 1: Mesh and Nodes','NumberTitle','off')
%Plot Tank
patch( 'vertices', p0, 'faces', TC(:,1:3), 'facecolor', [180/255, 180/255, 180/255],'LineWidth',0.05)
hold on
%Plot Restricted nodes
plot(p(nodes_out_bottom,1),p(nodes_out_bottom,2),'s','MarkerEdgeColor','k','MarkerFaceColor','r')
%Plot Fluid
x = [wall 2*R-wall 2*R-wall wall];
y = [wall wall H_fluid+wall H_fluid+wall];
patch(x,y,'c')
axis equal
grid on
title('FEM Mesh')
xlabel('X [m]')
ylabel('Z [m]')
legend('Tank Walls', 'Restricted Nodes','Fluid')
axis tight

%% Figure 2: Displacements of Walls

figure('Name','Figure 2: Displacements of Walls','NumberTitle','off')
% Figure 2.1: Displacement on the Left Wall of the Tank
subplot(2,2,1)
plot(vXnodes(nodes_in_left),vYnodes(nodes_in_left),'.k')
hold on
plot(vXnodes_Q(nodes_in_left),vYnodes_Q(nodes_in_left),'r.')
grid on
title('Displacement on the Left wall of the Tank')
xlabel('X [m]')
ylabel('Z [m]')
legend('Initial Position','Final Position')

% Figure 2.2: Displacement on the Right Wall of the Tank
subplot(2,2,2)
plot(vXnodes(nodes_in_right),vYnodes(nodes_in_right),'.k')
hold on
plot(vXnodes_Q(nodes_in_right),vYnodes_Q(nodes_in_right),'r.')
grid on
title('Displacement on the Left wall of the Tank')
xlabel('X [m]')
ylabel('Z [m]')
legend('Initial Position','Final Position','Location','northwest')

% Figure 2.3: Displacement on the Bottom of the Tank
subplot(2,2,3:4)
plot(vXnodes(nodes_in_bottom),vYnodes(nodes_in_bottom),'.k')
hold on
plot(vXnodes_Q(nodes_in_bottom),vYnodes_Q(nodes_in_bottom),'r.')
grid on
title('Displacement on the Bottom of the Tank')
xlabel('X [m]')
ylabel('Z [m]')
legend('Initial Position','Final Position')

%% Figure 3: SigmaX,SigmaY,TaoXY, SigmaO plots

figure('Name','STRESSES in XY plane','NumberTitle','off')

v=p0; %Vertices
f=TC(:,1:3); % Face

%Figure 3.1: SigmaX
subplot(2,2,1)
col=SIGMA_V(:,3);
patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','LineStyle','None');
c=colorbar;
c.Label.String = 'Stress (MPa)';
colormap jet
title('Stress in X Direction (Sigma X in MPa)')
xlabel('X [m]')
ylabel('Z [m]')
grid on
axis tight

%Figure 3.2: SigmaY
subplot(2,2,2)
col=SIGMA_V(:,4);
patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','LineStyle','None');
c=colorbar;
c.Label.String = 'Stress (MPa)';
colormap jet
title('Stress in Y Direction (Sigma Y in MPa)')
xlabel('X [m]')
ylabel('Z [m]')
grid on
axis tight

%Figure 3.3: TaoXY
subplot(2,2,3)
col=SIGMA_V(:,5);
patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','LineStyle','None');
c=colorbar;
c.Label.String = 'Stress (MPa)';
colormap jet
title('Shear Stress in XY Plane (Tao XY in MPa)')
xlabel('X [m]')
ylabel('Z [m]')
grid on
axis tight

%Figure 3.4: SigmaO
subplot(2,2,4)
col=SIGMA_V(:,6);
patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','LineStyle','None');
c=colorbar;
c.Label.String = 'Stress (MPa)';
colormap jet
title('Stress in Radial Direction (Sigma O in MPa)')
xlabel('X [m]')
ylabel('Z [m]')
grid on
axis tight

%% Figure 4: Sigma1,Sigma2,TaoMAX plots

figure('Name','Principal Stresses','NumberTitle','off')

%Figure 4.1: Sigma1
subplot(2,2,1)
col=SIGMA_V(:,7);
patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','LineStyle','None');
c=colorbar;
c.Label.String = 'Stress (MPa)';
colormap jet
title('Principal Stress 1 (Sigma 1 in MPa)')
xlabel('X [m]')
ylabel('Z [m]')
grid on
axis tight

%Figure 4.2: Sigma2
subplot(2,2,2)
col=SIGMA_V(:,8);
patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','LineStyle','None');
c=colorbar;
c.Label.String = 'Stress (MPa)';
colormap jet
title('Principal Stress 2 (Sigma 2 in MPa)')
xlabel('X [m]')
ylabel('Z [m]')
grid on
axis tight

%Figure 4.3: TaoMAX
subplot(2,2,3)
col=SIGMA_V(:,9);
patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','LineStyle','None');
c=colorbar;
c.Label.String = 'Stress (MPa)';
colormap jet
title('Max Shear Stress (TaoMax in MPa)')
xlabel('X [m]')
ylabel('Z [m]')
grid on
axis tight


%% BENCHMARKING

cprintf('[210/255,0,0]          ','| |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| |\n')
cprintf('[210/255,0,0]          ','| |                    BENCHMARKING                  | |\n')
cprintf('[210/255,0,0]          ','| |               (DO NOT CHANGE VALUES)             | |\n')
cprintf('[210/255,0,0]          ','| |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| |\n')

dlmwrite('Data.txt',sortrows([vXnodes_Q(nodes_in_right) vYnodes_Q(nodes_in_right)],2),'delimiter','\t','precision',9) %write output
C=[0.0000 .2986 .6192 .7086 .6440 .5277 .4078 .2912 .1720 .0585 0.0000]; %Constant for a fix fix from book
x_l=([0.0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0]*3)+0.1; % y axis

F=(R*Density_fluid*g_a*H)/(2*pi); %qr in the book
N=C*F; % Hoop force from table
y=(N*R)/(E*wall);

Data=load('Data.txt');
figure('Name','Benchmark','NumberTitle','off')

plot((Data(:,1)-4.9)/2,Data(:,2),'Color',[210/255 0 0])
hold on
plot(y,x_l,'o','Color',[10/255,30/255,150/255])
grid on
title('Benchmark Test')
xlabel('X [m]')
ylabel('Z [m]')
legend('TANKEF FEA SOLUTION','GHALI TABLE')

%% END
cprintf('\n')
cprintf('[10/255,30/255,150/255]','| |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| |\n')
cprintf('[10/255,30/255,150/255]','| |                      THE END                     | |\n')
cprintf('[10/255,30/255,150/255]','| |--------------------------------------------------| |\n')
cprintf('[10/255,30/255,150/255]','| |--------------------------------------------------| |\n')
cprintf('\n')