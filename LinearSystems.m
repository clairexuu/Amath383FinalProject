clear all;
close all;

%% Plotting Solutions of Linear System

A = [2 1; 1 2];
ode_RHS = @(t,y) A*[y(1) ; y(2)];
T = 1;
num_steps = 100;
h = T/num_steps;
init_cond = [1 ; -0.5];

[t_soln,y_soln] = ode45(ode_RHS,0:h:T,init_cond);

figure(1);
p1 = plot(t_soln,y_soln(:,1),'LineWidth',3,'Color',[0 0.6 0.7]);
hold on;
p2 = plot(t_soln,y_soln(:,2),'LineWidth',3,'Color',[0.6 0 0.7]);
set(gca,'FontSize',20);
xlim([0 T]);
ylim([-6 6]);
xticks(linspace(0,T,5));
grid on;
box on;
legend([p1 p2],{'P(t)','Q(t)'},'Location','SouthWest','FontSize',25);

%% Plotting Corresponding Trajectory in Phase Plane

figure(2);
plot([-4 ; 4],[0; 0],'k','LineWidth',1.5);
hold on;
plot([0 0],[-6 6],'k','LineWidth',1.5);
plot(y_soln(:,1),y_soln(:,2),'LineWidth',3,'Color',[0 0.6 0.7]);
set(gca,'FontSize',20);
xlabel('P(t)','FontSize',25);
ylabel('Q(t)','FontSize',25);
xlim([-4 4]);
ylim([-6 6]);
box on;
grid on;

%% Plotting Trajectory with Vector Field in Phase Plane

[X,Y] = meshgrid(-4:4,-6:6);

U = A(1,1)*X + A(1,2)*Y;
V = A(2,1)*X + A(2,2)*Y;

% Normalize the vectors to unit length

U_unit = U./sqrt(U.^2 + V.^2);
V_unit = V./sqrt(U.^2 + V.^2);

figure(3);
plot([-4 ; 4],[0; 0],'k','LineWidth',1.5);
hold on;
plot([0 0],[-6 6],'k','LineWidth',1.5);
quiver(X,Y,U_unit,V_unit,'off','Color',[0.6 0 0.7]);
plot(y_soln(:,1),y_soln(:,2),'LineWidth',3,'Color',[0 0.6 0.7]);
set(gca,'FontSize',20);
xlabel('P(t)','FontSize',25);
ylabel('Q(t)','FontSize',25);
xlim([-4 4]);
ylim([-6 6]);
box on;
grid on;

%% Depending on A, Trajectories/Vector Field in Phase Plane Are Different.

A = [0 1 ; 1 0];
ode_RHS = @(t,y) A*[y(1) ; y(2)];
T = 1;
num_steps = 100;
h = T/num_steps;
init_cond = [1 ; -0.5];

[t_soln,y_soln] = ode45(ode_RHS,0:h:T,init_cond);

[X,Y] = meshgrid(-4:4,-6:6);

U = A(1,1)*X + A(1,2)*Y;
V = A(2,1)*X + A(2,2)*Y;

U_unit = U./sqrt(U.^2 + V.^2);
V_unit = V./sqrt(U.^2 + V.^2);

figure(4);
plot([-4 ; 4],[0; 0],'k','LineWidth',1.5);
hold on;
plot([0 0],[-6 6],'k','LineWidth',1.5);
quiver(X,Y,U_unit,V_unit,'Color',[0.6 0 0.7]);
plot(y_soln(:,1),y_soln(:,2),'LineWidth',3,'Color',[0 0.6 0.7]);
set(gca,'FontSize',20);
xlabel('P(t)','FontSize',25);
ylabel('Q(t)','FontSize',25);
xlim([-4 4]);
ylim([-6 6]);
box on;

%% Depending on A, Trajectories/Vector Field in Phase Plane Are Different.

A = [0 1 ; -1 0];
ode_RHS = @(t,y) A*[y(1) ; y(2)];
T = 2*pi;
num_steps = 100;
h = T/num_steps;
init_cond = [1 ; -0.5];

[t_soln,y_soln] = ode45(ode_RHS,0:h:T,init_cond);

[X,Y] = meshgrid(-4:4,-6:6);

U = A(1,1)*X + A(1,2)*Y;
V = A(2,1)*X + A(2,2)*Y;

U_unit = U./sqrt(U.^2 + V.^2);
V_unit = V./sqrt(U.^2 + V.^2);

figure(5);
plot([-4 ; 4],[0; 0],'k','LineWidth',1.5);
hold on;
plot([0 0],[-6 6],'k','LineWidth',1.5);
quiver(X,Y,U_unit,V_unit,'Color',[0.6 0 0.7]);
plot(y_soln(:,1),y_soln(:,2),'LineWidth',3,'Color',[0 0.6 0.7]);
set(gca,'FontSize',20);
xlabel('P(t)','FontSize',25);
ylabel('Q(t)','FontSize',25);
xlim([-4 4]);
ylim([-6 6]);
box on;

%% Depending on A, Trajectories/Vector Field in Phase Plane Are Different.

A = [1 1 ; -1 0];
ode_RHS = @(t,y) A*[y(1) ; y(2)];
T = 2*pi;
num_steps = 100;
h = T/num_steps;
init_cond = [1 ; -0.5];

[t_soln,y_soln] = ode45(ode_RHS,0:h:T,init_cond);

[X,Y] = meshgrid(-4:4,-6:6);

U = A(1,1)*X + A(1,2)*Y;
V = A(2,1)*X + A(2,2)*Y;

U_unit = U./sqrt(U.^2 + V.^2);
V_unit = V./sqrt(U.^2 + V.^2);

figure(6);
plot([-4 ; 4],[0; 0],'k','LineWidth',1.5);
hold on;
plot([0 0],[-6 6],'k','LineWidth',1.5);
quiver(X,Y,U_unit,V_unit,'Color',[0.6 0 0.7]);
plot(y_soln(:,1),y_soln(:,2),'LineWidth',3,'Color',[0 0.6 0.7]);
set(gca,'FontSize',20);
xlabel('P(t)','FontSize',25);
ylabel('Q(t)','FontSize',25);
xlim([-4 4]);
ylim([-6 6]);
box on;

%% Depending on A, Trajectories/Vector Field in Phase Plane Are Different.

A = [-1 1 ; -1 0];
ode_RHS = @(t,y) A*[y(1) ; y(2)];
T = 4*pi;
num_steps = 100;
h = T/num_steps;
init_cond = [1 ; -0.5];

[t_soln,y_soln] = ode45(ode_RHS,0:h:T,init_cond);

[X,Y] = meshgrid(-4:4,-6:6);

U = A(1,1)*X + A(1,2)*Y;
V = A(2,1)*X + A(2,2)*Y;

U_unit = U./sqrt(U.^2 + V.^2);
V_unit = V./sqrt(U.^2 + V.^2);

figure(7);
plot([-4 ; 4],[0; 0],'k','LineWidth',1.5);
hold on;
plot([0 0],[-6 6],'k','LineWidth',1.5);
quiver(X,Y,U_unit,V_unit,'Color',[0.6 0 0.7]);
plot(y_soln(:,1),y_soln(:,2),'LineWidth',3,'Color',[0 0.6 0.7]);
set(gca,'FontSize',20);
xlabel('P(t)','FontSize',25);
ylabel('Q(t)','FontSize',25);
xlim([-4 4]);
ylim([-6 6]);
box on;

