% clear
clc
close all
fclose all;
delete('*.asv');

l=2;
CFL=.4;
c=.1;
t_end=l/c;
global epsilon
epsilon=0.03;

N_points=201;
Delta_x=(l-0)/(N_points-1);
x=0:Delta_x:l;
Delta_t=CFL*Delta_x/c;
N_timestep=ceil(t_end/Delta_t)+1;


phi_0=InitialValue(x);

phi_Exact(1,:)=phi_0;
phi_EE(1,:)=phi_0;
phi_RK2(1,:)=phi_0;
phi_LW(1,:)=phi_0;
phi_MLW(1,:)=phi_0;

figure1 = figure;
axes1 = axes('Parent',figure1);


for n=1:N_timestep-1
    
    phi_Exact(n+1,:) = ExactSolution(x, (n-1) * Delta_t,c);
%     phi_Exact(n+1,:) = InitialValue(x- c*(n-1) * Delta_t);

    phi_EE(n+1,:) = ExplicitEuler1stOBackward(phi_EE(n,:),Delta_t,Delta_t,CFL);
    phi_RK2(n+1,:) = RK2SecondOBackward(phi_RK2(n,:),Delta_t,Delta_t,CFL);
    phi_LW(n+1,:) = LaxWendroff(phi_LW(n,:),Delta_t,Delta_t,CFL);
    phi_MLW(n+1,:) = ModifiedLaxWendroff(phi_MLW(n,:),Delta_t,Delta_t,CFL);
    
    plot1 = plot(x,phi_Exact(n,:),'Parent',axes1,'Color',[0 0 0],...
        'DisplayName','Exact');
    hold(axes1,'on');
    plot2 = plot(x,phi_EE(n,:),'Parent',axes1,'Color',[0.85,0.33,0.10],...
        'DisplayName','EE-1stOB');
    plot3 = plot(x,phi_RK2(n,:),'Parent',axes1,'Color',[0.00,0.45,0.74],...
        'DisplayName','RK2-2ndOB');
    plot4 = plot(x,phi_LW(n,:),'Parent',axes1,'Color',[0.93,0.69,0.13],...
        'DisplayName','Lax-Wendroff');
    plot5 = plot(x,phi_MLW(n,:),'Parent',axes1,'Color',[0.49,0.18,0.56],...
        'DisplayName','Modified-LW');
    
    
    title(strcat('$t=',num2str((n-1)*Delta_t),' \mathrm{s}$'),'Interpreter','latex');
    ylabel({'$\phi(x)$'},'Interpreter','latex');
    xlabel({'$x$'},'Interpreter','latex');
    axis([0 l -1.5 1.5]);
    box(axes1,'on');
    set(axes1,'TickLabelInterpreter','latex');
    legend1 = legend([plot1, plot2, plot3, plot4, plot5]);
    set(legend1,'EdgeColor',[0 0 0],'Orientation','vertical','Location','Best',...
        'FontSize',9,'Interpreter','latex','FontName','Times New Roman');
    
    
    F(n) = getframe(figure1,[1 1 560 420]);
%         drawnow
    hold(axes1,'off');
end

n=N_timestep;

plot1 = plot(x,phi_Exact(n,:),'Parent',axes1,'Color',[0 0 0],...
    'DisplayName','Exact');
hold(axes1,'on');
plot2 = plot(x,phi_EE(n,:),'Parent',axes1,'Color',[0.85,0.33,0.10],...
    'DisplayName','EE-1stOB');
plot3 = plot(x,phi_RK2(n,:),'Parent',axes1,'Color',[0.00,0.45,0.74],...
    'DisplayName','RK2-2ndOB');
plot4 = plot(x,phi_LW(n,:),'Parent',axes1,'Color',[0.93,0.69,0.13],...
    'DisplayName','Lax-Wendroff');
plot5 = plot(x,phi_MLW(n,:),'Parent',axes1,'Color',[0.49,0.18,0.56],...
    'DisplayName','Modified-LW');


title(strcat('$t=',num2str((n-1)*Delta_t),' \mathrm{s}$'),'Interpreter','latex');
ylabel({'$\phi(x)$'},'Interpreter','latex');
xlabel({'$x$'},'Interpreter','latex');
axis([0 l -1.5 1.5]);
box(axes1,'on');
set(axes1,'TickLabelInterpreter','latex');
legend1 = legend([plot1, plot2, plot3, plot4, plot5]);
set(legend1,'EdgeColor',[0 0 0],'Orientation','vertical','Location','Best',...
    'FontSize',9,'Interpreter','latex','FontName','Times New Roman');


F(n) = getframe(figure1,[1 1 560 420]);
%     drawnow
hold(axes1,'off');


writerObj = VideoWriter('Wave.avi');
writerObj.FrameRate = floor(1/Delta_t);
open(writerObj);
writeVideo(writerObj, F);
close(writerObj);
