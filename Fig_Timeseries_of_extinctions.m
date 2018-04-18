function Fig_Timeseries_of_extincitons()

sub_CreateFIGURE_TEMP(1)

function [Finished,ThisComm] = sub_CreateFIGURE_TEMP(TargetSpp,ThisComm)

% Simulate the increase in colonisation ability of an invasive species
NumSpp = 15;
CoexThreshold = 1E-4; % This is the abundance we consider "persisting"

K_PET_VEC = [0.1 0.2 0.05];
K_PET_VEC = [0.025 0.05 0.1];

%% ========= First, load a pre-constructed dispersal vector ========= 
load PersistentCommunities *ommunities

if nargin == 1; 
   ThisComm = ceil(rand*NumCommunities(NumSpp));
end

c = Communities{NumSpp,ThisComm};
m = 0.05.*ones(NumSpp,1); % Natural mortality rates

%% ======== Solve for the species equil ========
EqP_0 = zeros(NumSpp,1);
EqP_0(1) = 1 - m(1)/c(1);
for n = 2:NumSpp
   EqP_0(n,1) = 1 - m(n)/c(n) - sum(EqP_0(1:n-1).*(1 + c(1:n-1)/c(n)));
end
%% =============================================

Y = [50 450];

% First perturb for Y1 years
k = zeros(NumSpp,1); k(TargetSpp) = K_PET_VEC(2);
[Perturb_p1,EqP,TV1] = ForwardSimulate(EqP_0,m,c,k,Y(1));

% Return to normal for Y2 years
k = zeros(NumSpp,1); 
[Perturb_p2,EqP,TV2] = ForwardSimulate(Perturb_p1(:,end),m,c,k,Y(2));

%% ========= Plot timeline of abundances ========= 
figure(1), clf; FS = 14; MS = 24; LW = 1.5; YL = 1.5e-6; YM = 1;
hold on; box on; set(gca,'linewidth',1.5)
Ptc = patch([0 Y(1) Y(1) 0],[YL YL YM YM],'k'); 
set(Ptc,'facealpha',0.3,'edgecolor','none')

% Plot all species
plot([-100 1],[EqP_0 EqP_0],'linewidth',LW,'color','b')
plot(TV1,Perturb_p1','linewidth',LW,'color','b')
plot(TV1(end)+TV2,Perturb_p2','linewidth',LW,'color','b')
% Plot the targeted invasive species
plot([-100 1],[EqP_0(TargetSpp) EqP_0(TargetSpp)],'linewidth',LW+0.5,'color','r')
plot(TV1,Perturb_p1(TargetSpp,:),'linewidth',LW+0.5,'color','r')
plot(TV1(end)+TV2,Perturb_p2(TargetSpp,:),'linewidth',LW+0.5,'color','r')

set(gca,'yscale','log','xtick',[-100:100:400],'yminortick','off')
ylim([YL YM]); xlim([-50 350])
xlabel('Time','fontsize',FS)
ylabel('Proportional abundance','fontsize',FS)

% % Print the figure
% set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 3 22 14])
% print('-dpdf',['Figures/TEMP_ExtinctionTimeline_c' num2str(ThisComm) '_s' num2str(TargetSpp) '.pdf'])
















