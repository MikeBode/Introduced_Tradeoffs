clear all

%% === Only run the following section the first time you start to generate communities ===
% NumCommunities = zeros(1,50);
% for i = 1:50
%    Communities{i,1} = [];
% end
% save PersistentCommunities *ommunities
%% =======================================================================================

% Keep track of the number of successfully coexisting metacommunities (there are two types)
Success = 0;

% How many species do we want?
NumSpp = 15;

% Set the mortality rate (uniform across species)
m = 0.05;

% Set the number of coexisting metacommunities that we want to find
for R = 1:1000
   
   COEX = 0; Trial = 0;
   while COEX == 0;

      % Define a vector of colonisation abilities, ascending in strength
      c = sort(2*rand(NumSpp,1),'ascend');
      
      % Solve for the equilbrium abundance of the first (i.e., competitively-dominant) species
      Eq = zeros(NumSpp,1);
      Eq(1) = 1 - m/c(1);
      
      % Sequentially solve for the equilibrium abundance of subsequent species
      for n = 2:NumSpp
         Eq(n,1) = 1 - m/c(n) - sum(Eq(1:n-1).*(1 + c(1:n-1)/c(n)));
      end
      
      % Note that we've gone through another trial
      Trial = Trial + 1;
      
      % Display the progress occasionally
      if mod(Trial,10000) == 0;
         disp([Trial Success])
      end
      
      % If all species are coexisting
      if min(Eq) > 0;
         
         % Track the number of coexisting communities
         Success = Success + 1;
         
         % Terminate the while loop
         COEX = 1;
         
      end
   end
   
   % Append this community to the list of successfully coexisting communities
   load PersistentCommunities *ommunities
   NumCommunities(NumSpp) = NumCommunities(NumSpp) + 1;
   Communities{NumSpp,NumCommunities(NumSpp)} = c;
   save PersistentCommunities *ommunities
end








