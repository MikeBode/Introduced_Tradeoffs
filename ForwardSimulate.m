%% ======================================================
%% This function forward simulates the community
function [p,EqP,TVec] = ForwardSimulate(p,m,c,k,Tmax)

NG = nargin;
if NG == 4
   Tmax = 5000;
end
dt = 1; NumSpp = length(c);
for t = 1:Tmax./dt

   for s = 1:NumSpp
      p(s,t+1) = max(0, p(s,t) ...
                               + dt.*( (c(s).*p(s,t) + k(s) ).*(1 - sum(p(1:s,t))) ...
                               - m(s).*p(s,t) ...
                               - sum( c(1:s-1).*p(s,t).*p(1:s-1,t) )));
      
   end
   
   if mod(t,100) == 1
      Diff = sum(abs(p(:,t+1)-p(:,t)));
   end
   
   if NG == 4 && Diff < 1e-6
      break
   end
end
TVec = [0:dt:dt*t]; EqP = p(:,end);