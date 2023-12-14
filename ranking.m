% RANKING.M      (RANK-based fitness assignment)
%
% This function performs ranking of individuals.
%
% Syntax:  FitnV = ranking(ObjV, RFun, SUBPOP)
%
% This function ranks individuals represented by their associated
% cost, to be *minimized*, and returns a column vector FitnV
% containing the corresponding individual fitnesses. For multiple
% subpopulations the ranking is performed separately for each
% subpopulation.
%
% Input parameters:
%    ObjV      - Column vector containing the objective values of the
%                individuals in the current population (cost values).
%    RFun      - (optional) If RFun is a scalar in [1, 2] linear ranking is
%                assumed and the scalar indicates the selective pressure.
%                If RFun is a 2 element vector:
%                RFun(1): SP - scalar indicating the selective pressure
%                RFun(2): RM - ranking method
%                         RM = 0: linear ranking
%                         RM = 1: non-linear ranking
%                If RFun is a vector with length(Rfun) > 2 it contains
%                the fitness to be assigned to each rank. It should have
%                the same length as ObjV. Usually RFun is monotonously
%                increasing.
%                If RFun is omitted or NaN, linear ranking
%                and a selective pressure of 2 are assumed.
%    SUBPOP    - (optional) Number of subpopulations
%                if omitted or NaN, 1 subpopulation is assumed
%
% Output parameters:
%    FitnV     - Column vector containing the fitness values of the
%                individuals in the current population.
%                
%
% Author:     Hartmut Pohlheim (Carlos Fonseca)
% History:    01.03.94     non-linear ranking
%             10.03.94     multiple populations
%             21.01.03     updated for MATLAB v6 by Alex Shenfield

function FitnV = ranking(ObjV, RFun, SUBPOP);

% Identify the vector size (Nind)
   [Nind,ans] = size(ObjV);

   if nargin < 2, RFun = []; end  %这里传入参数只有一个Objv，即RMSEP%
   if nargin > 1, if isnan(RFun), RFun = []; end, end
   if prod(size(RFun)) == 2 %如果 RFun空矩阵，则size(RFun)=[0 0],prod之后=0。
       if RFun(2) == 1, NonLin = 1;
       elseif RFun(2) == 0, NonLin = 0; 
       else error('Parameter for ranking method must be 0 or 1'); end
       RFun = RFun(1);
       if isnan(RFun), RFun = 2; end
   elseif prod(size(RFun)) > 2,
       if prod(size(RFun)) ~= Nind, error('ObjV and RFun disagree'); end
   elseif prod(size(RFun)) < 2, NonLin = 0; %进入到这里
   end
   
   if nargin < 3, SUBPOP = 1; end
   if nargin > 2,
       if isempty(SUBPOP), SUBPOP = 1;
       elseif isnan(SUBPOP), SUBPOP = 1;
       elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar'); end
   end

   if (Nind/SUBPOP) ~= fix(Nind/SUBPOP), error('ObjV and SUBPOP disagree'); end
   Nind = Nind/SUBPOP;  % Compute number of individuals per subpopulation
   %Nind=20
% Check ranking function and use default values if necessary
   if isempty(RFun),
      % linear ranking with selective pressure 2
         RFun = 2*[0:Nind-1]'/(Nind-1); %RFun=(0:19)/19*2
   elseif prod(size(RFun)) == 1
      if NonLin == 1,
         % non-linear ranking
         if RFun(1) < 1, error('Selective pressure must be greater than 1');
         elseif RFun(1) > Nind-2, error('Selective pressure too big'); end
         Root1 = roots([RFun(1)-Nind [RFun(1)*ones(1,Nind-1)]]);
         RFun = (abs(Root1(1)) * ones(Nind,1)) .^ [(0:Nind-1)'];
         RFun = RFun / sum(RFun) * Nind;
      else
         % linear ranking with SP between 1 and 2
         if (RFun(1) < 1 | RFun(1) > 2),
            error('Selective pressure for linear ranking must be between 1 and 2');
         end
         RFun = 2-RFun + 2*(RFun-1)*[0:Nind-1]'/(Nind-1);
      end
   end;

   FitnV = [];

% loop over all subpopulations
for irun = 1:SUBPOP %只执行一次，irun=1
   % Copy objective values of actual subpopulation
      ObjVSub = ObjV((irun-1)*Nind+1:irun*Nind); %ObjSub=ObjV
   % Sort does not handle NaN values as required. So, find those...
      NaNix = isnan(ObjVSub); %全0
      Validix = find(~NaNix); %1:20
   % ... and sort only numeric values (smaller is better).
      [ans,ix] = sort(-ObjVSub(Validix));

   % Now build indexing vector assuming NaN are worse than numbers,
   % (including Inf!)...
      ix = [find(NaNix) ; Validix(ix)];
   % ... and obtain a sorted version of ObjV
      Sorted = ObjVSub(ix); %对ObjV进行降序排序

   % Assign fitness according to RFun.
      i = 1;
      FitnVSub = zeros(Nind,1);%Nind=20 %RFun=(0:19)/19*2
      for j = [find(Sorted(1:Nind-1) ~= Sorted(2:Nind)); Nind]' %找到ObjV降序排序后前19个和后19个对应不上的。因为rmsep总会不一样，所以这里应该是j=1:19
         FitnVSub(i:j) = sum(RFun(i:j)) * ones(j-i+1,1) / (j-i+1); %这一步i=j,最终FitnVSub=RFun
         i =j+1;
      end

   % Finally, return unsorted vector.
      [ans,uix] = sort(ix);
      FitnVSub = FitnVSub(uix);

   % Add FitnVSub to FitnV
      FitnV = [FitnV; FitnVSub];
end

% End of function