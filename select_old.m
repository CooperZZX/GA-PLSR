% SELECT.M          (universal SELECTion)
%
% This function performs universal selection. The function handles
% multiple populations and calls the low level selection function
% for the actual selection process.
%
% Syntax:  SelCh = select(SEL_F, Chrom, FitnV, GGAP, SUBPOP)
%
% Input parameters:
%    SEL_F     - Name of the selection function
%    Chrom     - Matrix containing the individuals (parents) of the current
%                population. Each row corresponds to one individual.
%    FitnV     - Column vector containing the fitness values of the
%                individuals in the population.
%    GGAP      - (optional) Rate of individuals to be selected
%                if omitted 1.0 is assumed
%    SUBPOP    - (optional) Number of subpopulations
%                if omitted 1 subpopulation is assumed
%
% Output parameters:
%    SelCh     - Matrix containing the selected individuals.
%
% Author:     Hartmut Pohlheim
% History:    10.03.94     file created
%             22.01.03     tested under MATLAB v6 by Alex Shenfield

function SelCh = select(SEL_F, Chrom, FitnV, GGAP, SUBPOP);

% Check parameter consistency
   if nargin < 3, error('Not enough input parameter'); end

   % Identify the population size (Nind)
   [NindCh,Nvar] = size(Chrom);
   [NindF,VarF] = size(FitnV);
   if NindCh ~= NindF, error('Chrom and FitnV disagree'); end
   if VarF ~= 1, error('FitnV must be a column vector'); end
  
   if nargin < 5, SUBPOP = 1; end
   if nargin > 4,
      if isempty(SUBPOP), SUBPOP = 1;
      elseif isnan(SUBPOP), SUBPOP = 1;
      elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar'); end
   end

   if (NindCh/SUBPOP) ~= fix(NindCh/SUBPOP), error('Chrom and SUBPOP disagree'); end %fix：删除每个数的小数部分
   Nind = NindCh/SUBPOP;  % Compute number of individuals per subpopulation

   if nargin < 4, GGAP = 1; end
   if nargin > 3,
      if isempty(GGAP), GGAP = 1;
      elseif isnan(GGAP), GGAP = 1;
      elseif length(GGAP) ~= 1, error('GGAP must be a scalar');
      elseif (GGAP < 0), error('GGAP must be a scalar bigger than 0'); end
   end

% Compute number of new individuals (to select)
   NSel=max(floor(Nind*GGAP+.5),2); %选择GGAP比例的染色体用于遗传的染色体

% Select individuals from population
   SelCh = [];
   for irun = 1:SUBPOP %输入参数<5，SUBPOP=1，相当于挑选前 FitnVSub=FitnV
      FitnVSub = FitnV((irun-1)*Nind+1:irun*Nind);
      ChrIx=feval(SEL_F, FitnVSub, NSel)+(irun-1)*Nind;%feval功能是调用函数体'SEL_F'进行计算 相当于从FitnVSub中用SEL_F方法挑选NSel个染色体
      SelCh=[SelCh; Chrom(ChrIx,:)];%选择第ChrIx行的染色体（Chrom的行数=染色体总数，列数=波段总数）
   end


% End of function