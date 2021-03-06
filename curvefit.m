function [fhandle,exitflag,message] = curvefit(points,nCurves,curveOrder,varargin)
  % Help for curvefit function
  % Inputs: points (nP,2) [y,x], point to fit with curve
  %       : nCurves, number of curves 
  %       : curveOrder, curve polynomial order
  % Options
  %
  % Outputs
  % fhandle: Function handle for evaluating the curve fit
  % exitflag: > 0 everything is ok
  %           1: Curvefit sucessfull
  %           0: Input error
  % message: exit message from algorithm

  fhandle = [];
  exitflag = 1;
  message = '';
  
  % check points structure
  if  size(points,2) ~= 2
    exitflag = 0;
    message = sprintf('curvefit: input points needs to be a nP x 2 array. Number of columns for input is %3i', size(points,2));
    return
  end
  
  % get options strucutre containing user defined options
  [options,exitflag,message] = setOptions(points,nCurves,curveOrder,varargin);
  if exitflag < 1
    return
  end
  
  % get problem structure containing information nessary to build optimization problem
  [prob,exitflag,message]  = problemformulation(options);
  if exitflag < 1
    return
  end
  
  % Assemble linear constraints
  [A,b,Aeq,beq] = getLinearConstraints(prob,options);
  
  % Specify function handles
  fun = @(x) getObjectiveFunction(x,prob);
  nonlcon = @(x) getNonLinearConstraints(x,prob,options);
  % Call optimizer(s)
  opti = fminslp(fun,prob.x0,A,b,Aeq,beq,prob.xL,prob.xU,nonlcon,'display','iter','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',true,'solver','linprog','MoveLimit',0.5);
  [xval,fval,exitflag,output] = opti.solve;
  message = output.message;
%  if exitflag < 1
%    return
%  end
  
  % Make minimal structure for function evaluation, save memory
  settings = struct('curveStart',[],...
                    'curveEnd',[],...
                    'nCurves',[],...
                    'curveNo2DVNo',[],...
                    'curveOrder',[],...
                    'cuvedfCof',[]);
  settings.curveStart = prob.curveStart;
  settings.curveEnd = prob.curveEnd;
  settings.nCurves = prob.nCurves;
  settings.curveNo2DVNo = prob.curveNo2DVNo;
  settings.curveOrder = prob.curveOrder;
  settings.cuvedfCof = prob.cuvedfCof;
  
  if options.floating
    fDVNo = prob.curveNo2FloatDVNo(1,2);
    settings.curveEnd(1) = xval(fDVNo);
    for curveNo = 2:prob.nCurves-1
      fDVNo = prob.curveNo2FloatDVNo(curveNo,1);
      settings.curveStart(curveNo) = xval(fDVNo);
      fDVNo = prob.curveNo2FloatDVNo(curveNo,2);
      settings.curveEnd(curveNo) = xval(fDVNo);
    end
    fDVNo = prob.curveNo2FloatDVNo(end,1);
    settings.curveStart(end) = xval(fDVNo);
  end
  
  
  % Make output function handle for evaluation of curve(s)
  fhandle = @(x,varargin) evalFit(x,settings,xval,varargin);
  
  % plot result or not
  if options.plot
    plotFit(xval,settings,points);
  end
  
  
end

function [options,exitflag,message] = setOptions(points,nCurves,curveOrder,input)
  % Initialize output variables
  exitflag = 1;
  message = [];
  
  % Here you can add new options if needed
  p = inputParser;
  p.CaseSensitive = false;
  % Helper functions for input parser
  checkEmpetyOrChar = @(x)(isempty(x) || ischar(x));
  checkEmptyOrNumericPositive = @(x)(isempty(x) || ((isnumeric(x) || isscalar(x)) && all(x > 0)));
  checkScalarNumPos = @(x)(isnumeric(x) || isscalar(x)) && (x > 0);
  checkScalarNumGTeqZero = @(x)(isnumeric(x) || isscalar(x)) && (x >=0 );
  checkScalarNum = @(x)(isnumeric(x) || isscalar(x));
  checkEmptyOrScalarNum = @(x)(isempty(x) || isnumeric(x) || isscalar(x));
  checkPoints =@(x) isnumeric(x) && size(x,2)==2 && size(x,1)>1;

  % Settings specifying the curve(s)

  % The order of the required inputs matter
  p.addRequired('points',checkPoints);
  p.addRequired('nCurves',checkScalarNumPos);
  p.addRequired('curveOrder',checkScalarNumPos);
  % 
  p.addParameter('startpos',[],checkEmptyOrScalarNum);
  p.addParameter('endpos',[], @(x)checkEmptyOrScalarNum(x));
  p.addParameter('curveContinuity',max(curveOrder-1,0), @(x)checkScalarNumGTeqZero(x));

  % Settings for floating curve intersection points
  p.addParameter('floating',false,  @(x)islogical(x));
  p.addParameter('beta',100,  @(x)checkScalarNumPos(x));

  % General settings
  p.addParameter('display','off',  @(x)checkEmpetyOrChar(x));
  p.addParameter('plot',false,  @(x)islogical(x));

  % Additional constraints
  p.addParameter('c0',[],  @(x)checkEmptyOrScalarNum(x) && size(x,2)==2 && size(x,1)>0); % local c0
  p.addParameter('c1',[],  @(x)checkEmptyOrScalarNum(x) && size(x,2)==2 && size(x,1)>0); % local c1
  p.addParameter('pointlb',[],  @(x)checkEmptyOrScalarNum(x)); % for each point, we can specify a lower bound value, curve must be above this value
  p.addParameter('pointub',[],  @(x)checkEmptyOrScalarNum(x)); % for each point, we can specify an upper bound value, curve must be below this value
  p.addParameter('kappa',[],  @(x)checkEmptyOrNumericPositive(x)); % local curvature constraints, the curve must not have an abs(geometric curvature) above this value

  % Solver settings
  p.addParameter('method','bound',  @(x)checkEmpetyOrChar(x));
  p.addParameter('solver','fminslp',  @(x)checkEmpetyOrChar(x));

  % pars input
  if nargin < 4 || isempty(input)
      parse(p,points,nCurves,curveOrder);
  else
      parse(p,points,nCurves,curveOrder,input{:});
  end

  % Output results to options structure
  options = p.Results;
  
  if isempty(options.startpos)
    options.startpos = min(points(:,2));
  end
  
  % Determine if curve order has been specified for each curve or just one value which should apply to all curves
  if numel(options.curveOrder) ~= options.nCurves
    if numel(options.curveOrder) == 1
      options.curveOrder = repmat(options.curveOrder,[options.nCurves,1]);
    else
       message = sprintf('Please only specify curve order either by one value, or for each curve');
       exitflag = 0;
    end
  end
  
    % Determine if curve continuity has been specified for each curve segment or just one value which should apply to all curve segments
  if numel(options.curveContinuity) ~= options.nCurves-1
    if numel(options.curveContinuity) == 1
      options.curveContinuity = repmat(options.curveContinuity,[options.nCurves-1,1]);
    else
       message = sprintf('Please only specify curve continuity either by one value, or for each curve segment');
       exitflag = 0;
    end
  end
  
  % Determine number of points in the fit
  options.nP = size(options.points,1);

end

function [prob,exitflag,message] = problemformulation(options)
  % Initialize output variables
  exitflag = 1;
  message = [];

  % Load some key data from options structure to prob structure
  prob.nCurves = options.nCurves;
  prob.curveOrder = options.curveOrder;
  prob.curveContinuity = options.curveContinuity;
  prob.method = options.method;
  
  % Determine how the user has specified the problem
  if numel(options.startpos) == 1
    % distribute start positions for the curves between the specified start and end position
    if isempty(options.endpos)
      options.endpos = max(options.points(:,2));
    elseif numel(options.endpos) > 1
      message = sprintf('Please only specify one end position');
      exitflag = 0;
      return
    end

    if options.startpos > options.endpos
      message = sprintf('Specified start position is greater than specified end position. %3.4f > %3.4f', options.startpos,options.endpos);
      exitflag = 0;
      return
    end
    temp = linspace(options.startpos,options.endpos,options.nCurves+1);
    prob.curveStart = temp(1:end-1);
    prob.curveEnd = temp(2:end);
  else
    if numel(options.startpos) ~= options.nCurves
      message = sprintf('Please specify start positions for each curve interval');
      exitflag = 0;
      return
    end
    for cvNo = 1:options.nCurves-1
      if options.startpos(cvNo) >= options.startpos(cvNo+1)
        message = sprintf('Specified start positions must come in increasing order. startPos(%i) >= startPos(%i)',options.startpos(cvNo),options.startpos(cvNo+1));
        exitflag = 0;
        return
      end
    end
    prob.curveStart = options.startpos(:);
    prob.curveEnd = [prob.curveStart(2:end);options.endpos];
  end
  
  maxOrder = max(prob.curveOrder+1); % We need +1 when having floating start and end positions
  prob.cuvedfCof = zeros(maxOrder,maxOrder+1);
  
  for dfNo = 1:maxOrder
    for orderNo = dfNo:maxOrder
      % calculate orderNo*(orderNo-1)*(orderNo-2)*...*(orderNo-n)
      devPart = orderNo;
      for ii = 1:dfNo-1
        devPart = devPart*(orderNo-ii);
      end
      prob.cuvedfCof(dfNo,orderNo+1) = devPart;
    end
  end
  
  
  % Initialize counters for design variables and constraints
  prob.nCurveDV = 0; % Number of design variables related to curve formulation
  prob.nBoundDV = 0; % Number of bound design variables, one for each curve segment
  prob.nFloatDV = 0; % Number of floating start positions, one for each curve
  prob.nDV      = 0; % Total number of design variables
  prob.nA       = 0; % Number of linear in-equality constraints
  prob.nAeq     = 0; % Number of linear equality constraints
  prob.nC       = 0; % Number of non-linear in-equality constraints
  prob.nCeq     = 0; % Number of non-linear equality constraints
  prob.nG       = 0; % Total number of constraints
  
  if options.floating
    [prob,exitflag,message] = setupFloatingFormulation(prob,options);
  else
    [prob,exitflag,message] = setupStandardFormulation(prob,options);
  end % not floating formulation
  
  % Total number of design variables
  prob.nDV = prob.nCurveDV + prob.nBoundDV + prob.nFloatDV;
 
  % Total number of constraints
  prob.nG = prob.nA + prob.nAeq + prob.nC + prob.nCeq;
  
  % Allocate design variable data structures
  prob.xL = repmat(-10000,[prob.nDV,1]); % Initial lower bound for design variables
  prob.xU = repmat(10000,[prob.nDV,1]); % Initial upper bound for design variables
  if strcmpi(options.method,'bound')
    % Bookkeeping to get bound design variable numbers wrt., point numbers
    prob.pointNo2BoundDVNo = (prob.nCurveDV+1):(prob.nCurveDV+prob.nBoundDV);
    % Modify upper and lower bounds for bound design variables
    prob.xL(prob.pointNo2BoundDVNo) = 0;
    prob.xU(prob.pointNo2BoundDVNo) = 1e6;
  end
  % Initialize design variables
  prob.x0 = zeros(prob.nDV,1);
  
  if options.floating
    for curveNo = 2:prob.nCurves
      fDVNo = prob.curveNo2FloatDVNo(curveNo,1);
      prob.xL(fDVNo) = prob.xa;
      prob.xU(fDVNo) = prob.xb;
      prob.x0(fDVNo) = prob.curveStart(curveNo);
    end
  end
  
  
end

function [prob,exitflag,message] = setupStandardFormulation(prob,options)
  exitflag = 1;
  message = '';
  
  % Setup bookkeeping between points and curves
  prob.nPointsPerCurve = zeros(prob.nCurves,1);
  prob.pointNo2CurveNo = zeros(options.nP,1);
  
  % Temp array to keep track of which points have been assigned to curves
  % This is nessary when points lie exactly between two curves start and end positions
  assigned = false(options.nP,1);
  for curveNo = 1:prob.nCurves
    for pNo = 1:options.nP
      pos = options.points(pNo,2);
      if pos >= prob.curveStart(curveNo) && pos <= prob.curveEnd(curveNo) && ~assigned(pNo)
        prob.nPointsPerCurve(curveNo) = prob.nPointsPerCurve(curveNo) + 1;
        assigned(pNo) = true;
      end
    end
  end
  
  % check if we have assigned all points to curves
  if ~all(assigned)
    nAssigned = sum(assigned);
    message = sprintf('Not all points have been assigned to a curve, this should not be possible. \n Total number of points assigned as %i, and total number of points is %i',nAssigned,options.nP);
    exitflag = 0;
    return
  end
  
  % Allocate data structure containing point ids for each curve
  for curveNo = 1:prob.nCurves
    prob.curvePoints{curveNo} = zeros(prob.nPointsPerCurve(curveNo),1);
  end
  
  % Load data into allocated cell array
  assigned = false(options.nP,1);
  for curveNo = 1:prob.nCurves
    iCount = 0;
    for pNo = 1:options.nP
      pos = options.points(pNo,2);
      if pos >= prob.curveStart(curveNo) && pos <= prob.curveEnd(curveNo) && ~assigned(pNo)
        iCount = iCount + 1;
        prob.curvePoints{curveNo}(iCount) = pNo;
        prob.pointNo2CurveNo(pNo) = curveNo;
        assigned(pNo) = true;
      end
    end
  end
  
  % Count number of curve design variables i.e., coefficients
  % In the design variable system/bookkeeping, these come first.
  prob.curveNo2DVNo = cell(prob.nCurves,1);
  for curveNo = 1:prob.nCurves
    prob.curveNo2DVNo{curveNo} = (prob.nCurveDV+1):(prob.curveOrder(curveNo) + prob.nCurveDV + 1);
    prob.nCurveDV = prob.nCurveDV + prob.curveOrder(curveNo) + 1;
  end
  
  if strcmpi(options.method,'bound')
    % Count number of bound design variables. One for each curve segment
    prob.nBoundDV = options.nP;
    % Count number of bound constraints. These are non-linear in-equality constraints.
    % We have to bound from below and above i.e., two constraints pr. point.
    prob.nC = prob.nC + sum(prob.nPointsPerCurve)*2; 
  end
  
  % Count number of continuity constraints between each curve segment
  % These are linear equality
  prob.nAeq = prob.nAeq + sum(options.curveContinuity+1); % plus 1 due to c^0 part
  
  % Count number of local c0 constraints, these are linear equality
  if ~isempty(options.c0)
    prob.nAeq = prob.nAeq + size(options.c0,1);
  end
  
  % Count number of local c1 constraints, these are linear equality
  if ~isempty(options.c1)
    prob.nAeq = prob.nAeq + size(options.c1,1);
  end
  
  % Count number of lower bound point constraints
  % These are linear in-equality constraints
  if ~isempty(options.pointlb)
    if options.nP ~= size(options.pointlb,1)
      message = sprintf('Number of lower bound point constraints(%i) is not equal to the number of points specified in the fit(%i)',size(options.pointlb,1),options.nP);
      exitflag = 0;
      return
    else
      prob.nA = prob.nA + options.nP;
    end
  end
  
  % Count number of upper bound point constraints
  % These are linear in-equality constraints
  if ~isempty(options.pointub)
    if options.nP ~= size(options.pointub,1)
      message = sprintf('Number of upper bound point constraints(%i) is not equal to the number of points specified in the fit(%i)',size(options.pointub,1),options.nP);
      exitflag = 0;
      return
    else
      prob.nA = prob.nA + options.nP;
    end
  end
  
end

function [prob,exitflag,message] = setupFloatingFormulation(prob,options)
  exitflag = 1;
  message = '';
  
  % Count number of curve design variables i.e., coefficients
  % In the design variable system/bookkeeping, these come first.
  prob.curveNo2DVNo = cell(prob.nCurves,1);
  for curveNo = 1:prob.nCurves
    prob.curveNo2DVNo{curveNo} = (prob.nCurveDV+1):(prob.curveOrder(curveNo) + prob.nCurveDV + 1);
    prob.nCurveDV = prob.nCurveDV + prob.curveOrder(curveNo) + 1;
  end
  
  if strcmpi(options.method,'bound')
    % Count number of bound design variables. One for each curve segment
    prob.nBoundDV = options.nP;
    % Count number of bound constraints. These are non-linear in-equality constraints.
    % We have to bound from below and above i.e., two constraints pr. point.
    % For the floating formulation, each curve has to be "presented" with all points. Consequently, we get a-lot of constraints
    prob.nC = prob.nC + prob.nCurves*options.nP*2; 
  end
  
  % Specify how many floating / curve start position design variables 
  prob.nFloatDV = prob.nCurves-1; % each curve segment has two variables, for the first segment, the start dv is fixed, for the last segment, the end dv is fixed
  prob.floatNo2DVNo = (prob.nCurveDV + prob.nBoundDV+1:prob.nCurveDV + prob.nBoundDV+prob.nFloatDV)';
  
  % We need to normalize the input x-coordinate between 0-1. We will call this normalized coordinate s(x)
  prob.xa = min(options.points(:,2));
  prob.xb = max(options.points(:,2));
  prob.xc = prob.xb-prob.xa;
  % make function for easy conversion between normalized s(x) to x(s) and vice-versa
  prob.xs = @(s) s.*prob.xc+prob.xa;
  prob.sx = @(x) (x-prob.xa)./prob.xc;
  prob.dsdx = 1/prob.xc;
  % Normalize all point coordinates
  prob.sxPoints = prob.sx(options.points(:,2));
  
  % get float dv no from curve no, the first curve is assigned a variable, but it will be fixed through it's bounds
  prob.curveNo2FloatDVNo = zeros(prob.nCurves,2);
  prob.curveNo2FloatDVNo(1,1) = 0; % first segment has no start variable
  prob.curveNo2FloatDVNo(1,2) = prob.floatNo2DVNo(1); % first segment, end variable
  for curveNo = 2:prob.nCurves-1
    prob.curveNo2FloatDVNo(curveNo,1) = prob.curveNo2FloatDVNo(curveNo-1,2);
    prob.curveNo2FloatDVNo(curveNo,2) = prob.curveNo2FloatDVNo(curveNo,1)+1;
  end
  % Last curve only has a start variable
  prob.curveNo2FloatDVNo(prob.nCurves,1) = prob.curveNo2FloatDVNo(prob.nCurves-1,2);
  
  % We need to ensure that curve start positions don't "overlap/ cross" i.e., pos1 <= pos2,...,posn-1 <= posn
  prob.nA = prob.nA + prob.nFloatDV-1;
  
  % Count number of continuity constraints between each curve segment
  % These are non-linear equality
  for curveNo = 1:prob.nCurves-1
    % Update constraint counter
     prob.nCeq =  prob.nCeq + 1;
    % Specify continuity constraints from c^1 to c^n: d^nf1(x)/dx^n - d^nf2(x)/dx^n = 0 
    for cc = 1:prob.curveContinuity(curveNo)
      % Update constraint counter
       prob.nCeq =  prob.nCeq + 1;
    end % curveContinuity
  end % nCurves
  
  % Count number of local c0 constraints, these are non-linear equality
  if ~isempty(options.c0)
    prob.nCeq = prob.nCeq + size(options.c0,1);
  end
  
  % Count number of local c1 constraints, these are non-linear equality
  if ~isempty(options.c1)
    prob.nCeq = prob.nCeq + size(options.c1,1);
  end
  
  % Count number of lower bound point constraints
  % These are non-linear in-equality constraints
  if ~isempty(options.pointlb)
    if options.nP ~= size(options.pointlb,1)
      message = sprintf('Number of lower bound point constraints(%i) is not equal to the number of points specified in the fit(%i)',size(options.pointlb,1),options.nP);
      exitflag = 0;
      return
    else
      prob.nC = prob.nC + options.nP;
    end
  end
  
  % Count number of upper bound point constraints
  % These are non-linear in-equality constraints
  if ~isempty(options.pointub)
    if options.nP ~= size(options.pointub,1)
      message = sprintf('Number of upper bound point constraints(%i) is not equal to the number of points specified in the fit(%i)',size(options.pointub,1),options.nP);
      exitflag = 0;
      return
    else
      prob.nC = prob.nC + options.nP;
    end
  end
  
end

function [fval,df] = getObjectiveFunction(xval,prob)
  % Main function for evaluating the objective function
  switch prob.method
    case 'bound'
      fval = sum(xval(prob.pointNo2BoundDVNo));
      if nargout > 1
        df = zeros(prob.nDV,1);
        df(prob.pointNo2BoundDVNo) = 1;
      end
  end
end


function [c,ceq,dc,dceq] = getNonLinearConstraints(xval,prob,options)
   % Main function for evaluating the non-linear constraints
  if prob.nC > 0
    c = zeros(prob.nC,1);
  else
    c = [];
  end
  
  if prob.nCeq > 0
    ceq = zeros(prob.nCeq,1);
  else
    ceq = [];
  end
  
  % If no non-linear constraints have been specified, then we exit the function
  if isempty(ceq) && isempty(c)
    dc = [];
    dceq = [];
    return
  end
  
  % Initialize counters
  cNo = 0;
  ceqNo = 0;
  
  if strcmpi(options.method,'bound') && ~options.floating
    % Evaluate bound constraints
    for curveNo = 1:prob.nCurves
      for No = 1:prob.nPointsPerCurve(curveNo)
        pNo = prob.curvePoints{curveNo}(No);
        x = options.points(pNo,2);
        y = options.points(pNo,1);
        
        % Extract first coefficient in polynomial
        f1 = xval(prob.curveNo2DVNo{curveNo}(1));
        for orderNo = 1:prob.curveOrder(curveNo)
          f1 = f1 + xval(prob.curveNo2DVNo{curveNo}(1+orderNo))*x^orderNo;
        end
        
        % update constraint counter
        cNo = cNo + 1;
        % f(x) - target - bound <= 0
        c(cNo) = (f1-y) - xval(prob.pointNo2BoundDVNo(pNo));
        
        % update constraint counter
        cNo = cNo + 1;
        % target - f(x) - bound <= 0
        c(cNo) = (y-f1) - xval(prob.pointNo2BoundDVNo(pNo));
      end
    end
  elseif strcmpi(options.method,'bound') && options.floating
    % Evaluate bound constraints
    for curveNo = 1:prob.nCurves
      if curveNo ~= 1 && curveNo~=prob.nCurves
        x1 = xval(prob.curveNo2FloatDVNo(curveNo,1));
        x2 = xval(prob.curveNo2FloatDVNo(curveNo,2));
        sx1 = prob.sx(x1);
        sx2 = prob.sx(x2);
        pwf = getPointWeightsForCurveInterval(prob.sxPoints,sx1,sx2,options.beta); % get point weight factor
      elseif curveNo == 1
        x2 = xval(prob.curveNo2FloatDVNo(curveNo,2));
        sx2 = prob.sx(x2);
        pwf = getFixedFloatingStartPointWeightsForCurveInterval(prob.sxPoints,sx2,options.beta); % get point weight factor
      elseif curveNo==prob.nCurves
        x1 = xval(prob.curveNo2FloatDVNo(curveNo,1));
        sx1 = prob.sx(x1);
        pwf = getFixedFloatingEndPointWeightsForCurveInterval(prob.sxPoints,sx1,options.beta); % get point weight factor
      end
      
      for pNo = 1:options.nP
        x = options.points(pNo,2);
        y = options.points(pNo,1);
        
        % Extract first coefficient in polynomial
        f1 = xval(prob.curveNo2DVNo{curveNo}(1));
        for orderNo = 1:prob.curveOrder(curveNo)
          f1 = f1 + xval(prob.curveNo2DVNo{curveNo}(1+orderNo))*x^orderNo;
        end
        
        % update constraint counter
        cNo = cNo + 1;
        % f(x) - target - bound <= 0
        c(cNo) = (f1-y)*pwf(pNo) - xval(prob.pointNo2BoundDVNo(pNo));
        
        % update constraint counter
        cNo = cNo + 1;
        % target - f(x) - bound <= 0
        c(cNo) = (y-f1)*pwf(pNo) - xval(prob.pointNo2BoundDVNo(pNo));
      end
    end
  end
  
  if options.floating
    for curveNo = 1:prob.nCurves-1
      % Extract end position of the current curve
      fDVNo = prob.curveNo2FloatDVNo(curveNo,2);
      x = xval(fDVNo);
      
      % Specify c^0 continuity: f1(x) - f2(x) = 0
      % f1 part
      % Extract first design variable for curve
      DVNo = prob.curveNo2DVNo{curveNo}(1);
      f1 = xval(DVNo);
      for orderNo = 1:prob.curveOrder(curveNo)
        DVNo = prob.curveNo2DVNo{curveNo}(1+orderNo);
        f1 = f1 + xval(DVNo)*x^(orderNo);
      end
      
      % f2 part
      % Extract first design variable from the curve infront of current curve
      DVNo = prob.curveNo2DVNo{curveNo+1}(1);
      f2 = xval(DVNo);
      for orderNo = 1:prob.curveOrder(curveNo+1)
        DVNo = prob.curveNo2DVNo{curveNo+1}(1+orderNo);
        f2 = f2 + xval(DVNo)*x^(orderNo);
      end
     
      % Update constraint counter
      ceqNo = ceqNo + 1;
      ceq(ceqNo) = (f1-f2);
     
      % Specify continuity constraints from c^1 to c^n: d^nf1(x)/dx^n - d^nf2(x)/dx^n = 0 
      for cc = 1:prob.curveContinuity(curveNo)
        % f1 part
        f1 = 0;
        for orderNo = cc:prob.curveOrder(curveNo)
          DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
          f1 = f1 + prob.cuvedfCof(cc,orderNo+1)*xval(DVNo)*x^(orderNo-cc);
        end
        
        % f2 part  
        f2 = 0;
        for orderNo = cc:prob.curveOrder(curveNo+1)
          DVNo = prob.curveNo2DVNo{curveNo+1}(orderNo+1);
          f2 = f2 + prob.cuvedfCof(cc,orderNo+1)*xval(DVNo)*x^(orderNo-cc);
        end
        
        % Update constraint counter
        ceqNo = ceqNo + 1;
        ceq(ceqNo) = (f1-f2);
        
      end % curveContinuity
    end % nCurves
  end
  
  % Error checks
  if cNo ~= prob.nC
    err_msg = sprintf('Total number of non-linear in-equality constraints(%i) does not match specified(%i)',cNo,prob.nC);
    error(err_msg);
  end
  
  if ceqNo ~= prob.nCeq
    err_msg = sprintf('Total number of non-linear equality constraints(%i) does not match specified(%i)',ceqNo,prob.nCeq);
    error(err_msg);
  end
  
  % DSA
  if nargout > 2
    [dc,dceq] = getNonLinearConstraintsDSA(xval,c,ceq,prob,options);
  end
  
  
end

function [dc,dceq] = getNonLinearConstraintsDSA(xval,c,ceq,prob,options)
   % Main function for evaluating the gradients of the non-linear constraints
  if prob.nC > 0
    dc = sparse(prob.nDV,prob.nC);
  else
    dc = [];
  end
  
  if prob.nCeq > 0
    dceq = sparse(prob.nDV,prob.nCeq);
  else
    dceq = [];
  end
  
  % Initialize counters
  cNo = 0;
  ceqNo = 0;
  
  if strcmpi(options.method,'bound') && ~options.floating
    % Evaluate bound constraints
    for curveNo = 1:prob.nCurves
      for No = 1:prob.nPointsPerCurve(curveNo)
        pNo = prob.curvePoints{curveNo}(No);
        x = options.points(pNo,2);
        y = options.points(pNo,1);
        
        % update constraint counter for the first constraint
        cNo = cNo + 1;
        % Set the first design variable
        dc(prob.curveNo2DVNo{curveNo}(1),cNo) = 1;
        dc(prob.curveNo2DVNo{curveNo}(1),cNo+1) = -1;
        % loop for the remaining design variables associated with the current curveNo
        for orderNo = 1:prob.curveOrder(curveNo)
          dc(prob.curveNo2DVNo{curveNo}(1+orderNo),cNo) = x^orderNo;
          dc(prob.curveNo2DVNo{curveNo}(1+orderNo),cNo+1) = -x^orderNo;
        end
        % add gradient for the bound variable
        dc(prob.pointNo2BoundDVNo(pNo),cNo) = -1;
        dc(prob.pointNo2BoundDVNo(pNo),cNo+1) = -1;
        % update constraint counter for the second constraint
        cNo = cNo + 1;
      end
    end
  elseif strcmpi(options.method,'bound') && options.floating
    % Evaluate bound constraints
    for curveNo = 1:prob.nCurves
     if curveNo ~= 1 && curveNo~=prob.nCurves
       x1 = xval(prob.curveNo2FloatDVNo(curveNo,1));
       x2 = xval(prob.curveNo2FloatDVNo(curveNo,2));
       sx1 = prob.sx(x1);
       sx2 = prob.sx(x2);
       pwf = getPointWeightsForCurveInterval(prob.sxPoints,sx1,sx2,options.beta); % get point weight factor
       dpwfs = getFloatingStartPointDSA(prob.sxPoints,sx1,sx2,options.beta).*prob.dsdx; % get DSA of sx1 wrt., all points
       dpwfe = getFloatingEndPointDSA(prob.sxPoints,sx1,sx2,options.beta).*prob.dsdx; % get DSA of sx2 wrt., all points
     elseif curveNo == 1
       x2 = xval(prob.curveNo2FloatDVNo(curveNo,2));
       sx2 = prob.sx(x2);
       pwf = getFixedFloatingStartPointWeightsForCurveInterval(prob.sxPoints,sx2,options.beta); % get point weight factor
       dpwfs = getFixedFloatingStartPointWeightsForCurveIntervalDSA(prob.sxPoints,sx2,options.beta).*prob.dsdx; % get DSA of sx1 wrt., all points
       dpwfe = getFloatingEndPointDSA(prob.sxPoints,0,sx2,options.beta).*prob.dsdx; % get DSA of sx2 wrt., all points
     elseif curveNo==prob.nCurves
       x1 = xval(prob.curveNo2FloatDVNo(curveNo,1));
       sx1 = prob.sx(x1);
       pwf = getFixedFloatingEndPointWeightsForCurveInterval(prob.sxPoints,sx1,options.beta); % get point weight factor
       dpwfs = getFloatingStartPointDSA(prob.sxPoints,sx1,1,options.beta).*prob.dsdx; % get DSA of sx1 wrt., all points
       dpwfe = getFixedFloatingEndPointWeightsForCurveIntervalDSA(prob.sxPoints,sx1,options.beta).*prob.dsdx; % get DSA of sx2 wrt., all points
     end
     
      for pNo = 1:options.nP
        x = options.points(pNo,2);
        y = options.points(pNo,1);
        
        % update constraint counter for the first constraint
        cNo = cNo + 1;
        % Set the first design variable
        dc(prob.curveNo2DVNo{curveNo}(1),cNo) = 1*pwf(pNo);
        dc(prob.curveNo2DVNo{curveNo}(1),cNo+1) = -1*pwf(pNo);
        % loop for the remaining design variables associated with the current curveNo
        for orderNo = 1:prob.curveOrder(curveNo)
          dc(prob.curveNo2DVNo{curveNo}(1+orderNo),cNo) = (x^orderNo)*pwf(pNo);
          dc(prob.curveNo2DVNo{curveNo}(1+orderNo),cNo+1) = -(x^orderNo)*pwf(pNo);
        end
        % add gradient for the bound variable xval()
        dc(prob.pointNo2BoundDVNo(pNo),cNo) = -1;
        dc(prob.pointNo2BoundDVNo(pNo),cNo+1) = -1;
        % add gradient from floating start point
        if curveNo ~= 1
          dc(prob.curveNo2FloatDVNo(curveNo,1),cNo) = (c(cNo)+xval(prob.pointNo2BoundDVNo(pNo)))*dpwfs(pNo); % we have to add the bound variable, to get the "real" constraint value
          dc(prob.curveNo2FloatDVNo(curveNo,1),cNo+1) = (c(cNo+1)+xval(prob.pointNo2BoundDVNo(pNo)))*dpwfs(pNo); % we have to add the bound variable, to get the "real" constraint value
        end
        if curveNo ~= prob.nCurves
          % add gradient from floating end point
          dc(prob.curveNo2FloatDVNo(curveNo,2),cNo) = (c(cNo)+xval(prob.pointNo2BoundDVNo(pNo)))*dpwfe(pNo); % we have to add the bound variable, to get the "real" constraint value
          dc(prob.curveNo2FloatDVNo(curveNo,2),cNo+1) = (c(cNo+1)+xval(prob.pointNo2BoundDVNo(pNo)))*dpwfe(pNo); % we have to add the bound variable, to get the "real" constraint value
        end
        % update constraint counter for the second constraint
        cNo = cNo + 1;
      end
    end
  end
  
  if options.floating
    for curveNo = 1:prob.nCurves-1
      % Extract end position of the current curve
      fDVNo = prob.curveNo2FloatDVNo(curveNo,2);
      x = xval(fDVNo);
      
      % Specify c^0 continuity: f1(x) - f2(x) = 0
      
      % Update constraint counter
      ceqNo = ceqNo + 1;
      % f1 part
      % Extract first design variable for curve
      DVNo = prob.curveNo2DVNo{curveNo}(1);
      dceq(DVNo,ceqNo) = 1;
      for orderNo = 1:prob.curveOrder(curveNo)
        DVNo = prob.curveNo2DVNo{curveNo}(1+orderNo);
        dceq(DVNo,ceqNo) = x^(orderNo);
      end
      
      % Take care of the floating DV for f1 part
      DVNo = prob.curveNo2DVNo{curveNo}(2);
      dceq(fDVNo,ceqNo) = xval(DVNo);
      for orderNo = 2:prob.curveOrder(curveNo)
        DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
        dceq(fDVNo,ceqNo) = dceq(fDVNo,ceqNo) + orderNo*xval(DVNo)*x^(orderNo-1);
      end
      
      % f2 part
      % Extract first design variable from the curve infront of current curve
      DVNo = prob.curveNo2DVNo{curveNo+1}(1);
      dceq(DVNo,ceqNo) = -1;
      for orderNo = 1:prob.curveOrder(curveNo+1)
        DVNo = prob.curveNo2DVNo{curveNo+1}(1+orderNo);
        dceq(DVNo,ceqNo) = -x^(orderNo);
      end
      
      % Take care of the floating DV for f2 part
      DVNo = prob.curveNo2DVNo{curveNo+1}(2);
      dceq(fDVNo,ceqNo) = dceq(fDVNo,ceqNo) - xval(DVNo);
      for orderNo = 2:prob.curveOrder(curveNo+1)
        DVNo = prob.curveNo2DVNo{curveNo+1}(orderNo+1);
        dceq(fDVNo,ceqNo) = dceq(fDVNo,ceqNo) - orderNo*xval(DVNo)*x^(orderNo-1);
      end
     
      % Specify continuity constraints from c^1 to c^n: d^nf1(x)/dx^n - d^nf2(x)/dx^n = 0 
      for cc = 1:prob.curveContinuity(curveNo)
        % Update constraint counter
        ceqNo = ceqNo + 1;
        % f1 part
        for orderNo = cc:prob.curveOrder(curveNo)
          DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
          dceq(DVNo,ceqNo) = prob.cuvedfCof(cc,orderNo+1)*x^(orderNo-cc);
        end
        
        % f2 part  
        for orderNo = cc:prob.curveOrder(curveNo+1)
          DVNo = prob.curveNo2DVNo{curveNo+1}(orderNo+1);
          dceq(DVNo,ceqNo) = -prob.cuvedfCof(cc,orderNo+1)*x^(orderNo-cc);
        end
        
        % Take care of the floating DV
        ccc = cc + 1;
        for orderNo = ccc+1:prob.curveOrder(curveNo)
          DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
          dceq(fDVNo,ceqNo) = dceq(fDVNo,ceqNo) + prob.cuvedfCof(ccc,orderNo+1)*xval(DVNo)*x^(orderNo-ccc);
        end
        for orderNo = ccc+1:prob.curveOrder(curveNo+1)
          DVNo = prob.curveNo2DVNo{curveNo+1}(orderNo+1);
          dceq(fDVNo,ceqNo) = dceq(fDVNo,ceqNo) - prob.cuvedfCof(ccc,orderNo+1)*xval(DVNo)*x^(orderNo-ccc);
        end
        
      end % curveContinuity
    end % nCurves
  end
  
  % Error checks
  if cNo ~= prob.nC
    err_msg = sprintf('Number of added non-linear in-equality constraints (%i) does not match expected (%i), check inputs',cNo,prob.nC);
    error(err_msg);
  end
  
  % Error checks
  if ceqNo ~= prob.nCeq
    err_msg = sprintf('Number of added non-linear in-equality constraints (%i) does not match expected (%i), check inputs',cNo,prob.nC);
    error(err_msg);
  end
  
end

function [A,b,Aeq,beq] = getLinearConstraints(prob,options)
  % Main function for assembling the linear constraints
  
  % Initialize output arrays, default to zero
  if prob.nA > 0
    A = sparse(prob.nA,prob.nDV);
    b = zeros(prob.nA,1);
  else
    A = [];
    b = [];
  end
  
  if prob.nAeq > 0
    Aeq = sparse(prob.nAeq,prob.nDV);
    beq = zeros(prob.nAeq,1);
  else
    Aeq = [];
    beq = [];
  end
  
  % Initialize constraint counter
  ANo = 0; 
  AeqNo = 0; 
  
  % Setup curve continuity constraints
  if ~options.floating
    for curveNo = 1:prob.nCurves-1
      % Extract end position of the current curve
      x = prob.curveEnd(curveNo);
      
      % Specify c^0 continuity: f1(x) - f2(x) = 0
      
      % Update constraint counter
      AeqNo = AeqNo + 1;
      % Set RHS
      beq(AeqNo) = 0;
      % f1 part
      % Extract first design variable for curve
      DVNo = prob.curveNo2DVNo{curveNo}(1);
      Aeq(AeqNo,DVNo) = 1;
      for orderNo = 1:prob.curveOrder(curveNo)
        DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
        Aeq(AeqNo,DVNo) = x^(orderNo);
      end
      
      % f2 part
      % Extract first design variable from the curve infront of current curve
      DVNo = prob.curveNo2DVNo{curveNo+1}(1);
      Aeq(AeqNo,DVNo) = -1;
      for orderNo = 1:prob.curveOrder(curveNo+1)
        DVNo = prob.curveNo2DVNo{curveNo+1}(orderNo+1);
        Aeq(AeqNo,DVNo) = -x^(orderNo);
      end
      
      % Specify continuity constraints from c^1 to c^n: d^nf1(x)/dx^n - d^nf2(x)/dx^n = 0 
      for cc = 1:prob.curveContinuity(curveNo)
        % Update constraint counter
        AeqNo = AeqNo + 1;
        % Set RHS
        beq(AeqNo) = 0;
        % f1 part
        DVNo = prob.curveNo2DVNo{curveNo}(1);
        Aeq(AeqNo,DVNo) = 0;
        for orderNo = cc:prob.curveOrder(curveNo)
          DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
          Aeq(AeqNo,DVNo) = prob.cuvedfCof(cc,orderNo+1)*x^(orderNo-cc);
        end
        % f2 part  
        DVNo = prob.curveNo2DVNo{curveNo+1}(1);
        Aeq(AeqNo,DVNo) = 0;
        for orderNo = cc:prob.curveOrder(curveNo+1);
          DVNo = prob.curveNo2DVNo{curveNo+1}(orderNo+1);
          Aeq(AeqNo,DVNo) = - prob.cuvedfCof(cc,orderNo+1)*x^(orderNo-cc);
        end
        
      end % curveContinuity
    end % nCurves
    
    % Local c0 constraints
    if ~isempty(options.c0)
      nC0 = size(options.c0,1);
      for cNo = 1:nC0
        xTarget = options.c0(cNo,1);
        yTarget = options.c0(cNo,2);
        
        % Determine which curve "covers" the specified constraint
        if xTarget <= prob.curveStart(1) % Extrapolate from first curve
          curveNo = 1;
        elseif xTarget >= prob.curveEnd(end) % Extrapolate from last curve
          curveNo = prob.nCurves;
        else % Position lies in a curve segment
          curveNos = find(xTarget >= prob.curveStart & xTarget <= prob.curveEnd);
          if ~isempty(curveNos)
            curveNo = curveNos(1);
          else
            err_msg = sprintf('Specified c0 point not located in curve intervals, this should not be possible, check inputs: c0(xTarget,yTarget) = (%d,%d)',xTarget,yTarget);
            error(err_msg);
          end
        end
        % Update constraint counter
        AeqNo = AeqNo + 1;
        % Set rhs
        beq(AeqNo) = yTarget;
        % extract coefficients
        % Set first value
        DVNo = prob.curveNo2DVNo{curveNo}(1);
        Aeq(AeqNo,DVNo) = 1;
        % Set the remaining coefficients
        for orderNo = 1:prob.curveOrder(curveNo)
          DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
          Aeq(AeqNo,DVNo) = xTarget^orderNo;
        end
      end
    end
    
    % Local c1 constraints
    if ~isempty(options.c1)
      nC1 = size(options.c1,1);
      for cNo = 1:nC1
        xTarget = options.c1(cNo,1);
        yTarget = options.c1(cNo,2);
        
        % Determine which curve "covers" the specified constraint
        if xTarget <= prob.curveStart(1) % Extrapolate from first curve
          curveNo = 1;
        elseif xTarget >= prob.curveEnd(end) % Extrapolate from last curve
          curveNo = prob.nCurves;
        else % Position lies in a curve segment
          curveNos = find(xTarget >= prob.curveStart & xTarget <= prob.curveEnd);
          if ~isempty(curveNos)
            curveNo = curveNos(1);
          else
            err_msg = sprintf('Specified c1 point not located in curve intervals, this should not be possible, check inputs: c0(xTarget,yTarget) = (%d,%d)',xTarget,yTarget);
            error(err_msg);
          end
        end
        % Update constraint counter
        AeqNo = AeqNo + 1;
        % Set rhs
        beq(AeqNo) = yTarget;
        % extract coefficients
        % Set first value
        DVNo = prob.curveNo2DVNo{curveNo}(1);
        Aeq(AeqNo,DVNo) = 0;
        % Set the remaining coefficients
        for orderNo = 1:prob.curveOrder(curveNo)
          DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
          Aeq(AeqNo,DVNo) = orderNo*xTarget^(orderNo-1);
        end
      end
    end
    
    % Lower and upper bounds on points
    if ~isempty(options.pointlb) || ~isempty(options.pointub)
      for pNo = 1:options.nP
        xTarget = options.points(pNo,2);
        curveNo = prob.pointNo2CurveNo(pNo);
        
        % Lower bound
        if ~isempty(options.pointlb) 
          % Update constraint counter
          ANo = ANo + 1;
          % Set rhs for lb
          b(ANo) = -options.pointlb(pNo);
          % extract coefficients
          % Set first value
          DVNo = prob.curveNo2DVNo{curveNo}(1);
          A(ANo,DVNo) = -1;
          % Set the remaining coefficients
          for orderNo = 1:prob.curveOrder(curveNo)
            DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
            A(ANo,DVNo) = -xTarget^(orderNo);
          end
        end
        
        % Upper bound
        if ~isempty(options.pointub) 
          % Update constraint counter
          ANo = ANo + 1;
          % Set rhs for lb
          b(ANo) = options.pointub(pNo);
          % extract coefficients
          % Set first value
          DVNo = prob.curveNo2DVNo{curveNo}(1);
          A(ANo,DVNo) = 1;
          % Set the remaining coefficients
          for orderNo = 1:prob.curveOrder(curveNo)
            DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
            A(ANo,DVNo) = xTarget^(orderNo);
          end
        end
      end % for np
    end % pointlb and pointub
    
  end % not floating
  
  if options.floating
    % Setup linear in-equality constraints to ensure that start position variables increase from first to last i.e., f1<=f2 --> f1-f2<=0
    for fNo = 1:prob.nFloatDV-1
      % Extract end position of the current curve
      fDVNo1 = prob.floatNo2DVNo(fNo);
      fDVNo2 = prob.floatNo2DVNo(fNo+1);
      ANo = ANo + 1;
      A(ANo,fDVNo1) = 1;
      A(ANo,fDVNo2) = -1;
    end
    
  end
  
  
  % Error checks
  if ANo ~= prob.nA
    err_msg = sprintf('Number of added linear in-equality constraints (%i) does not match expected (%i), check inputs',ANo,prob.nA);
    error(err_msg);
  end
  
  if AeqNo ~= prob.nAeq
    err_msg = sprintf('Number of added linear equality constraints (%i) does not match expected (%i), check inputs',AeqNo,prob.nAeq);
    error(err_msg);
  end
  
end

function [fval,df] = evalFit(xIn,settings,xval,varargin)
  % This is the main part of the function used to evaluate the fit at various locations (xIn). 
  % The function is passed to the user as output for easy evaluation of the fit.
  
  p = inputParser;
  p.CaseSensitive = false;
  % Helper functions for input parser
  checkEmptyOrNumericPositive = @(x)(isempty(x) || ((isnumeric(x) || isscalar(x)) && all(x > 0)));
  checkPoints =@(x) isnumeric(x) && min(size(x))==1;
  % The order of the required inputs matter
  p.addRequired('x',checkPoints);
  % 
  p.addParameter('plot',false,@(x)islogical(x));
  p.addParameter('df',[], @(x)checkEmptyOrNumericPositive(x));
  
    % pars input
  if nargin < 4 || isempty(varargin)
      parse(p,xIn);
  else
      temp = varargin{:};
      parse(p,xIn,temp{:});
  end
  % Output results to options structure
  options = p.Results;
  
  [x,idx] = sort(xIn);
  nP = numel(x);
  fval = zeros(size(x)); % same output size as input
  
  if ~isempty(options.df)
    ndf = numel(options.df);
    df = zeros(nP,ndf);
  else
    df = []; % Initialize to empty
    ndf = 0;
  end
  
  for pNo = 1:nP
    
    % Determine which curve "covers" the specified constraint
    if x(pNo) <= settings.curveStart(1) % Extrapolate from first curve
      curveNo = 1;
    elseif x(pNo) >= settings.curveEnd(end) % Extrapolate from last curve
      curveNo = settings.nCurves;
    else % Position lies in a curve segment
      curveNos = find(x(pNo) >= settings.curveStart & x(pNo) <= settings.curveEnd);
      if ~isempty(curveNos)
        curveNo = curveNos(1); % if the point lies at an "intersection" between two curves, pick the first curve
      else
        err_msg = sprintf('Specified point not located in curve intervals, this should not be possible, check inputs: x(%i) = %d',x(idx(pNo)),x(pNo));
        error(err_msg);
      end
    end
    
    DVNo = settings.curveNo2DVNo{curveNo}(1);
    fval(pNo) = xval(DVNo);
    for orderNo = 1:settings.curveOrder(curveNo)
      DVNo = settings.curveNo2DVNo{curveNo}(1+orderNo);
      fval(pNo) = fval(pNo) + xval(DVNo)*x(pNo)^orderNo;
    end
    
    if ~isempty(options.df)
      for dfNo = 1:ndf
        cc = options.df(dfNo);
        for orderNo = cc:settings.curveOrder(curveNo)
          DVNo = settings.curveNo2DVNo{curveNo}(orderNo+1);
          df(pNo,dfNo) = df(pNo,dfNo) + settings.cuvedfCof(cc,orderNo+1)*xval(DVNo)*x(pNo)^(orderNo-cc);
        end
      end
    end
  end
  
  if options.plot
    nPlots = 1 + ndf;
    % Get some nice colors from linspecer, see  Jonathan C. Lansey (2020). Beautiful and distinguishable line colors + colormap 
    % (https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap)
    %  MATLAB Central File Exchange. Retrieved February 18, 2020. 
    colors = linspecer(nPlots); %
    figure('Name','curvefit plot');
    subplot(nPlots,1,1)
    plot(x,fval,'color',colors(1,:));
    xlabel('x')
    ylabel('f(x)')
    grid on
    % Plot all the derivatives, if any
    for plotNo = 2:nPlots
      subplot(nPlots,1,plotNo)
      plot(x,df(:,plotNo-1),'color',colors(plotNo,:));
      xlabel('x')
      if options.df(plotNo-1) > 1
        ylabel(sprintf('d^%if/dx^%i',plotNo-1,plotNo-1))
      else
        ylabel('df/dx')
      end
      grid on
    end
  end
  
  
end

function plotFit(xval,settings,points)
  % This function is used to generate a plot which shows the fit and the "target" points.
  % A comparison between the fit and the target points is also shown.
  
  x = unique([points(:,2);linspace(min(points(:,2)),max(points(:,2)),5*numel(points(:,2)))']);
  nP = numel(x);
  fval = zeros(nP,1);
  for pNo = 1:nP
    % Determine which curve "covers" the specified constraint
    if x(pNo) <= settings.curveStart(1) % Extrapolate from first curve
      curveNo = 1;
    elseif x(pNo) >= settings.curveEnd(end) % Extrapolate from last curve
      curveNo = settings.nCurves;
    else % Position lies in a curve segment
      curveNos = find(x(pNo) >= settings.curveStart & x(pNo) <= settings.curveEnd);
      if ~isempty(curveNos)
        curveNo = curveNos(1); % if the point lies at an "intersection" between two curves, pick the first curve
      else
        err_msg = sprintf('Specified point not located in curve intervals, this should not be possible, check inputs: x(%i) = %d',x(idx(pNo)),x(pNo));
        error(err_msg);
      end
    end
    
    DVNo = settings.curveNo2DVNo{curveNo}(1);
    fval(pNo) = xval(DVNo);
    for orderNo = 1:settings.curveOrder(curveNo)
      DVNo = settings.curveNo2DVNo{curveNo}(1+orderNo);
      fval(pNo) = fval(pNo) + xval(DVNo)*x(pNo)^orderNo;
    end
  end
  
  % get max and min values from points and fval
  minVal = min([fval(:);points(:,1)]);
  maxVal = max([fval(:);points(:,1)]);
  nCurvePos = (settings.nCurves + 1);
  curvePos = zeros(nCurvePos,1);
  curvePos(1:end-1) = settings.curveStart(:);
  curvePos(end) = settings.curveEnd(end);
  colors = linspecer(3); %
  figure('Name','curvefit: Compare fit to target');
  subplot(2,1,1)
  hold on
  plot(x,fval,'color',colors(1,:));
  scatter(points(:,2),points(:,1),'o');
  for ii = 1:nCurvePos
    plot([curvePos(ii),curvePos(ii)],[minVal,maxVal],'color',colors(3,:));
  end
  legend('Fit','Target points','Curve Start/End')
  xlabel('x')
  ylabel('f(x)')
  grid on
  subplot(2,1,2) % show difference in value between fit and input points
  % We can interpolate to get exact values as we know the positions have been evaluated
  fp = interp1(x,fval,points(:,2));
  plot(points(:,2),fp-points(:,1),'color',colors(2,:));
  xlabel('x')
  ylabel('f(x)-target')
  legend('Difference')
  grid on
end

function pwf = getFixedFloatingStartPointWeightsForCurveInterval(pointPos,endpos,beta) 
  pwf = (1-HS(endpos,pointPos,beta));
end

function dpwf = getFixedFloatingStartPointWeightsForCurveIntervalDSA(pointPos,endpos,beta) 
  dpwf = 1-HSDSA(endpos,pointPos,beta);
end

function pwf = getFixedFloatingEndPointWeightsForCurveInterval(pointPos,startpos,beta) 
  pwf = HS(startpos,pointPos,beta);
end

function dpwf = getFixedFloatingEndPointWeightsForCurveIntervalDSA(pointPos,startpos,beta) 
  dpwf = HSDSA(startpos,pointPos,beta);
end

function [pwf] = getPointWeightsForCurveInterval(pointPos,startpos,endpos,beta) 
  pwf = HS(startpos,pointPos,beta).*(1-HS(endpos,pointPos,beta));
end

function [dpwfs] = getFloatingStartPointDSA(pointPos,startpos,endpos,beta) 
  % Derivative of getPointWeightsForCurveInterval wrt., startpos
  [dwfs] = HSDSA(startpos,pointPos,beta);
  [wfe] = HS(endpos,pointPos,beta);
  dpwfs = dwfs.*(1-wfe);
end

function [dpwfe] = getFloatingEndPointDSA(pointPos,startpos,endpos,beta) 
  % Derivative of getPointWeightsForCurveInterval wrt., endpos
  [dwfe] = HSDSA(endpos,pointPos,beta);
  [wfs] = HS(startpos,pointPos,beta);
  dpwfe = -wfs.*dwfe;
end

function [wf] = HS(sx,ps,beta)
  % Unit step function / projection filter
  % sx is the normalized position of a point which we try to fit the curve to. 
  % ps is the normalized design variable associated with either a start or end position of a curve.
  % beta is a slope parameter for this particular unit step approximation. 
  % The formulation is based on Wang,  F.,  Lazarov,  B.,  and  Sigmund,  O 2011 and Soerensen, R, and Lund, E 2015
    wf = (tanh(beta*sx)+tanh(beta*(ps-sx)))./(tanh(beta*sx)+tanh(beta*(1-sx)));
    wf(isnan(wf))=0;
 end
 
 function [dwf] = HSDSA(sx,ps,beta)
   % Derivative of HS wrt., the position ps
    dwf =(((1-tanh(beta*sx).^2)-(1-tanh(beta*(ps-sx)).^2)*beta)/(tanh(beta*sx)+tanh(beta*(1-sx))) ...
    - (tanh(beta*ps)+tanh(beta*(sx-ps)))*beta*((1-tanh(beta*sx).^2)-(1-tanh(beta*(1-sx)).^2))./(tanh(beta*ps)+tanh(beta*(1-sx))).^2);
  
  dwf(isnan(dwf))=0;
 end
