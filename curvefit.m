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

  [options,exitflag,message_out] = setOptions(points,nCurves,curveOrder,varargin);
  message = sprintf('%s \n',message_out);
  [prob,exitflag,message_out]  = problemformulation(options);
  message = sprintf('%s \n',message_out);
  
  [A,b,Aeq,beq] = getLinearConstraints(prob,options);
  
  fun = @(x) getObjectiveFunction(x,prob);
  nonlcon = @(x) getNonLinearConstraints(x,prob,options);
  
  opti = fminslp(fun,prob.x0,A,b,Aeq,beq,prob.xL,prob.xU,nonlcon,'display','iter','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',false,'solver','glpk');
  [xval,fval,exitflag,output] = opti.solve;
  message = output.message;
  
  % Make minimal structure for function evaluation
  settings = struct('curveStart',[],...
                    'curveEnd',[],...
                    'nCurves',[],...
                    'curveNo2DVNo',[],...
                    'curveOrder',[]);
  settings.curveStart = prob.curveStart;
  settings.curveEnd = prob.curveEnd;
  settings.nCurves = prob.nCurves;
  settings.curveNo2DVNo = prob.curveNo2DVNo;
  settings.curveOrder = prob.curveOrder;
  fhandle = @(x,varargin) evalFit(x,settings,xval,varargin);
  
  if options.plot
    plotFit(xval,settings,points);
  end
  
  
end

function [options,exitflag,message]=setOptions(points,nCurves,curveOrder,input)
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
  
  if ~options.floating
    [prob,exitflag,message] = setupStandardFormulation(prob,options);
  else
    
  end % not floating formulation
  
  % Total number of design variables
  prob.nDV = prob.nCurveDV + prob.nBoundDV + prob.nFloatDV;
 
  % Total number of constraints
  prob.nG = prob.nA + prob.nAeq + prob.nC + prob.nCeq;
  
  % Allocate design variable data structures
  prob.xL = repmat(-10000,[prob.nDV,1]); % Initial lower bound for design variables
  prob.xU = repmat(10000,[prob.nDV,1]); % Initial upper bound for design variables
  if strcmpi(options.method,'bound')
    % Modify upper and lower bounds for bound design variables
    prob.xL(prob.nCurveDV+1:prob.nCurveDV+prob.nBoundDV) = 0;
    prob.xU(prob.nCurveDV+1:prob.nCurveDV+prob.nBoundDV) = 1e6;
  end
  % Initialize design variables
  prob.x0 = zeros(prob.nDV,1);
  
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
    prob.nBoundDV = prob.nCurves;
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

function [fval,df] = getObjectiveFunction(xval,prob)
  switch prob.method
    case 'bound'
      fval = sum(xval(prob.nCurveDV+1:prob.nCurveDV+prob.nBoundDV));
      if nargout > 1
        df = zeros(prob.nDV,1);
        df(prob.nCurveDV+1:prob.nCurveDV+prob.nBoundDV) = 1;
      end
  end
end


function [c,ceq,dc,dceq] = getNonLinearConstraints(xval,prob,options)
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
        c(cNo) = (f1-y) - xval(prob.nCurveDV+curveNo);
        
        % update constraint counter
        cNo = cNo + 1;
        % target - f(x) - bound <= 0
        c(cNo) = (y-f1) - xval(prob.nCurveDV+curveNo);
      end
    end
  end
  
  % Error checks
  if cNo ~= prob.nC
    err_msg = sprintf('Total number of non-linear in-equality constraints(%i) does not match specified(%i)',cNo,prob.nC)
    error(err_msg);
  end
  
  if ceqNo ~= prob.nCeq
    err_msg = sprintf('Total number of non-linear equality constraints(%i) does not match specified(%i)',ceqNo,prob.nCeq)
    error(err_msg);
  end
  
  % DSA
  if nargout > 2
    [dc,dceq] = getNonLinearConstraintsDSA(xval,prob,options);
  end
  
  
end

function [dc,dceq] = getNonLinearConstraintsDSA(xval,prob,options)
  
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
        dc(prob.nCurveDV+curveNo,cNo) = -1;
        dc(prob.nCurveDV+curveNo,cNo+1) = -1;
        % update constraint counter for the second constraint
        cNo = cNo + 1;
      end
    end
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
      % f1 part
      % Extract first design variable for curve
      DVNo = prob.curveNo2DVNo{curveNo}(1);
      Aeq(AeqNo,DVNo) = 1;
      for orderNo = 1:prob.curveOrder(curveNo)
        DVNo = prob.curveNo2DVNo{curveNo}(1+orderNo);
        Aeq(AeqNo,DVNo) = x^(orderNo);
      end
      
      % f2 part
      % Extract first design variable from the curve infront of current curve
      DVNo = prob.curveNo2DVNo{curveNo+1}(1);
      Aeq(AeqNo,DVNo) = -1;
      for orderNo = 1:prob.curveOrder(curveNo+1)
        DVNo = prob.curveNo2DVNo{curveNo+1}(1+orderNo);
        Aeq(AeqNo,DVNo) = -x^(orderNo);
      end
      
      % Specify continuity constraints from c^1 to c^n: d^nf1(x)/dx^n - d^nf2(x)/dx^n = 0 
      for cc = 1:prob.curveContinuity(curveNo)
        % Update constraint counter
        AeqNo = AeqNo + 1;
        % f1 part
        for orderNo = cc:prob.curveOrder(curveNo)
          % calculate orderNo*(orderNo-1)*(orderNo-2)*...*(orderNo-n)
          devPart = orderNo;
          for ii = 1:cc-1
            devPart = devPart*(orderNo-ii);
          end
          DVNo = prob.curveNo2DVNo{curveNo}(orderNo+1);
          Aeq(AeqNo,DVNo) = devPart*x^(cc-1);
        end
        % f2 part  
        for orderNo = cc:prob.curveOrder(curveNo+1);
          % calculate orderNo*(orderNo-1)*(orderNo-2)*...*(orderNo-n)
          devPart = orderNo;
          for ii = 1:cc-1
            devPart = devPart*(orderNo-ii);
          end
          DVNo = prob.curveNo2DVNo{curveNo+1}(orderNo+1);
          Aeq(AeqNo,DVNo) = -devPart*x^(cc-1);
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
      for cNo = 1:nC0
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
      for pNo = 1:options.nP;
        xTarget = options.points(pNo,1);
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
            A(ANo,DVNo) = -orderNo*xTarget^(orderNo-1);
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
            A(ANo,DVNo) = orderNo*xTarget^(orderNo-1);
          end
        end
        
      end % for np
    end % pointlb and pointub
  end % not floating
  
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
    % Here you can add new options if needed
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
          % calculate orderNo*(orderNo-1)*(orderNo-2)*...*(orderNo-n)
          devPart = orderNo;
          for ii = 1:cc-1
            devPart = devPart*(orderNo-ii);
          end
          DVNo = settings.curveNo2DVNo{curveNo}(orderNo+1);
          df(pNo,dfNo) = df(pNo,dfNo) + devPart*xval(DVNo)*x(pNo)^(orderNo-cc);
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
  colors = linspecer(2); %
  figure('Name','curvefit: Compare fit to target');
  subplot(2,1,1)
  hold on
  plot(x,fval,'color',colors(1,:));
  scatter(points(:,2),points(:,1),'o');
  legend('Fit','target points')
  xlabel('x')
  ylabel('f(x)')
  grid on
  subplot(2,1,2)
  % We can interpolate to get exact values as we know the positions have been evaluated
  fp = interp1(x,fval,points(:,2));
  plot(points(:,2),fp-points(:,1),'color',colors(2,:));
  xlabel('x')
  ylabel('f(x)-target')
  legend('Difference')
  grid on
end


function [dpwfs] = getFloatingStartPointDSA(pointPos,startpos,endpos,beta) 
  % Derivative of getPointWeightsForCurveInterval wrt., startpos
  [dwfs] = HSDSA(startpos,pointPos,beta);
  [wfe] = HS(endpos,pointPos,beta);
  dpwfs = dwfs.*(1-wfe);
end

function [pwf] = getPointWeightsForCurveInterval(pointPos,startpos,endpos,beta) 
  pwf = HS(startpos,pointPos,beta).*(1-HS(endpos,pointPos,beta));
end

function [wf] = HS(xs,ps,beta)
  % Unit step function / projection filter
  % xs is the normalized position of a point which we try to fit the curve to. 
  % ps is the normalized design variable associated with either a start or end position of a curve.
  % beta is a slope parameter for this particular unit step approximation. 
  % The formulation is based on Wang,  F.,  Lazarov,  B.,  and  Sigmund,  O 2011 and Soerensen, R, and Lund, E 2015
    wf = (tanh(beta*xs)+tanh(beta*(ps-xs)))./(tanh(beta*xs)+tanh(beta*(1-xs)));
 end
 
 function [dwf] = HSDSA(xs,ps,beta)
   % Derivative of HS wrt., the position ps
    dwf =(((1-tanh(beta*xs).^2)-(1-tanh(beta*(ps-xs)).^2)*beta)/(tanh(beta*xs)+tanh(beta*(1-xs))) ...
    - (tanh(beta*ps)+tanh(beta*(xs-ps)))*beta*((1-tanh(beta*xs).^2)-(1-tanh(beta*(1-xs)).^2))./(tanh(beta*ps)+tanh(beta*(1-xs))).^2);
 end