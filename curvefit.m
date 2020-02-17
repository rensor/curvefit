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
  
  x = linspace(0,1,1000);
  [pwf] = getPointWeightsForCurveInterval(x,0.25,0.75,100);
  [dpwfs] = getFloatingStartPointDSA(x,0.25,0.75,100);
  figure;
  subplot(2,1,1)
  plot(x,pwf)
  subplot(2,1,2)
  plot(x,dpwfs)
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
  p.addParameter('curveContinuity',curveOrder, @(x)checkScalarNumPos(x));

  % Settings for floating curve intersection points
  p.addParameter('floating',false,  @(x)islogical(x));
  p.addParameter('beta',100,  @(x)checkScalarNumPos(x));

  % General settings
  p.addParameter('display','off',  @(x)checkEmpetyOrChar(x));
  p.addParameter('plot',false,  @(x)islogical(x));

  % Additional constraints
  p.addParameter('c0',[],  @(x)checkEmptyOrScalarNum(x)); % local c0
  p.addParameter('c1',[],  @(x)checkEmptyOrScalarNum(x)); % local c1
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
  
    % Determine if curve continuity has been specified for each curve or just one value which should apply to all curves
  if numel(options.curveContinuity) ~= options.nCurves
    if numel(options.curveContinuity) == 1
      options.curveContinuity = repmat(options.curveContinuity,[options.nCurves,1]);
    else
       message = sprintf('Please only specify curve continuity either by one value, or for each curve');
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
  prob.nA       = 0; % Number of linear in-equality constraints
  prob.nAeq     = 0; % Number of linear equality constraints
  prob.nC       = 0; % Number of non-linear in-equality constraints
  prob.nCeq     = 0; % Number of non-linear equality constraints
  prob.nG       = 0; % Total number of constraints
  
  if ~options.floating
    [prob,exitflag,message] = setupStandardFormulation(prob,options);
  else
    
  end % not floating formulation
 
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
  prob.nAeq = prob.nAeq + sum(options.curveContinuity(1:end-1)+1);
  
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


function [c,ceq,dc,dceq] = getNonLinearConstraints(xval,prob)

end

function [fval,df] = getObjectiveFunction(xval,prob)

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