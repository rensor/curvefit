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

    % check points structure
    if  size(points,2) ~= 2
    exitflag = 0;
    message = sprintf('curvefit: input points needs to be a nP x 2 array. Number of columns for input is %3i', size(points,2));
    fhandle = [];
    return
    end

    options = setOptions(points,nCurves,curveOrder,varargin);
    prob = problemformulation(options);
    x = linspace(0,1,1000);
    [pwf] = getPointWeightsForCurveInterval(x,0.25,0.75,100);
    [dpwfs] = getFloatingStartPointDSA(x,0.25,0.75,100);
    figure;
    subplot(2,1,1)
    plot(x,pwf)
    subplot(2,1,2)
    plot(x,dpwfs)
end

function [options]=setOptions(points,nCurves,curveOrder,input)
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
    p.addParameter('curveContinuity',2,  @(x)checkScalarNumPos(x));

    % Settings for floating curve intersection points
    p.addParameter('floating',false,  @(x)islogical(x));
    p.addParameter('beta',100,  @(x)checkScalarNumPos(x));

    % General settings
    p.addParameter('display','off',  @(x)checkEmpetyOrChar(x));
    p.addParameter('plot',false,  @(x)islogical(x));

    % Additional constraints
    p.addParameter('c0',[],  @(x)checkEmptyOrScalarNum(x));
    p.addParameter('c1',[],  @(x)checkEmptyOrScalarNum(x));
    p.addParameter('pointlb',[],  @(x)checkEmptyOrScalarNum(x));
    p.addParameter('pointub',[],  @(x)checkEmptyOrScalarNum(x));
    p.addParameter('kappa',[],  @(x)checkEmptyOrNumericPositive(x));

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

end

function [prob,exitflag,message] = problemformulation(options)
    % Load options to prob structure
    prob = options;
    prob.nP = size(prob.points,1);

    % Determine how the user has specified the problem
    if numel(options.startpos) == 1
        % distribute start positions for the curves between the specified start and end position
        if isempty(options.endpos)
            options.endpos = max(prob.points(:,2));
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

    if ~prob.floating
    else
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
    dwf =((1-tanh(beta*xs).^2)-(1-tanh(beta*(ps-xs)).^2)*beta)/(tanh(beta*xs)+tanh(beta*(1-xs))) ...
    - (tanh(beta*ps)+tanh(beta*(xs-ps)))*beta*((1-tanh(beta*xs).^2)-(1-tanh(beta*(1-xs)).^2))./(tanh(beta*ps)+tanh(beta*(1-xs))).^2;
 end