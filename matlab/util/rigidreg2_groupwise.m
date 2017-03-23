function [Iregd,rigidparams,varargout] = rigidreg2_groupwise(Itemplate,Itarget,N, varargin)
% AFFINEREG2 Register a pair of images using a 2D affine transformation.
%
%   ITARG_REGD = AFFINEREG2(ITEMPLATE,ITARGET,N) registers ITARGET (moving)
%   to ITEMPLATE (stationary) by optimizing an affine transformation, and 
%   returns the registered image, ITARGET subject to the determined 
%   transformation. The number of graylevels *must* be specified in N.
%
%   [ITARG_REGD,RIGIDPARAMS] = AFFINEREG2(...) returns the determined
%   transformation parameters RIGIDPARAMS, a 5-by-1 array.  See rotmat.
%
%   [ITARG_REGD,RIGIDPARAMS,SIMVALS] = AFFINEREG2(...) returns the value of
%   the similarity measure at each iteration.
%
%   [ITARG_REGD,RIGIDPARAMS,SIMVALS,M] = AFFINEREG2(...) returns the linear
%   homogeneous coorinate transformation matrix M of size 4-by-4.
%
%   [ITARG_REGD,RIGIDPARAMS,SIMVALS,M,FINALMI] = AFFINEREG2(...) returns
%   the MI between ITARG_REGD and ITEMPLATE.
%
%   [ ... ] = AFFINEREG2(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control the registration.
%   Parameters are:
%
%   'InterpMethod' - The coordinate system interpolation method to use when
%      calculating intensities on the deformed target image. Possible 
%      values are 'nearest' (default), 'linear', and 'cubic'. Both 'linear'
%      and 'cubic' can return values not present in I, and 'cubic' can 
%      return values outside the range of values in I.
%
%   'Measure' - The similariy measure to drive transformation optimization.
%      Possible values are 'MI' (default), 'NMI', 'ECC', 'JE', and 'EOH'.
%
%   'BGValue' - The background intensity to use when pixels from outside of
%      the frame of the original image are required. An integer, default 0.
%
%   'NoDisp' - Whether or not to display information and images at each
%      transformation iteration.  Boolean value, default true.
%
%   'StartParams' - Starting transformation parameters (see imdef). Double 
%      array of size 5-by-1, representing [rz; tx; ty; sx; sy]. Default is
%      [0; 0; 0; 1; 1;]. Optionaly, scale can be omitted by using a 3-by-1 
%      array. A 4-by-1 array provides 1 scale parameter for isometric xy 
%      scale (i.e. [rz tx ty sxy]').
%
%   'ParamWeights' - Weight given to each parameter during optimization.
%      Double array of size 5-by-1, optionally length 3 or 4, with the same
%      meanings as in StartParams. Default is [125; 4000; 4000; .25; .25].
%      For example, if you thing there is a large ammount of translation
%      that needs to be corrected, increase the second and third elements. 
%
%   'ParamOffsets' - Constant values added to the parameter set at every
%      iteration. This is functionally different from 'StartParams' in that
%      these parameters never change, and the optimizer is blinded to them.
%
%   'scales' - Constrain scale between xyz dimensions. 
%         'scales' Value     |            Effect            
%         ---------------------------------------------------
%                0           |           no scaling         
%                1           |       isotropic (single xy scale)
%                2 (default) |  anisotropic (independent x,y scales)
%      IMPORTANT: StartParams and ParamWeights lengths override 'scales'.
%
%   'initM' - Initial transformation matrix (see imdef).  Must be 4-by-4.
%      Applied before any other parameters, either from 'StartParams' or
%      'ParamOffsets'.
%
%   'coordshift' - Offset from volume center, about which deformations are 
%      performed, [Y; X]. Positive values impose a shift toward origin.
%      Default coordshift is [0; 0]. MATLAB coordinates this time! This
%      parameter is not properly tested yet.
%
%   NOTE: This function enforces N as a required input for two reasons,
%   (1) the operator will be less likely to make a mistake about the range
%   of the data, and (2) determining N automatically would add at least 
%   some computation time.
%
%   Intensities of input images must must fall within 0 to N-1!! Scale your
%   images appropriately, or registration will fail or be unpredictable.
%
%   See also: imdef, rotmat, objdef2d, affinereg3, imdef3, rotmat3.
%
%JC

% Parse inputs
ip=inputParser;
% ip.addRequired('Itarget', @isnumeric);
% ip.addRequired('Itemplate', @isnumeric);
% ip.addRequired('N', @isnumeric);
%ip.addOptional('interpmethod', 'nn', @(x)any(strcmpi(x,{'nearest','linear','blurjit'})));
interpmethods={'nearest','linear','cubic','*nearest','*linear','*cubic'};
ip.addParamValue('interpmethod', 'nearest', @(x)any(strcmpi(x,interpmethods)));
measures={'MI','NMI','ECC','JE','EOH','MMI'};
ip.addParamValue('measure', 'MI', @(x)any(strcmpi(x,measures)));
ip.addParamValue('bgvalue', 0, @(x)isnumeric(x) && numel(x)==1);
ip.addParamValue('nodisp', true, @islogical);
nmaxparams=6; nminparams=3; % acceptable number of parameters in startparams and paramweights
ip.addParamValue('startparams', [0; 0; 0; 1; 1; 0;], @(x)isfloat(x) && numel(x)>=nminparams && numel(x)<=nmaxparams);
ip.addParamValue('paramweights', [125; 4000; 4000; .01; .01; 125], @(x)isfloat(x) && numel(x)>=nminparams && numel(x)<=nmaxparams);
ip.addParamValue('scales', 2, @(x)isnumeric(x) && numel(x)==1 && any(x==[0 1 2]));
ip.addParamValue('paramoffsets', zeros(6,1), @(x)isfloat(x) && numel(x)>=nminparams && numel(x)<=nmaxparams);
ip.addParamValue('initM', eye(4), @(x)isfloat(x) && numel(x)==16 && size(x,1)==4);
ip.addParamValue('coordshift', zeros(2,1), @(x)isfloat(x) && numel(x)==2);
ip.addParamValue('Itarget_bg', zeros(size(Itarget,1), size(Itarget,2)), @(x)ismatrix(x) && numel(size(Itarget))==2); 
ip.addParamValue('number_many_slices', 8, @(x)isnumeric(x) && numel(x)==1); 

% ip.parse(Itemplate,Itarget,N,varargin{:}); % get Itemplate,Itarget,N directly
ip.parse(varargin{:});

% Check input 'startparams' and 'paramweights', compare lengths
% defaultfields=intersect(ip.UsingDefaults,{'startparams','paramweights'});
customfields=setdiff({'startparams','paramweights'},ip.UsingDefaults);
ninputparams=numel(customfields);
if ninputparams,
    inputparamlengths=zeros(ninputparams,1);
    for k=1:ninputparams,
        customfield=char(customfields(k));
        inputparamlengths(k)=numel(ip.Results.(customfield));
    end
    if any(diff(inputparamlengths)),
        error('affinereg2:StartParamsParamWeightsLengthMismatch','Input StartParams and ParamWeights should be the same length.'); end
    noptimparams=inputparamlengths(1);
    
    % Based on the number of input parameters, decide how to scale
    gotR2=noptimparams==6; % first guess
    scales=noptimparams-nminparams-double(gotR2); % first guess

    % Verify that input scales agree with any input params
    if ~any(strcmp(ip.UsingDefaults,{'scales'})), % scales was input
        if ~gotR2 && scales==2 && ip.Results.scales==1, % 5 parameters, scales==1, last should be theta2
            scales=ip.Results.scales; % --> override scales guess
            gotR2=true; % --> override gotR2 guess
        elseif scales~=ip.Results.scales, % scales was input, and does not agree with other params
            scales
            ip.Results.scales
            error('affinereg2:InputScalesAndParamsMismatch','Length of input(s) StartParams/ParamWeights do not agree with scales. Omit scales or fix StartParams/ParamWeights.');
        end
        %warning('affinereg2:IgnoringInputScales','Input scales redundant with StartParams and/or ParamWeights.\n');
    end
else % no input ParamWeights or StartParams, only defaults
    scales=ip.Results.scales;
    gotR2=true;
    noptimparams=nminparams+scales+double(gotR2);
end

% Trim default 'startparams' and 'paramweights' to input length, inputs are
% already verified above to have the correct length
if numel(ip.Results.startparams)>noptimparams,
startparams=ip.Results.startparams([1:nminparams+scales 6*ones(1,gotR2)]);
else startparams=ip.Results.startparams(:); end
if numel(ip.Results.paramweights)>noptimparams,
paramweights=ip.Results.paramweights([1:nminparams+scales 6*ones(1,gotR2)]);
else paramweights=ip.Results.paramweights(:); end
startparams=startparams(:); paramweights=paramweights(:); % columns
paramweights(paramweights==0)=eps*8;

% Pad paramoffsets to nmaxparams, if necessary
paramoffsets=ip.Results.paramoffsets(:);
nparamoffsets=numel(paramoffsets);
% deal with 1 or more scale input when scales==1
if scales==1 && nparamoffsets==4, % assume only one scale given, no theta2
    paramoffsets(end+1)=paramoffsets(end); % extends paramoffsets by duplicating scale
end
if nparamoffsets==5 && scales==1, warning('Assuming two scales and no theta2 in paramoffsets!'); end
paramoffsets=[paramoffsets; zeros(nmaxparams-nparamoffsets,1)]; % extends paramoffsets by zero padding

% Output information
paramstring=repmat('%.4g, ',1,noptimparams);
fprintf(['StartParams=[' paramstring '\b\b].\n'],startparams);
fprintf(['ParamWeights=[' paramstring '\b\b].\n'],paramweights);
paramstring=repmat('%.4g, ',1,nmaxparams);
fprintf(['ParamOffsets=[' paramstring '\b\b].\n'],paramoffsets);

% Make some short variable names
interpmethod=ip.Results.interpmethod;
bgvalue=cast(ip.Results.bgvalue,class(Itarget));
preM=ip.Results.initM;

Itarget_bg = ip.Results.Itarget_bg;
number_many_slices = ip.Results.number_many_slices;

%%% Register
if ~ip.Results.nodisp,
    subplot(231)
    imdisp(Itemplate(:,:,number_many_slices+1))
    subplot(232)
    imdisp(Itemplate(:,:,number_many_slices))
    subplot(233)
    imdisp(Itarget_bg)
    %imdisp(Itemplate(:,:,7))
    
    subplot(234)
    imdisp(sum(Itemplate,3));
    subplot(235)
    imdisp(Itemplate(:,:,number_many_slices))
    
    title('Template slice')

    subplot(236)
    if ~isempty(intersect(ip.UsingDefaults,{'initM'})),
        Iinit=imdef(startparams,Itarget(:,:,1),interpmethod,ip.Results.bgvalue,...
            preM,ip.Results.coordshift);
    else
        Iinit=Itarget(:,:,1);
    end
    imdisp(Iinit+double(Itarget_bg))
    title('Target slice')
    
    optoptions = optimset('Display','iter','MaxIter',100);
else
    optoptions = optimset('Display','notify','MaxIter',100);
end

% Adjust input startparams by paramweights to give proper deformation inside objective function
startparams_objdef=startparams(:);
startparams_objdef(1:nminparams)=startparams_objdef(1:nminparams)./paramweights(1:nminparams); % rot1 and translations
if noptimparams>nminparams, % got some scales
    startparams_objdef(nminparams+1:noptimparams-gotR2) = ...
        (startparams_objdef(nminparams+1:noptimparams-gotR2)-1)./paramweights(nminparams+1:noptimparams-gotR2) + 1; % only scales
end
if gotR2, startparams_objdef(end)=startparams_objdef(end)./paramweights(end); end % rot2

% Do optimization with objective function:
fprintf('Measure: %s.  Interpolation: %s.\n',ip.Results.measure,ip.Results.interpmethod);

rigidparams=fminsearch(@objdef2drigid_groupwise, startparams_objdef, optoptions,...
	double(Itarget), double(Itemplate), double(Itarget_bg), ip.Results.bgvalue, N,...
    ip.Results.nodisp, paramweights, ip.Results.measure,...
    ip.Results.interpmethod, paramoffsets, preM, ip.Results.coordshift, scales==1, number_many_slices);

% Grab error record
if nargout>2,
    err_record=objdef2d;
    varargout{1}=-err_record(2:end); % MIiters
end

% Clear persistent variable in objective function
clear objdef2drigid_groupwise.m

% Apply the weightings that the objective function applied - 5 params
rigidparams(1)=rigidparams(1).*paramweights(1)+paramoffsets(1);
rigidparams(2:3)=rigidparams(2:3).*paramweights(2:3)+paramoffsets(2:3);
if numel(rigidparams)==6, rigidparams(6)=rigidparams(6)*paramweights(6)+paramoffsets(6); end
if numel(rigidparams)<4, % ==3 (no scale)
    rigidparams(4:5)=[1 1];
elseif numel(rigidparams)<5, % ==4 (isoscale)
    rigidparams(4:5)=(rigidparams([4 4])'-1).*paramweights([4 4])' + 1;
elseif numel(rigidparams)==5 && scales==1,
    rigidparams(6)=rigidparams(5)*paramweights(5)+paramoffsets(5);
    rigidparams(4:5)=(rigidparams([4 4])'-1).*paramweights([4 4])' + 1;
else % ==5 && scales==2 || >5 (anisoscale)
    rigidparams(4:5)=(rigidparams(4:5)'-1).*paramweights(4:5)' + 1;
end
rigidparams(4:5)=rigidparams(4:5)+paramoffsets(4:5);

rigidparams(4:5) = 1;

%%% Outputs
% Generate the registered image
[Iregd]=imdef(rigidparams, double(Itarget), interpmethod, bgvalue, preM, ip.Results.coordshift);
[IregdLIN,M]=imdef(rigidparams, double(Itarget), 'linear', bgvalue, preM, ip.Results.coordshift);

% Fix data range for cubic interpolation
Iregd(Iregd<0)=0; Iregd(Iregd>(N-1))=N-1;

% Image/volume sizes
[nrows,ncols,nslices,template_attribs]=size(Itemplate);
target_attribs=size(Itarget,4);

% Crop Iregd to Itemplate size
Iregd_crop=imshave(round(Iregd),[nrows ncols]);
Iregd_crop(Iregd_crop<0)=0; Iregd_crop(Iregd_crop>(N-1))=N-1;
% IregdNN_crop=imshave(round(IregdNN),[nrows ncols]);
IregdLIN_crop=imshave(round(IregdLIN),[nrows ncols]);

if ~isempty(ip.UsingDefaults) && ~any(strcmp(ip.UsingDefaults,'coordshift')),
    fprintf('Remember to apply your coordshift when you use rigidparams or M!!!.\n');
end

% Output M and calculate finalMI
if nargout>3,
    varargout{2}=M;
    if nargout>4,
        multiattrib=target_attribs>1 || template_attribs>1;
        if multiattrib,
            Ivect=reshape(Iregd_crop,prod([nrows ncols nslices]),target_attribs);
            Jvect=reshape(Itemplate,prod([nrows ncols nslices]),template_attribs);
            finalMI=CMIgen(double(Ivect),double(Jvect),N);
        else
            finalMI=MI(double(Itemplate),double(IregdLIN_crop),N);
        end
        varargout{3}=finalMI;
        % paramstring=repmat('%7.4f, ',1,numel(rigidparams));
        paramstring=repmat('%.4g, ',1,numel(rigidparams));
        fprintf(['Done: params=[' paramstring '\b\b], finalMI=%g.\n'],rigidparams,finalMI);
    else
        % paramstring=repmat('%7.4f, ',1,numel(rigidparams));
        paramstring=repmat('%.4g, ',1,numel(rigidparams));
        fprintf(['Done: params=[' paramstring '\b\b].\n'],rigidparams);
    end
end

