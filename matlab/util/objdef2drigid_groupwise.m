function err = objdef2drigid_groupwise(rigidparams,I,J,I_bg, fillwith,N,nodisp,paramweights,measure,interpmethod,paramoffsets,preM,coordshift,onescale, number_many_slices)
% objective function for MI based registration - I is target, J is template

persistent err_record regiter

if nargin==0,
    err=err_record(1:regiter);
    clear err_record regiter
    return;
end

if isempty(regiter),
    regiter=0;
    err_record=zeros(1,2e3); % allocate for 2000 iterations
end
regiter=regiter+1;

if nargin<13,
    onescale=false; % says how to interpret rigidparams and paramweights when numel(rigidparams)==5
    if nargin<12 || isempty(coordshift),
        coordshift=zeros(2,1);
        if nargin<11,
            preM=[];
        end
    end
end

% Weight parameters so optimization routine works - REMEMBER THESE VALUES
rot=rigidparams(1)*paramweights(1)+paramoffsets(1);

trans=rigidparams(2:3)'.*paramweights(2:3)'+paramoffsets(2:3)';

% Fix paramoffsets to length 6 if calling routine did not already
if numel(paramoffsets)==4,
    paramoffsets(5:6)=[paramoffsets(4); 0];
elseif numel(paramoffsets)<4,
    paramoffsets(4:6)=0; % paramoffsets=[paramoffsets; zeros(2,1)];
elseif numel(paramoffsets)~=6,
    error('paramoffsets must specify all scales!');
end

if numel(rigidparams)==6, rot2=rigidparams(6)*paramweights(6)+paramoffsets(6);
else rot2=0; end

% Scaling is optional
if numel(rigidparams)<4, % ==3 (no scale)
    scale=[1 1];
elseif numel(rigidparams)<5, % ==4 (isoscale)
    scale=(rigidparams([4 4])'-1).*paramweights([4 4])' + 1;
elseif numel(rigidparams)==5 && onescale,
    scale=(rigidparams([4 4])'-1).*paramweights([4 4])' + 1;
    rot2=rigidparams(5)*paramweights(5)+paramoffsets(5);
else % ==5 && scales==2 || >5 (anisoscale)
    scale=(rigidparams(4:5)'-1).*paramweights(4:5)' + 1;
end
scale=scale+paramoffsets(4:5)';

[Inrows,Incols,DIslices,DIattribs]=size(I);
[nrows,ncols,DJslices,DJattribs]=size(J);

scale = [1 1]; %only for rigid registration, must be removed otherwise. 

% Don't always convert to double, leave in integer to speed things up
Inew=imdef([rot trans scale rot2],double(I),interpmethod,fillwith,preM,coordshift);
Inew=round(Inew); % this should be a noop for 'nearest'

Inew = (Inew + I_bg)/2; %% add the background (an object that needs to be considered in the registrations)

imdisp(Inew);


if strfind(interpmethod,'cubic'), % match *cubic and cubic
    Inew(Inew<0)=0;
    Inew(Inew>(N-1))=N-1;
end

% Same size images - Target (I) may be bigger than template (J)
if any([nrows ncols]~=[Inrows,Incols]); Inew=imshave(Inew,[nrows ncols]); end

switch lower(measure)
    case 'mi'
        err=mimex(J,Inew,N); 
    case 'mmi'
        err = 0;
        addedElem = 1;
        K = number_many_slices;
        k = K;
        while (addedElem < number_many_slices && k>=1)
            if (sumall(J(:,:,k))>0)
                err = err + mimex(J(:,:,k),Inew,N)*exp(-(number_many_slices-k)*(number_many_slices-k)/(2*2));
                addedElem = addedElem + 1;
            end
            k = k - 1;
        end

        if (size(J,3)>K+1) % forward sum
            k = K+2;
            while (k<size(J,3))
                if (sumall(J(:,:,k))>0)
                    err = err + mimex(J(:,:,k),Inew,N)*exp(-(number_many_slices+2-k)*(number_many_slices+2-k)/(2*2));
                    addedElem = addedElem + 1;
                end
                k = k + 1;
            end
        end

        if (sumall(J(:,:,K+1))>0) % last slice has a different target
            err = err +  0.25*addedElem*mimex(J(:,:,number_many_slices+1),Inew,N); % use addedElement to weight the (K+1)th slice to have the same weight
        end
        %{
        err = 0;
        addedElem = 1;
        K = number_many_slices;
        k = K;
        while (addedElem < number_many_slices && k>=1)
            if (sumall(J(:,:,k))>0)
                err = err + mimex(J(:,:,k),Inew,N)*exp(-(number_many_slices-k)*(number_many_slices-k)/(2*2));
                addedElem = addedElem + 1;
            end
            k = k - 1;
        end

        if (size(J,3)>K+1) % forward sum
            k = K+2;
            while (k<size(J,3))
                if (sumall(J(:,:,k))>0)
                    err = err + mimex(J(:,:,k),Inew,N)*exp(-(number_many_slices+2-k)*(number_many_slices+2-k)/(2*2));
                    addedElem = addedElem + 1;
                end
                k = k + 1;
            end
        end

        if (sumall(J(:,:,K+1))>0) % last slice has a different target
            err = err +  0.25*addedElem*mimex(J(:,:,number_many_slices+1),Inew,N); % use addedElement to weight the (K+1)th slice to have the same weight
        end 
           %}

    otherwise
        error('Measure name not recognized for 1D intensity images.');

end



err=-err;
err_record(regiter)=err;

if ~nodisp,
    fprintf('Deforming: ');
    fprintf([repmat('%1.4f ',1,9) ''],[rot trans scale rot2]);
    fprintf('\n');
    % preserve title
    titlestring=get(get(gca,'title'),'String');
    imdisp(J(:,:,number_many_slices)-Inew(:,:,1))
    title(titlestring)
    drawnow
end
