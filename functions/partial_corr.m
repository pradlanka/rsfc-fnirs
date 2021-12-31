function [R,Pcorr,dfe] = partial_corr(data,modelorder,includeflag, robustflag,pcaflag)
% PARTIAL_CORR calculates the partial correlation of channel data
% controlling for short seperation, long channels or both.
%
% Inputs: 
%   data: RS-fNIRS channel data
%   modelorder: Model order to be used for the AR model for calcualting
%        innovation terms. Default value is 8*Fs (Sampling rate)
%         Eg: '8x' denotes 8 times sampling rate
%   includeflag: Can be either 'short', 'long' or 'both' to partial out the
%        the effect of either just short channels or remaining long
%        channels or both.
%   robustflag: Can either be true or false. true indicate using robust
%   regression & WLS to estimate parameters to correct for motion artifacts,
%   false uses OLS estimates
%   pcaflag: Can either be true or false. true uses PCA before partial out
%   the effects of the remaining channels, false uses the channel data
%   itself
%
% Outputs:
%   R: Matrix with the partial correlation estimates between channels
%   Pcorr: P-values corresponding to the partial correlation estimates
%   dfe: Effective degrees of freedom 

if(nargin<2 || isempty(modelorder))
    modelorder = round(8*data.Fs);
end
if(nargin<3)
    includeflag='short';
end
if(nargin<4)
    robustflag=true;
end
if(nargin<5)
    pcaflag=true;
end

if(isstr(modelorder))
    Pmax = round(data.Fs*str2num(modelorder(1:strfind(modelorder,'x')-1)));
else
    Pmax = modelorder;
end

SD = data.probe.link;
data = data.data;
data = nirs.math.innovations(data,Pmax);
nchan = size(data,2);
w = ones(size(data,1),1);

% Use a joint weighing function to downweigh motion corrupted timepoints
if robustflag 
data = zscore(data);
w  = wfun(vecnorm(data,2,2)/sqrt(nchan)-(mean(vecnorm(data,2,2)/sqrt(nchan))));
end
W = diag(w);

data=W*data; % Downweight the data 

if  or(ismember('ShortSeperation',SD.Properties.VariableNames),any(SD.ShortSeperation))
    shrtindx = SD.ShortSeperation;
     short= data(:, shrtindx);
     data = data(:, ~shrtindx);
     SD = SD(~shrtindx,:);
else
    includeflag ='long';
    SD.ShortSeperation=[];
end

R=eye(size(data,2)); % Initializing R
P=ones(size(data,2));% Initializing P

[SD,~,lst]=unique(SD);

SD = SD(:,1:2);
for i=1:height(SD)
    for j=i+1:height(SD)
       lst2=(lst==i| lst==j);
        Y=data(:,lst2);
        if(strcmp(includeflag,'short'))
            Ypartial= short; % Partial short sepration channels
        elseif(strcmp(includeflag,'long'))

           Ypartial=data(:,~lst2); % Partial other long channels
        else 
            Ypartial=[data(:,~lst2), short]; % Partial both
        end
        if pcaflag % If PCA falg is turned on first perform PCA before partial correlation
        m = mean(Ypartial,1);
        Ypartial = bsxfun(@minus, Ypartial, m);    
        [U,S,V]=nirs.math.mysvd(Ypartial);
        ncomp = find(cumsum(diag(S.^2))/sum(diag(S.^2))>=0.9, 1 ); 
        % Chosing min no of prinicpal components that explain at least 90% variance
         Ypartial = U(:,1:ncomp);
        end
        if robustflag
            % If robust flag is turned on use robust regression to partial
            % out and use robust correlation to estimate correlation
            % between residuals after partialing out
            beta = []; df  =[];
            for k=1:size(Y,2)
            [beta(:,k), stats] = nirs.math.robustfit(Ypartial,Y(:,k));
             df(:,k) = stats.w;
            end
            w1= min(df,[],2);
            Y = Y-[ones(size(Y,1),1),Ypartial]*beta;
            W1=diag(w1);
            for iter=1:10 % Running multiple times stabilizes the weights
           [R(lst2,lst2),P(lst2,lst2), w2]= nirs.math.robust_corrcoef2(W1*Y);
             W1=min(W1,diag(w2));
            end  
            dfe(lst2,lst2) = sum(min(diag(W),diag(W1)))-size(Ypartial,2)-1; %Effective df
        else
        Ypartial(:,end+1)=1;
        Y = Y-Ypartial*inv(Ypartial'*Ypartial)*Ypartial'*Y;    
        [R(lst2,lst2),P(lst2,lst2)]=corrcoef(Y);   
        dfe(lst2,lst2) = sum(w)-size(Ypartial,2)-1;
        end
    end
end
% Estimating test statistics and p-values
Tstat = R .* sqrt((dfe-2) ./ (1 - R.^2));
Pcorr = 2*tpvalue(-abs(Tstat),dfe-2);
dfe = mean(dfe,'All');
end


function w = wfun(r, tune)
% WFUN computes that weights to downweight outliers as specified by Tukey bisquare loss function.
if nargin<2
    tune = 4.685;
end
    s = mad(r, 0) / 0.6745;
    r = r/s/tune;
    
    w = (1 - r.^2) .* (r < 1 & r > -1);
end

% ------------------------------------------------
% Code used from NIRS toolbox
function p = tpvalue(x,v)
%TPVALUE Compute p-value for t statistic.

normcutoff = 1e7;
if length(x)~=1 && length(v)==1
   v = repmat(v,size(x));
end

% Initialize P.
p = NaN(size(x));
nans = (isnan(x) | ~(0<v)); % v == NaN ==> (0<v) == false

% First compute F(-|x|).
%
% Cauchy distribution.  See Devroye pages 29 and 450.
cauchy = (v == 1);
p(cauchy) = .5 + atan(x(cauchy))/pi;

% Normal Approximation.
normal = (v > normcutoff);
p(normal) = 0.5 * erfc(-x(normal) ./ sqrt(2));

% See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1.
gen = ~(cauchy | normal | nans);
p(gen) = betainc(v(gen) ./ (v(gen) + x(gen).^2), v(gen)/2, 0.5)/2;

% Adjust for x>0.  Right now p<0.5, so this is numerically safe.
reflect = gen & (x > 0);
p(reflect) = 1 - p(reflect);

% Make the result exact for the median.
p(x == 0 & ~nans) = 0.5;

end

