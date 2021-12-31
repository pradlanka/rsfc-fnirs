function [G, F, df1, df2, p] = mymvgc(Y,Ys,Pmax,grangertype,includeflag, includeZeroLag, robustflag, pcaflag,criterion)
%MYMVGC calculates the multivariate granger causality which assess the lagged 
%information transfer betweeen fNIRS channels 
% Input arguments:
%   Y: long channel data
%   Ys: short channel data, should be [] if no short channel data is present
%    Pmax: model order to be used in the AR models
%   grangertype: can be 'multivariate' or 'bivariate'
%   includeflag: can be 'short' or 'long' or all, 'short' uses short channel data if Ys
%   is non empty, long only uses the long channel data and 'all' uses both
%    Please not for bivariate GC with short channel selected, we try to find 
%   X -> Y, after controlling for the short channels, whereas bivariate GC
%   with long only tries to find the bivariate GC between X -> Y
%   Similarly, multivariate GC with short channel selected, we try to find 
%   X -> Y, after controlling for the short and the remaining long channels, 
%   whereas multivariate GC with long only tries to find the multivariate GC 
%    between X -> Y, after controlling for the remaining long channels
%   includeZeroLag: Includes the zero-lag term, though traditionally not 
%   included in GC, can be logical true or false
%   robustflag: can be true or false. Does robust regression if flag is true
%   pcaflag: Performs PCA on X before regression for the restricted and
%   unrestriced model
%   criterion: IC used for select optimal model order with a max of Pmax
% Output arguments:
%   G: Matrix of Granger causality metric. Calculated as th e log of ratio of MSEu an
%   and MSEr
%   F: F statistic matrix
%   p: p-value corresponding to the F statistic
%   df1: df1 associated with F-statisic, diff of degrees of freedeom
%   between restrcited and unrestricted models
%   df2: df2 associated with F-statisic,


% Sanity checks
if strcmp(grangertype, 'multivariate')
if (isempty(Ys) && ~(strcmp(includeflag,'long')))
    includeflag = 'long';
    disp('Short channel data unvailable. Only using long channel data');
end


% If too many lag terms are included may lead to an error as no effective
% dfs reminaining
% If too many lag terms are included may lead to an error
if (~pcaflag && ~strcmp(includeflag,'short') && (size(Y,1)<=(size(Y,2)+size(Ys,2)+1)*Pmax))
    disp('Max model order too high for GC model. Using PCA for dimentions reduction');
	pcaflag = true;
end
end


[m, n] = size(Y); % [ time x channel ]
[~, ns] = size(Ys);
maxlag =1;

% Center data
Y = bsxfun( @minus , Y , nanmean(Y) );


Y = nirs.math.innovations(Y,Pmax); % Caclulate the innovation terms
Y(1:Pmax,:)=[];

%If short channels available & required in MVGC
if ~strcmp(includeflag,'long')
Ys = bsxfun( @minus , Ys , nanmean(Ys) );
Ys = nirs.math.innovations(Ys,Pmax);
Ys(1:Pmax,:)=[];
end

% Downweighting timepoints for motion if robust flag is selected
if robustflag 
    if ~strcmp(includeflag,'long')
        data = [Y Ys];
    else
        data = Y;
    end
   data = zscore(data);
   nchan = size(data,2);
   w  = wfun(vecnorm(data,2,2)/sqrt(nchan)-(mean(vecnorm(data,2,2)/sqrt(nchan))));
    W = diag(w);

    Y=W*Y; 
    if ~strcmp(includeflag,'long')
     Ys= W*Ys;
    end
end

%% Calculate restricted and unrestricted OLS model fits  
% Initialzing variables
SSEr = nan(n,n); SSEu = nan(n,n); G= nan(n,n); F = nan(n,n); p = nan(n,n);
dfr = nan(n,n); dfu = nan(n,n); df1 = nan(n,n); nPC = nan(n,n); 
for i=1:n
 for j=1:n
     if i~=j
          Yreg = Y(:,i);
          Yreg(1:maxlag) = [];
        % Restricted model terms
        if strcmp(grangertype, 'multivariate')
           if (strcmp(includeflag,'long'))
             XvalsO  = Y(:,setdiff(1:n,[i,j]));
           elseif(strcmp(includeflag,'all'))
              XvalsO = [Y(:,setdiff(1:n,[i,j])),Ys];
           else
               XvalsO = Ys;
           end
           
         if pcaflag
            [~,pc_scores, ~, ~,pc_var,~] = pca(XvalsO ,'Economy',false);
            perc_exp = 0.9; % 90 percent variance explained by the PC
            nPC(i,j) =  find(perc_exp*100.0 <= round(cumsum(pc_var),5),1); 
            XvalsO = pc_scores(:, 1:nPC(i,j));
         end
        no = size(XvalsO,2);
        
        XvalsR = zeros(m-Pmax,no*(maxlag+1));
         
        for in = 1:no    
            inds = (in-1)*maxlag +in : in*(maxlag+1) ;
            XvalsR(:, inds) = nirs.math.lagmatrix( XvalsO(:,in), 0:maxlag);
        end
        if ~includeZeroLag
           XvalsR(:,1:maxlag+1:end) = [];
        end
        XvalsR(1:maxlag,:)= [];
     
        else
         XvalsR = [];
        end
        
        % Unrestricted model extra terms if modified MVGC desired with
        % zero-lag terms included
        XU =  nirs.math.lagmatrix( Y(:,j), 0:maxlag);
        if includeZeroLag
          XU = XU(maxlag+1:end,:);
        else
          XU = XU(maxlag+1:end,2:maxlag+1);
        end
        
        XvalsU = [XU XvalsR];
        
        if robustflag % Robust MVGC
           w= ones(size(XvalsR,1),1);
             W= diag(w);
            for iter =1:10
            [br, statsr] = robustfit(W*XvalsR, W*Yreg);
            [bu, statsu] = robustfit(W*XvalsU, W*Yreg);
             w = min([diag(W), statsr.w, statsu.w],[],2);
             W=diag(w);
            end
            dfu(i,j) = sum(diag(W))-length(bu);
            SSEu(i,j) = (statsu.s^2)*dfu(i,j);  
            dfr(i,j) = sum(diag(W))-length(br);
            SSEr(i,j) = (statsr.s^2)*dfr(i,j);

        else % Traditional MVGC
            br = [ones(m-Pmax-maxlag,1),XvalsR]\ Yreg;
            Er= Yreg - [ones(m-Pmax-maxlag,1),XvalsR]*br;
            SSEr(i,j) = Er'*Er;
            dfr(i,j) = m-Pmax-maxlag-length(br);
            bu = [ones(m-maxlag-Pmax,1),XvalsU]\ Yreg;
            Eu = Yreg - [ones(m-Pmax-maxlag,1),XvalsU]*bu;
        	SSEu(i,j) = Eu'*Eu;
            dfu(i,j) = m-maxlag-Pmax-length(bu);
        end

        % Compute granger stats
        df1(i,j)= dfr(i,j)-dfu(i,j);
        G(i,j)= max(0, log(SSEr(i,j)/dfr(i,j))-log(SSEu(i,j)/dfu(i,j)));% Due to numerical imprecision  
        F(i,j)= max(1e-6,((SSEr(i,j)/SSEu(i,j))-1) .* (dfu(i,j)./df1(i,j))); 
        p(i,j) = 1 - fcdf(F(i,j),df1(i,j),dfu(i,j));
     end
 end
end
df2 = mean(dfu,'all', 'omitnan'); % Get one values instead of a matrix
df1 = mean(df1,'all', 'omitnan'); % Get one values instead of a matrix
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