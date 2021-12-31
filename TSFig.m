% This code generates Fig 2 in the manuscript which shows the timeseries
% and the channel covariance
 
Curr_Folder = pwd; % Get path of working directory

% This code requires NIRS Brain AnalyzIR toolbox to run. The toolbox can be
% downloaded from  https://github.com/huppertt/nirs-toolbox. 


% Please download the toolbox and add to path

addpath(genpath('nirs-toolbox')); % Add toolbox to path

addpath(genpath('functions')); % Add to path folder containing important functions]

% setting parameters
pmax=10;  % model order to use for the AR model
t = (0:1/4:300)'; % Fs  = 4 Hz
probe =defaultProbe();
lags=true(1); % Zero lag relationshop
SNR1 = 1;  % ratio of brain signals to skin
SNR2 = 100;  % ratio of  brain signal to random noise

% Generating non-spatial noise
raw0 =nirs.testing.simARNoise_shortsep(probe,t, pmax, 0);
e0=raw0 .data;
% Generating spatially correlated global noise
raw150=nirs.testing.simARNoise_shortsep(probe,t, pmax, 150);
e150=raw150.data;
 
%%  Simulating the connectivity part
types=probe.types;
et=zeros(length(t),height(probe.link));

for id=1:length(types)
    lst=find(~probe.link.ShortSeperation & ...
        probe.link.type==types(id));

        flag = 1;
        while(flag~=0)
         nrow = length(lst);
         truth = repmat(zeros(nrow),1,1,length(lags));
           for nlag = 1:length(lags)
            if(lags(nlag))
            fract=.05; % Percentage of true positives 

            rand_vec1 = rand(0.5*nrow*(nrow-1),1)> (1-fract);
            truth1  = zeros(nrow);
            truth1(logical(triu(ones(nrow),1))) = rand_vec1;
            rand_vec2 = rand(0.5*nrow*(nrow-1),1)> (1-fract);
            truth1(logical(tril(ones(nrow),-1))) = rand_vec2;
            truth(:,:,nlag)=truth1;
           end
           end
         truth(:,:,1) = triu(truth(:,:,1))+triu(truth(:,:,1))'; 

         sum_mtx = ones(length(lst))+ (ones(length(lst),1)*(sum(truth+permute(truth,[2,1,3]),[2,3]))'); 
         truth(:,:,1) = 2*truth(:,:,1)./sqrt((sum_mtx).*(sum_mtx(1,:)'*(ones(1,length(lst)))));
         truth(:,:,2:end)  = truth(:,:,2:end)./sqrt((sum_mtx).*(sum_mtx(1,:)'*(ones(1,length(lst)))));


       for n = 0:2*length(lags)-2
           if rem(n,2)==0
               truth_mtx = triu(truth(:,:,n/2+1));
           else
               truth_mtx = tril(truth(:,:,n));
           end
            [~,flag1]=chol(truth_mtx+truth_mtx'+eye(length(lst)));
            flag = and(flag, flag1);
       end
       end
  
                    
        % Zero lag relationship    
    et(:,lst)  = mvnrnd(zeros(length(lst),1),eye(length(lst)), length(t) );
    
    % Adding lagged relationships if required
    for ilag=1:length(lags) 
        if(lags(ilag))
            for i=1:length(lst)-1
                for j=i+1:length(lst)
                    if truth(i,j,ilag)
                     e1=randn(length(t),1);
                      et(:,lst(i))= et(:,lst(i))+e1;
                      et(:,lst(j))= et(:,lst(j))+ circshift(e1, -(ilag-1));
                    end
                    if truth(j,i,ilag)
                       e2=randn(length(t),1);
                       et(:,lst(j))= et(:,lst(j))+e2;
                       et(:,lst(i))= et(:,lst(i))+ circshift(e2, -(ilag-1));
                    end
                 end
            end
         end 
    end
     % Passing the data through an AR filter
    a = randAR( pmax );
    for l = 1:length(lst)
        et(:,lst(l)) = filter(1, [1; -a], et(:,lst(l)));
    end
end

% Scaling the neural components according to specified SNR
lst=find(~probe.link.ShortSeperation);
et(:,lst)=SNR1*6*et(:,lst)./(ones(length(t),1)*std(et(:,lst),[],1));

% Scaling random white noise according to specified SNR
ez = randn(length(t),height(probe.link));
ez = (SNR1/SNR2)*6*ez./(ones(length(t),1)*std(ez,[],1));

% Total signal  (Summation of neural connectivity component, global spatial
% component & white noise)
y0= (e0+et+ez)+200;

y0=y0./(ones(size(y0,1),1)*mean(y0,1));

% Simulated data containing just temporally autocorrelated noise
raw1 = nirs.core.Data();
raw1.data  = 100 * exp( -y0);
raw1.probe  = probe;
raw1.time   = t;


% Simulated data containing temporally autocorrelated noise & spatial
% correlated global noise
y150= (e150+et+ez)+200;

y150=y150./(ones(size(y150,1),1)*mean(y150,1));

raw2 = nirs.core.Data();
raw2.data   = 100 * exp( -y150);
raw2.probe  = probe;
raw2.time   = t;

% Adding motion artifacts to create a simulated data with temporally autocorrelated noise 
% & spatial correlated global noise & motion artifacts.

raw3 = nirs.testing.simMotionArtifact(raw2); 

Fs = raw1.Fs; % Sampling rate

% % Visualizaing raw data for the three conditions
title1 = 'W temporal autocorr but wo global phys noise or motion artifacts';
figure('Name',title1,'NumberTitle','off','Position',[100 100 1200 200]); plot(raw1.data((1*Fs:round(Fs):200*Fs), 5),'LineWidth',2);ylim([30,45]);
hold on;   yline(35, 'r--', 'LineWidth', 3); hold on;   yline(40, 'r--', 'LineWidth', 3);
xlabel('Time (t) in sec'), ylabel('Intensity'); ax =gca; ax.FontSize = 16; ax.FontWeight = 'bold';
title2 = 'W temporal autocorr and global phys noise but wo motion artifacts';
figure('Name',title2,'NumberTitle','off','Position',[100 100 1200 200]); plot(raw2.data((1*Fs:round(Fs):200*Fs), 5),'LineWidth',2);ylim([30,45]);
hold on;   yline(35, 'r--', 'LineWidth', 3); hold on;   yline(40, 'r--', 'LineWidth', 3);
xlabel('Time (t) in sec'), ylabel('Intensity'); ax =gca; ax.FontSize = 16; ax.FontWeight = 'bold';
title3 = 'W temporal autocorr, global phys noise and motion artifacts';
figure('Name',title3,'NumberTitle','off', 'Position',[100 100 1200 200]); plot(raw3.data((1*Fs:round(Fs):200*Fs), 5),'LineWidth',2);ylim([0,80]);
hold on;   yline(30, 'r--', 'LineWidth', 3); hold on;   yline(35, 'r--', 'LineWidth', 3);
xlabel('Time (t) in sec'), ylabel('Intensity'); ax =gca; ax.FontSize = 16; ax.FontWeight = 'bold';

% Visualzing the channel correlation matrix
figure; imagesc(corr(raw1.data), [-1,1]); hold on; colorbar; ax =gca; ax.FontSize = 15; ax.FontWeight = 'bold';
figure; imagesc(corr(raw2.data), [-1,1]); hold on; colorbar; ax =gca; ax.FontSize = 15; ax.FontWeight = 'bold';
figure; imagesc(corr(raw3.data), [-1,1]); hold on; colorbar; ax =gca; ax.FontSize = 15; ax.FontWeight = 'bold';




function a = randAR( P )
% random Pth order AR coef
a = flipud( cumsum( rand(P, 1) ) );
a = a / sum(a) * 0.99;
end

function probe = defaultProbe(lambda)
% Default probe generates the default probe for simulating data
if(nargin==0)
    lambda=[690 830];
end


srcPos(:,1) = (-80:20:80)';
srcPos(:,2:3) = 0;

detPos(:,1) = (-70:20:70)';
detPos(:,2) = 25;
detPos(:,3) = 0;
detPos=[detPos; srcPos-ones(size(srcPos,1),1)*[-5 5 0]];

probe = nirs.core.Probe(srcPos,detPos);

link = [1	1	0
    2	1   0
    2	2   0
    3	2   0
    3	3   0
    4	3   0
    4	4   0
    5	4	0
    5	5	0
    6	5	0
    6	6	0
    7	6	0
    7	7	0
    8	7	0
    8	8	0
    9	8	0
    1   9   1
    2   10  1
    3   11  1
    4   12  1
    5   13  1
    6   14  1
    7   15  1
    8   16  1
    9   17  1];

link=[repmat(link,length(lambda),1) reshape(repmat(lambda(:)',size(link,1),1),1,[])'];

link = sortrows(link);
short = false(size(link,1),1);
short(end-size(srcPos,1)+1:end)=true;
probe.link = table(link(:,1), link(:,2), link(:,4),(link(:,3)==1), ...
    'VariableNames', {'source', 'detector', 'type','ShortSeperation'});


end

function [r,T]=mvnrnd(mu,sigma,cases)
% Generates samples from Multivariate normal distribution

[T,err] = cholcov(sigma);
r = randn(cases,size(T,1)) * T + ones(cases,1)*mu';

end


