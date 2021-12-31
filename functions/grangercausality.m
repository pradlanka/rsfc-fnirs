function [G, F, df1, df2, p] = grangercausality(data, modelorder, grangertype, includeflag, includeZeroLag, robustflag, pcaflag, criterion)

if(~isempty(strfind(class(data),'.core.Data')))
    Fs=data.Fs;
else
    Fs=1;
end
if(nargin<8)
    criterion='AICc';
end

if(nargin<7)
    pcaflag = true;
end

if(nargin<6)
    robustflag = false;
end

if(nargin<5)
    includeZeroLag =false;
end

if(nargin < 4 && any(data.probe.link.ShortSeperation))
     includeflag ='short';
end

if(nargin<3)
    grangertype = 'multivariate';
end

if(nargin<2 || isempty(modelorder))
    modelorder = round(8*Fs);
end

if(isstr(modelorder))
    Pmax = round(Fs*str2num(modelorder(1:strfind(modelorder,'x')-1)));
else
    Pmax = modelorder;
end

if  any(data.probe.link.ShortSeperation)
    shrtindx = data.probe.link.ShortSeperation;
    short= data.data(:, shrtindx);
    data = data.data(:, ~shrtindx);
else
    short = [];
    data = data.data;
    includeflag ='long';
end

[G, F, df1, df2, p] = mymvgc(data,short,Pmax, grangertype,includeflag, includeZeroLag, robustflag, pcaflag, criterion);
 
end
