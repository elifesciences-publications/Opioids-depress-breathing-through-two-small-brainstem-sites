function [wCoefs]=cwt(data,scales,wavelet)
%this function is based on MassSpecWavelet code written by and based on
%"Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching "
%(Du et al 2006);

if strcmp(wavelet,"mexh")
    %generate continuous mexican hat wavelet
    psi_xval=linspace(-8,8,1024);
    mexihatfunc=(2/sqrt(3)*pi^(-0.25))*(1-psi_xval.^2).*exp(-psi_xval.^2/2);
else 
    error("Other wavelet formats not currently supported cuz Iris didn't have time for that ish! Mexh only!")
end

datalen=length(data);
base=2;nLevel=nan;
addLength=extendNBase(data,nLevel,base);
[data,padcond]=extendlength(data,addLength);
newdatalen=length(data);
nbscales=length(scales);
psi_xval=psi_xval-psi_xval(1);%now 0:16
dxval=psi_xval(2); %sampling interval
xmax=psi_xval(end);

wCoefs=zeros(datalen,length(scales));
for i=1:length(scales)
    scale=scales(i);
    j= 1 + floor((0:(scale * xmax))/(scale * dxval)); %ie numsamps=scale*16+1, subsamples mexican hat wavelet (numsamps will be equivalent to timescale of wavelet)
    if length(j)==1
        j=[1 1];
    end
    lenWave=length(j);
    f(1:lenWave)=fliplr(mexihatfunc(j))-mean(mexihatfunc(j));
    if length(f)>newdatalen
        disp("scale is too large")
        break
    end
    temp_wCoefs=conv(data,f,'same'); %idk 
    padding=addLength/2;
    if padcond=="bilateral"
        temp_wCoefs=temp_wCoefs((padding+1):(length(temp_wCoefs)-padding));
    elseif padcond=="unilateral"
        temp_wCoefs=temp_wCoefs((addLength+1):end);
    end
    wCoefs(:,i)=temp_wCoefs;
end

end

function [addLength]=extendNBase(data,nLevel,base)
nRow=length(data);

if isnan(nLevel)
    nR1=nextn(nRow,base);
else
    error("other formats not available yet cuz iris didn't have time for this")
end

addLength=nR1-nRow;
end


function [nR1]=nextn(nRow,base)
    P=nextpow2(nRow);
    nR1=2^P;
end

function [data,cond]=extendlength(data,addLength)
%padd data with zeros

data=reshape(data,[length(data),1]); %turn into collumn vector
padlen=addLength/2;
if floor(padlen)==padlen %bidirectional pad if padlen is an int
    data=vertcat(zeros(padlen,1),data,zeros(padlen,1));
    cond="bilateral";
else
    data=vertcat(zeros(addLength,1),data);
    cond="unilateral";
end
end
