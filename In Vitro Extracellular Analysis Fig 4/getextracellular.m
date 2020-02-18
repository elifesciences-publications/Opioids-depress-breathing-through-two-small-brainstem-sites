function [pathtodata,figpath,sumdatapath,datatype,inputchs,qualchs,genotypes,filtered]=getextracellular(recID,varargin)
%this directory returns relevant path and experimental information for
%extracellular recordings of the respiratory rhythm in vitro, in slices
%containing a particular genotype background. Often two genotypes (slices)
%are recorded at once (ex: vglut and wt) in order to have an internal 
%experimental control, and/or to increase experimental efficiency.

switch recID
    
    case 'rec38'
        pathtodata='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\rec38_2019_04_16_vglut2_wt_oprflox_DAMGO\Raw Data'; %path to raw data, in .mat or .txt format
        figpath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\rec38_2019_04_16_vglut2_wt_oprflox_DAMGO\Figures'; 
        sumdatapath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\rec38_2019_04_16_vglut2_wt_oprflox_DAMGO\matFiles'; %folder to store summary data structures
        datatype='*.mat'; %data type for raw data
        inputchs={'c001_Time','c002_IN_0','c003_IN_2'}; %all channels recorded in this session
        qualchs={'c002_IN_0','c003_IN_2'}; %all 'qualified' channels that will be used for analysis
        genotypes={'wt','vglut2'}; %genotypes corresponding to the recordings
        filtered=false; % online paynter filtering, ie data is prefiltered 
    case 'rec40'
        pathtodata='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\rec40_2019_04_30_gad_wt_oprflox_DAMGO\Raw Data';
        figpath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\rec40_2019_04_30_gad_wt_oprflox_DAMGO\Figures';
        sumdatapath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\rec40_2019_04_30_gad_wt_oprflox_DAMGO\matFiles';
        datatype='*.mat';
        inputchs={'c001_Time','c002_Vm_prime','c003_IN_2'};
        qualchs={'c002_Vm_prime','c003_IN_2'};
        genotypes={'wt','gad2'};
        filtered=false;
    case'rec54'
        pathtodata='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\rec54_2019_10_03_dbx1_wt_oprflox_DAMGO\Raw Data';
        figpath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\rec54_2019_10_03_dbx1_wt_oprflox_DAMGO\Figures';
        sumdatapath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vitro Example Extracellular Data\rec54_2019_10_03_dbx1_wt_oprflox_DAMGO\matFiles';
        datatype='*.mat';
        inputchs={'c001_Time','c002_IN_0','c003_IN_2','c004_IN_4','c005_IN_5'};
        qualchs={'c002_IN_0','c003_IN_2'};
        genotypes={'wt','dbx1'};
        filtered=false;
    case ''
        pathtodata='';
        figpath='';
        sumdatapath='';
        datatype='*.mat';
        inputchs={'c001_Time',''};
        qualchs={''};
        genotypes={''};
        filtered=false;
end

end
