function [genotype,pathtodata,figpath,surgconds,expconds]=getmousepath(mouseID)
%directory for plethysmograph breath data

%note: injectBA1/2 stands for injection of virus into brain areas 1/ 2
switch mouseID
    
    case 'ei_155_2'
        genotype='oprflox';
        pathtodata='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vivo Example Breathing Data\155_2';
        figpath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vivo Example Breathing Data\155_2\Figures';
        surgconds={'beforeinject','injectBA1','injectBA2'};
        expconds={'Baseline','PBC mOR KO','PBC KF mOR KO'};
    case 'ei_156_30'
        genotype='oprflox';
        pathtodata='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vivo Example Breathing Data\156_30';
        figpath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vivo Example Breathing Data\156_30\Figures';
        surgconds={'beforeinject','injectBA1','injectBA2'};
        expconds={'Baseline','KF mOR KO','PBC KF mOR KO'};
    case 'i_171_0'
        genotype='oprflox';
        pathtodata='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vivo Example Breathing Data\171_0';
        figpath='C:\Users\Iris Bachmutsky\Desktop\Bachmutskyetal.2020\In Vivo Example Breathing Data\171_0\Figures';
        surgconds={'beforeinject','injectBA1'};
        expconds={'Baseline','PBC control Inj'};
    case ''
        genotype='';
        pathtodata='';
        figpath='';
        surgconds={'beforeinject'};
        expconds={'Baseline'};
end

end