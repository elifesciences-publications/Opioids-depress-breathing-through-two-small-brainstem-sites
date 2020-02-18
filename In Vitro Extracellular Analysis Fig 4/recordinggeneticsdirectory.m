function reclist= recordinggeneticsdirectory(subset)
%this directory returns ids for recordings containing a particular genotype
%background. Often two genotypes are recorded at once (ex: vglut and wt) in
%order to have an internal experimental control, or to increase
%experimental efficiency, so some recids are listed under two genotypes.

switch subset
    
    case 'all'
        reclist={'rec38','rec40','rec54'};
    case 'wt'
        reclist={'rec38','rec40','rec54'};
    case 'vglut2'
        reclist={'rec38'};
    case 'dbx1'
        reclist={'rec54'};
    case 'foxp2'
        reclist={};
    case 'vgat'
        reclist={};
    case 'gad2'
        reclist={'rec40'};
    case 'glyt2'
        reclist={};
    case ''
        reclist={};
        
end
