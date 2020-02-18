function mouselist= mousegeneticsdirectory(subset)

switch subset
    
    case 'all'
        mouselist={'ei_155_2','ei_156_30','i_171_0'};
    case'pbckos'
        mouselist={'ei_155_2'};
    case 'kfkos'
        mouselist={'ei_156_30'};
    case 'pbckfkos'
        mouselist={'ei_155_2','ei_156_30'};
    case 'cohort1'
        mouselist={'ei_155_2'};
    case 'cohort2'
        mouselist={'ei_156_30'};
    case 'controls'
        mouselist={'i_171_0'};
end

end
