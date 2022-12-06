function out = pcfCohort(path_to_cohort_file,options,summarize,load_options)


load(path_to_cohort_file,"total_runs","ids")


for si = total_runs:-1:1
    [out_temp,t_temp] = pcfTimeSeries(sprintf("../data/%s",ids(si)),options,summarize,load_options);
    if summarize
        out(si) = struct("t",t_temp,"avg",out_temp.avg,"std",out_temp.std);
    else
        out(si) = struct("t",t_temp,"cell_pcfs",out_temp);
    end
end