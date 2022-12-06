function [t,source_locations,target_locations] = pcfLoadData(path_to_data_file,options)

load(path_to_data_file,"time")
t = time;

load(path_to_data_file,"tumor_locations","immune_locations")
if isstring(options)
    switch options
        case "default"
            source_locations = tumor_locations;
            target_locations = immune_locations;
        case "tumor_to_active"
            source_locations = tumor_locations;

            load(path_to_data_file,"immune_type")
            target_locations = immune_locations(immune_type~=-1,:);

        case "high_antigen_to_immune"
            load(path_to_data_file,"tumor_is_ha")
            source_locations = tumor_locations(tumor_is_ha,:);
            target_locations = immune_locations;

        case "high_antigen_mut_to_immune"
            load(path_to_data_file,"tumor_is_ha","tumor_is_fgfr3_mut")
            source_locations = tumor_locations(tumor_is_ha & tumor_is_fgfr3_mut,:);
            target_locations = immune_locations;

        case "high_antigen_nonmut_to_immune"
            load(path_to_data_file,"tumor_is_ha","tumor_is_fgfr3_mut")
            source_locations = tumor_locations(tumor_is_ha & ~tumor_is_fgfr3_mut,:);
            target_locations = immune_locations;

        case "low_antigen_to_immune"
            load(path_to_data_file,"tumor_is_ha")
            source_locations = tumor_locations(~tumor_is_ha,:);
            target_locations = immune_locations;

        case "low_antigen_mut_to_immune"
            load(path_to_data_file,"tumor_is_ha","tumor_is_fgfr3_mut")
            source_locations = tumor_locations(~tumor_is_ha & tumor_is_fgfr3_mut,:);
            target_locations = immune_locations;

        case "low_antigen_nonmut_to_immune"
            load(path_to_data_file,"tumor_is_ha","tumor_is_fgfr3_mut")
            source_locations = tumor_locations(~tumor_is_ha & ~tumor_is_fgfr3_mut,:);
            target_locations = immune_locations;
    end


end

source_locations = double(source_locations);
target_locations = double(target_locations);
