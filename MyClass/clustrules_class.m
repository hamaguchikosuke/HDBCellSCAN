classdef clustrules_class < handle
% clustrules class defining clustering rule to be used in HDBCellScan.
% % set the following 
% Compact = (3,default)
% diameter= (8 pixel, default)
% npix_fraclow (0.5, default)
% npix_frachigh (50, default)
% MinClust is a dependent property: pi*(diameter)^2 which will be used in HDBSCAN minimal cluster size parameter.
% ex)
% clear clustrules
% clustrules.Compact=4;
% clustrules.diameter=10;
% clustrules.npix_fraclow=0.1;
% clustrules.npix_frachigh = 40;
% cl = clustrules_class(clustrules);
% 
% Kosuke Hamaguchi 20240302

properties
    Compact = 3;
    diameter = 8;
    npix_fraclow = 0.5;
    npix_frachigh= 50;
    parent =[];
end

properties (Dependent)
    MinClust
    MinNpix
    MaxNpix
end

methods 

    function obj = clustrules_class(input_struct)
        props = properties(obj);
        for ii=1:length(props)
            pname=props{ii};
            try 
                obj.(pname)=input_struct.(pname);
            end
        end
    end

    function MinClust = get.MinClust(obj)
        MinClust = round(pi*(obj.diameter/2)^2);
    end
    function MinNpix = get.MinNpix(obj)
        MinNpix = round(pi/4 * obj.diameter^2 *obj.npix_fraclow);
    end
    function MaxNpix = get.MaxNpix(obj)
        MaxNpix = round(pi/4 * obj.diameter^2 *obj.npix_frachigh);
    end
end


end
