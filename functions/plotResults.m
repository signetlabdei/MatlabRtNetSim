function plotResults(results,dontPlot)

if ~exist("dontPlot","var")
    dontPlot = "";
end

resFields = string(fields(results));
for i = 1:length(resFields)
    field = resFields(i);
    if any(field == dontPlot)
        continue
    end
    
    data = getData(results,field);
    if ~isvector(data) ||...
        all(isnan(data))
        continue
    end
    
    % plot
    figure
    ax(1) = subplot(2,1,1);
    ecdf(data,"bounds","on")
    
    ax(2) = subplot(2,1,2);
    histogram(data,"Normalization","pdf"); hold on
    ksdensity(data)
    
    ylabel("p(x)")
    xlabel("x = " + strrep(field,"_","\_"))
    linkaxes(ax,'x')
end

end


function data = getData(res,field)
data = {res.(field)};
data = cell2mat(data(:));
end