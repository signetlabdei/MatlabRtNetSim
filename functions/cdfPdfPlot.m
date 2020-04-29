function cdfPdfPlot(data,Nfig,xAxis,boundsFlag,saveFlag,savePath)

if ~exist("saveFlag","var")
    saveFlag = false;
end
if ~exist("savePath","var")
    savePath = "figs";
end

figure(Nfig)
subplot(2,1,1)
ecdf(data,"Bounds",boundsFlag)
grid on
xlabel("")
subplot(2,1,2)
histogram(data,"Normalization","pdf")
grid on
xlabel(sprintf("x = %s",xAxis))
ylabel("p(x)")

if saveFlag
    savefig(sprintf(savePath+ "/%d.fig",Nfig));
end

end