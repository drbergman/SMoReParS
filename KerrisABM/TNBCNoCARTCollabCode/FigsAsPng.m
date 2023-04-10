cd Data/AA_1__AB_20/Vasculature_Images

for ii = 10:10:300
    disp(ii);
    F = openfig(strcat('Figure_',num2str(ii),'.fig'));
    saveas(F,strcat('PFig_',num2str(ii),'.png'));
    close all
end

cd ../../..