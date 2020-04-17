function ZeldaDataAnalysis


%Analysis and plots for the Zelda experiment. The idea is to reproduce
%the plots from the Current Biology paper for the wild-type and Zelda
%germline clones.


%% Load the data and set up parameters

%Parameters:
MinEmbryos=3;       %Minimum number of embryos to be able to do an average
MinParticles=3;     %Minimum number of particles to be able to do an average

%This is Ana's data taken at 55uW. Do they all have the same rotation? Do
%we want to retake this all at lower power, which is the standard now?
DataEve2=LoadMS2Sets('eve2_HMM_noML_35uW');
%DataP2Pv5=LoadMS2Sets('P2P-MS2v5-lacZ-36uW');

%% Checking sample quality

%Wild-type

%Check the rotation settings 
[DataEve2.ImageRotation]
%From these rotations, it looks like we should expect less than 3%
%differences between the current data sets. I think we're good to go.


%% Integrated fluorescence in each nc

close all

%Plot the integrated fluroescence for each data set

%Eve2
[TotalProdEve2,TotalProdErrorEve2,TotalProdNEve2,...
    MeanTotalProdEve2,SDTotalProdEve2,SETotalProdEve2]=...
    IntegratemRNA(DataEve2,MinParticles,MinEmbryos);
% 
% %nc13
% nc=13;
% figure(1)
% hold all
% LegendLabel={};
% for i=1:length(DataEve2)
%     errorbar(DataEve2(i).APbinID,TotalProdEve2(i,:,nc),...
%         TotalProdErrorEve2(i,:,nc),'.-')
%     LegendLabel={LegendLabel{:},DataEve2(i).SetName(10:end-1)};
% end
% hold off
% box on
% title(['nc',num2str(nc)])
% legend(LegendLabel)

%nc14
nc=14;
figure(2)
hold all
LegendLabel={};
for i=1:length(DataEve2)
    errorbar(DataEve2(i).APbinID,TotalProdEve2(i,:,nc),...
        TotalProdErrorEve2(i,:,nc),'.-')
    LegendLabel={LegendLabel{:},DataEve2(i).SetName(10:end-1)};
end
hold off
box on
title(['nc',num2str(nc)])
legend(LegendLabel)

% 
% %P2Pv5
% [TotalProdP2Pv5,TotalProdErrorP2Pv5,TotalProdNP2Pv5,...
%     MeanTotalProdP2Pv5,SDTotalProdP2Pv5,SETotalProdP2Pv5]=...
%     IntegratemRNA(DataP2Pv5,MinParticles,MinEmbryos);
% 
% %nc13
% nc=13;
% figure(3)
% hold all
% LegendLabel={};
% for i=1:length(DataP2Pv5)
%     errorbar(DataP2Pv5(i).APbinID,TotalProdP2Pv5(i,:,nc),...
%         TotalProdErrorP2Pv5(i,:,nc),'.-')
%     LegendLabel={LegendLabel{:},DataP2Pv5(i).SetName(10:end-1)};
% end
% hold off
% box on
% title(['nc',num2str(nc)])
% legend(LegendLabel)
% 
% %nc14
% nc=14;
% figure(4)
% hold all
% LegendLabel={};
% for i=1:length(DataP2Pv5)
%     errorbar(DataP2Pv5(i).APbinID,TotalProdP2Pv5(i,:,nc),...
%         TotalProdErrorP2Pv5(i,:,nc),'.-')
%     LegendLabel={LegendLabel{:},DataP2Pv5(i).SetName(10:end-1)};
% end
% hold off
% box on
% title(['nc',num2str(nc)])
% legend(LegendLabel)
% 
% 

% %Compare P2P vs. P2Pv5
% 
% 
% %nc14
% nc=14;
% figure(5)
% PlotHandle=errorbar(DataEve2(1).APbinID,MeanTotalProdEve2(:,nc),SETotalProdEve2(:,nc),'.-k');
% hold on
% PlotHandle(end+1)=errorbar(DataP2Pv5(1).APbinID,MeanTotalProdP2Pv5(:,nc),SETotalProdP2Pv5(:,nc),'.-r');
% hold off
% title(['nc',num2str(nc)])
% legend('P2P-MS2','P2P-MS2v5')
% xlim([0.1,0.75])
% xlabel('AP position')
% ylabel('Integrated mRNA')
% StandardFigure(PlotHandle,gca)
% 
% %nc13
% nc=13;
% figure(6)
% PlotHandle=errorbar(DataEve2(1).APbinID,MeanTotalProdEve2(:,nc),SETotalProdEve2(:,nc),'.-k');
% hold on
% PlotHandle(end+1)=errorbar(DataP2Pv5(1).APbinID,MeanTotalProdP2Pv5(:,nc),SETotalProdP2Pv5(:,nc),'.-r');
% hold off
% title(['nc',num2str(nc)])
% legend('P2P-MS2','P2P-MS2v5')
% xlim([0.1,0.75])
% xlabel('AP position')
% ylabel('Integrated mRNA')
% StandardFigure(PlotHandle,gca)
% 
% 
% 
%% Mean vector AP

figure(1)
plot(DataEve2(3).ElapsedTime,DataEve2(3).MeanVectorAP(:,10),'.-k')



%% Rate of initiation
%close all


%Average the data from multiple embryos - Maybe do this inside LoadMS2Sets

%Wild-type P2P
%Pull out the information from individual embryos
IndividualRateP2P=nan(length(DataEve2),length(DataEve2(1).APbinID),3);
IndividualSERateP2P=nan(length(DataEve2),length(DataEve2(1).APbinID),3);


for i=1:length(DataEve2)
   for j=1:length(DataEve2(1).APbinID)
        for nc=1:3   %This is nc12 - nc14
            if DataEve2(i).MeanFits(j,nc).Approved>0
                IndividualRateP2P(i,j,nc)=DataEve2(i).MeanFits(j,nc).RateFit;
                IndividualSERateP2P(i,j,nc)=DataEve2(i).MeanFits(j,nc).SDRateFit;
            end
       end
   end
end


%Do a weighted average
%Which AP positions have more than MinEmbryos data points?
FilterFits=sum(~isnan(IndividualRateP2P))>=MinEmbryos;
RateP2P=nansum(IndividualRateP2P.*(1./IndividualSERateP2P).^2)./...
    nansum((1./IndividualSERateP2P).^2);
RateP2P(~FilterFits)=nan;
RateP2P=squeeze(RateP2P);

ErrorRateP2P=1./nansum((1./IndividualSERateP2P).^2);
ErrorRateP2P(~FilterFits)=nan;
ErrorRateP2P=squeeze(ErrorRateP2P);



%Look at individual WT sets
figure(1)
clf
nc=2;
hold all
LegendLabel={};
PlotHandle=[];
for i=1:length(DataEve2)
    PlotHandle=[PlotHandle,errorbar(DataEve2(1).APbinID,IndividualRateP2P(i,:,nc),...
        IndividualSERateP2P(i,:,nc),'.-')];
    LegendLabel={LegendLabel{:},DataEve2(i).SetName(10:end-1)};
end
legend(LegendLabel)
xlabel('AP position')
ylabel('Rate of initiation (au/min)')
hold off
box on
xlim([0.2,0.8])
ylim([0,600])
StandardFigure(PlotHandle,gca)



