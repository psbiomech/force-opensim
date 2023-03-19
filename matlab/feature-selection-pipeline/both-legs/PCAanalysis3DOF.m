function PCAanalysis3DOF

warning off;
outputDatabase      = 'D:/ACLR_project/outputDatabase3DOF/';
if exist([outputDatabase 'PCAanalysis3DOF_inter.mat'],'file')
%     aux             = load([outputDatabase 'PCAanalysis3DOF_inter.mat'],...
%                            'trainingdata','varlabels','isACLR','weights',...
%                            'massheight','pcacoeffs');
%     trainingdata    = aux.trainingdata;
%     varlabels       = aux.varlabels;
%     isACLR          = aux.isACLR;
%     weights         = aux.weights;
%     pcacoeffs       = aux.pcacoeffs;
% %     massheight      = aux.massheight;
%     clear aux;
    load([outputDatabase 'PCAanalysis3DOF_inter.mat']);
else
    % Check whether the data has been interpolated or not
    if ~exist([outputDatabase 'experimentCTRLvsACLR3DOF_interp.mat'],'file')
        experimentCTRLvsACLR3DOF;
    end
    
    figuredir           = 'C:/Users/mariom/Documents/MATLAB/GaitProject/testgraphs/';
    aux                 = load([outputDatabase 'experimentCTRLvsACLR3DOF_interp.mat']);
    dataCTRL            = aux.dataCTRL;
    idxCTRL             = aux.idxCTRL;
    dataACLR            = aux.dataACLR;
    idxACLR             = aux.idxACLR;
    labels              = aux.labels;
    clear aux;
    weigths             = [calculateWeights(idxCTRL(:,2));  % Estimating the weights of each trial
                           calculateWeights(idxACLR(:,2))]; % based on the number of trials per subject
    weights             = weigths./sum(weigths);
    
    aux                 = load('C:\Users\mariom\Documents\MATLAB\GaitProject\toolbox\subjectmassheight.mat');
    subjectmassheight   = aux.subjectmassheight;
    clear aux;
    massheight          = zeros(length(weights),4);
    massheight(:,1:2)   = [idxCTRL(:,1:2);
                           idxACLR(:,1:2)];
    for i=1:size(subjectmassheight,1)
        idx                 = all(bsxfun(@eq,massheight(:,1:2),subjectmassheight(i,1:2)),2);
        massheight(idx,3:4) = repmat(subjectmassheight(i,3:4),sum(idx),1);
    end
    massheight          = massheight(:,3:4);
    
    dataCTRL.F_muscle   = dataCTRL.F_muscle(2:8,:,:);
    dataCTRL.phi_joint  = dataCTRL.phi_joint(1:13,:,:);%([1:7 11],:,:);
    dataCTRL.M_joint    = dataCTRL.M_joint([1:7 11:13],:,:);%(1:7,:,:);
    dataACLR.F_muscle   = dataACLR.F_muscle(2:8,:,:);
    dataACLR.phi_joint  = dataACLR.phi_joint(1:13,:,:);%([1:7 11],:,:);
    dataACLR.M_joint    = dataACLR.M_joint([1:7 11:13],:,:);%(1:7,:,:);
    labels.muscleLabels = labels.muscleLabels(2:8);
    labels.momentLabels = labels.angleLabels([1:7 11:13]);%(1:7);
    labels.angleLabels  = labels.angleLabels(1:13);%([1:7 11]);
    downtime            = linspace(0,100,100);
    fields              = fieldnames(dataCTRL);
    fields              = fields(3:5);
    isACLR              = [false(size(idxCTRL,1),1);
        true(size(idxACLR,1),1)];
    ylabelsplots        = {'Force (BW)','Angle (deg)','Moment (BW-HT)'};
    totalvalidscores    = 0; % Total number of features that make the difference in the data
    pcexplained         = {};
    inc                 = 1;
    for i=1:3
        for j=1:size(dataCTRL.(fields{i}),1)
            rawdata                     = [squeeze(dataCTRL.(fields{i})(j,:,:))';
                                           squeeze(dataACLR.(fields{i})(j,:,:))'];
            codedata                    = zeros(size(rawdata,1),length(downtime));
            for k=1:size(rawdata,1)
                % Downsampling the data
                codedata(k,:)           = interp1(dataCTRL.time,rawdata(k,:),downtime);
            end
            % Extracting the Principal Components
            [coeff,score,~,~,explained] = pca(codedata,'Weights',weights);
            switch i
                case 1
                    infieldlabel        = labels.muscleLabels{j};
                case 2
                    infieldlabel        = labels.angleLabels{j};
                case 3
                    infieldlabel        = labels.momentLabels{j};
            end
            
            % % Finding at which score can we represent 95% of the variance
            % numvalidscores                       = find(cumsum(explained)>=95,1,'first');
            
            % Using parallel analysis to identify the number of components
            % to keep
            if ~exist([outputDatabase 'PAdata.mat'],'file')
                [PAmeans,PApercentiles] = parallelanalysis(size(codedata,2),size(codedata,1));  %#ok<ASGLU>
                save([outputDatabase 'PAdata.mat'],'PAmeans','PApercentiles');
            else
                aux           = load([outputDatabase 'PAdata.mat'],'PApercentiles');
                PApercentiles = aux.PApercentiles;
                clear aux;
            end
            
            latent                               = sort(eig(weightedcorrs(codedata,weights)),'descend');
            numvalidscores                       = find(latent>PApercentiles,1,'last');
            disp(numvalidscores);
            pcexplained{inc,1}                   = [fields{i} '_' infieldlabel]; %#ok<*AGROW>
            for bb=1:numvalidscores
                pcexplained{inc,1+bb}            = explained(bb)./100;
            end
            inc                                  = inc+1;
            
            pcacoeffs.(fields{i}).(infieldlabel) = coeff;
            pcascores.(fields{i}).(infieldlabel) = score(:,1:numvalidscores);
            pcaquantl.(fields{i}).(infieldlabel) = quantile(pcascores.(fields{i}).(infieldlabel),[.025 .25 .50 .75 .975]);
            if numvalidscores>1
                pcatopidx.(fields{i}).(infieldlabel) = bsxfun(@le,pcascores.(fields{i}).(infieldlabel),...
                                                                  pcaquantl.(fields{i}).(infieldlabel)(2,:)); % This is 25%
                pcabotidx.(fields{i}).(infieldlabel) = bsxfun(@ge,pcascores.(fields{i}).(infieldlabel),...
                                                                  pcaquantl.(fields{i}).(infieldlabel)(4,:)); % This is 75%
            else
                pcatopidx.(fields{i}).(infieldlabel) = bsxfun(@le,pcascores.(fields{i}).(infieldlabel),...
                                                                  pcaquantl.(fields{i}).(infieldlabel)(2));
                pcabotidx.(fields{i}).(infieldlabel) = bsxfun(@ge,pcascores.(fields{i}).(infieldlabel),...
                                                                  pcaquantl.(fields{i}).(infieldlabel)(4));
            end
            totalvalidscores                     = totalvalidscores + numvalidscores;
            % Generating graphs
            for k=1:numvalidscores
                figurename             = [figuredir fields{i} '_' infieldlabel '_PC_' num2str(k) '.png'];
                if ~exist(figurename,'file')
                    % If the graph has not been created before...
                    meanlo             = mean(rawdata(pcabotidx.(fields{i}).(infieldlabel)(:,k),:),1);
                    stdlo              =  std(rawdata(pcabotidx.(fields{i}).(infieldlabel)(:,k),:),[],1);
                    meanhi             = mean(rawdata(pcatopidx.(fields{i}).(infieldlabel)(:,k),:),1);
                    stdhi              =  std(rawdata(pcatopidx.(fields{i}).(infieldlabel)(:,k),:),[],1);
                    plotlabels.title   = [fields{i} ' ' infieldlabel ' PC ' num2str(k)];
                    plotlabels.ylabel  = ylabelsplots{i};
                    plotlabels.legends = {'25%Q','75%Q'};
                    plotGroups(dataCTRL.time,meanlo,stdlo,meanhi,stdhi,plotlabels,i==1);
                    saveas(gcf,figurename,'png');
                end
            end
        end
    end
    
    % Extracting the trainingdata matrix
    startcol            = 1;
    trainingdata        = zeros(length(isACLR),totalvalidscores);
    varlabels           = cell(totalvalidscores,1);
    for i=1:3
        for j=1:size(dataCTRL.(fields{i}),1)
            switch i
                case 1
                    infieldlabel        = labels.muscleLabels{j};
                case 2
                    infieldlabel        = labels.angleLabels{j};
                case 3
                    infieldlabel        = labels.momentLabels{j};
            end
            endcol                          = startcol + size(pcascores.(fields{i}).(infieldlabel),2) - 1;
            trainingdata(:,startcol:endcol) = pcascores.(fields{i}).(infieldlabel);
            aux                             = 1;
            for k=startcol:endcol
                varlabels{k}   = [fields{i} ' ' infieldlabel ' PC' num2str(aux)];
                aux            = aux+1;
            end
            startcol                        = endcol + 1;
        end
    end
    % Storing the results
    save([outputDatabase 'PCAanalysis3DOF_inter.mat'],'pcacoeffs','pcascores',...
         'pcaquantl','pcatopidx','pcabotidx','trainingdata','varlabels','isACLR',...
         'weights','massheight','pcexplained');
end

if isempty(who('-file',[outputDatabase 'PCAanalysis3DOF_inter.mat'],'idxMAT'))
    % Second approach: using other techniques such as forward selection to
    % minimize the data. I would produce interpretable data at least in
    % principle.
    Signif              = IndFeat(trainingdata,isACLR);
    idx                 = Signif>=2.0; % Find out about the treshold
    numVar              = find(idx);   % Find the indexes of the valid variables
    opts                = statset('display','iter');
    % Backup of classifiers for forward selection
    %fun                 = @(XT,yT,Xt,yt)(sum(~eq(yt,classify(Xt,XT,yT,'mahalanobis'))));
    %fun                 = @(XT,yT,Xt,yt)(sum(~eq(yt,predict(ClassificationKNN.fit(XT,yT,'NSMethod','exhaustive','Distance','mahalanobis'),Xt))));
    %fun                 = @(XT,yT,Xt,yt)(sum(~eq(yt,svmclassify(svmtrain(XT,yT,'kernel_function','rbf'),Xt))));
    %fun                 = @(XT,yT,Xt,yt)(sum(~eq(yt,predict(NaiveBayes.fit(XT,yT),Xt)))); % Using Naive Bayes
    %fs                  = sequentialfs(fun,trainingdata(:,idx),isACLR,'cv',c,'NFeatures',10,...
    %                                   'Direction','forward','options',opts);
    %idx(numVar(~fs))    = 0;           % Remove those variables which are not valid
    %numVar              = find(idx);   % Find the indexes of the valid variables
    numIter             = 1000;
    idxMAT              = repmat(idx,numIter,1);
    fun                 = @(XT,yT,Xt,yt)(sum(~eq(yt,predict(NaiveBayes.fit(XT,yT),Xt)))); % Using Naive Bayes
    for i=1:numIter
        % Repeating the experiment 100 times to reduce the variability due to
        % the random initialization of the Naive Bayes classifier
        disp(['----> Start of the iteration No. ' num2str(i)]);
        c                     = cvpartition(isACLR,'k',10);
        fs                    = sequentialfs(fun,trainingdata(:,idx),isACLR,'cv',c,...
                                             'NFeatures',10,'Direction','forward',...
                                             'options',opts);
        idxMAT(i,numVar(~fs)) = 0; % Remove those variables which are not valid
        disp(['----> End of the iteration No. ' num2str(i)]);
    end
    timesSelected       = sum(idxMAT,1); % Times that each variable was selected
    [~,sortSel]         = sort(timesSelected);
    numVar              = sort(sortSel(end-9:end));
    statsignif          = zeros(length(numVar),7);
    for i=1:length(numVar)
        statsignif(i,1)                   =   mean(trainingdata(~isACLR,numVar(i)));
        statsignif(i,2)                   =    std(trainingdata(~isACLR,numVar(i)));
        statsignif(i,3)                   =   mean(trainingdata( isACLR,numVar(i)));
        statsignif(i,4)                   =    std(trainingdata( isACLR,numVar(i)));
        [statsignif(i,5),statsignif(i,6)] = ttest2(trainingdata(~isACLR,numVar(i)),...
                                                   trainingdata( isACLR,numVar(i)));
        mesvals                           = mes(trainingdata(~isACLR,numVar(i)),...
                                                trainingdata( isACLR,numVar(i)),...
                                                'hedgesg');
        statsignif(i,7)                   = mesvals.hedgesg;
    end
    T2                  = pdist2(trainingdata(:,numVar),zeros(1,10),'mahalanobis');
    X                   = linspace(0,10,1000);
    Fctrl               = ksdensity(T2(~isACLR),X);
    Faclr               = ksdensity(T2( isACLR),X);
    plot(X,Fctrl,X,Faclr,'LineWidth',2);
    save([outputDatabase 'PCAanalysis3DOF_inter.mat'],'-append','Signif','idx','idxMAT',...
         'numVar','statsignif','T2','X','Fctrl','Faclr');
else
    aux        = load([outputDatabase 'PCAanalysis3DOF_inter.mat'],'numVar','pcabotidx',...
                      'pcatopidx','statsignif');
    numVar     = aux.numVar;
    pcabotidx  = aux.pcabotidx;
    pcatopidx  = aux.pcatopidx;
    statsignif = aux.statsignif;
    clear aux;
end

%if exist('C:/Users/mariom/Documents/MATLAB/GaitProject/pca_graphs/main_10_features_pca_comps.fig','file')==0
    aux                 = load([outputDatabase 'experimentCTRLvsACLR3DOF_interp.mat'],...
                               'dataCTRL','dataACLR','labels');
    dataCTRL            = aux.dataCTRL;
    dataACLR            = aux.dataACLR;
    labels              = aux.labels;
    clear aux;
    downtime            = linspace(0,100,100);
    for i=1:length(numVar)
        idx     = find(isspace(varlabels{numVar(i)}));
        type    = varlabels{numVar(i)}(1:idx(1)-1);
        body    = varlabels{numVar(i)}(idx(1)+1:idx(2)-1);
        element = str2double(varlabels{numVar(i)}(idx(2)+3:end));
        switch(type)
            case 'F_muscle'
                bodyIdx           = find(strcmp(labels.muscleLabels,body));
                plotlabels.ylabel = 'Force (BW)';
                plotlabels.title   = ['F_{' body '} PC' num2str(element)];
            case 'phi_joint'
                bodyIdx           = find(strcmp(labels.angleLabels,body));
                plotlabels.ylabel = 'Angle (deg)';
                plotlabels.title   = ['\phi_{' body '} PC' num2str(element)];
            case 'M_joint'
                bodyIdx           = find(strcmp(labels.angleLabels,body));
                plotlabels.ylabel = 'Moment (BW-HT)';
                plotlabels.title   = ['M_{' body '} PC' num2str(element)];
        end
        rawdata            = [squeeze(dataCTRL.(type)(bodyIdx,:,:))';
                              squeeze(dataACLR.(type)(bodyIdx,:,:))'];
        codedata           = zeros(size(rawdata,1),length(downtime));
        for k=1:size(rawdata,1)
            % Downsampling the data
            codedata(k,:)  = interp1(dataCTRL.time,rawdata(k,:),downtime);
        end
        meanlo             = mean(rawdata(pcabotidx.(type).(body)(:,element),:),1);    % This is not low, but high
        stdlo              =  std(rawdata(pcabotidx.(type).(body)(:,element),:),[],1);
        meanhi             = mean(rawdata(pcatopidx.(type).(body)(:,element),:),1);    % This is not high, but low (sorry my mistake :P)
        stdhi              =  std(rawdata(pcatopidx.(type).(body)(:,element),:),[],1);
        normPC             = weightedcorrs([trainingdata(:,numVar(i)) codedata],weights);
        normPC             = normPC(2:end,1).^2;
        [~,maxEffect]    = max(normPC);
        
        plotlabels.legends = {'25%Q','75%Q'};
        subplot(2,2,1);
        plotGroups(dataCTRL.time,meanlo,stdlo,meanhi,stdhi,plotlabels,i==1); % Blue is HIGH
        ln1 = line([downtime(maxEffect) downtime(maxEffect)],[-100 100]);
        set(ln1,'LineWidth',1,'LineStyle','--','Color',[0 0 0]);
        subplot(2,2,3);
        [AX,H1,H2]         = plotyy(downtime,pcacoeffs.(type).(body)(:,element),downtime,100.*normPC);
        set(get(AX(1),'Ylabel'),'String','PC Coefficient (-)','Color',[0 0 0]);
        set(AX(1),'YColor',[0 0 0]);
        set(get(AX(2),'Ylabel'),'String','% of Explained Variance (--)','Color',[0 0 0]);
        set(AX(2),'YColor',[0 0 0]);
        xlabel('Landing Phase (%)');
        set(H1,'LineWidth',2,'LineStyle','-','Color',[0 0 0]);
        set(H2,'LineWidth',2,'LineStyle','--','Color',[0.4 0.4 0.4]);
        ln2 = line([downtime(maxEffect) downtime(maxEffect)],[-100 100]);
        set(ln2,'LineWidth',1,'LineStyle','--','Color',[0 0 0]);
        saveas(gcf,['C:/Users/mariom/Documents/MATLAB/GaitProject/pca_graphs/' plotlabels.title '.fig'],'fig');
        saveas(gcf,['C:/Users/mariom/Documents/MATLAB/GaitProject/pca_graphs/' plotlabels.title '.png'],'png');
    end
%end

if isempty(who('-file',[outputDatabase 'PCAanalysis3DOF_inter.mat'],'statsFSEL'))
    numIter                 = 1000;
    statsFSEL               = cell(numIter,1);
    normaxis                = linspace(0,1,numIter);
    fprFSEL                 = zeros(numIter);
    tprFSEL                 = zeros(numIter);
    oobVarSignif            = zeros(10,numIter);
    AUC                     = zeros(1,numIter);
    weights                 = (~isACLR).*sum(isACLR) + (isACLR).*sum(~isACLR);
    weights                 = weights./sum(weights);
    for i=1:numIter
        disp(['Training RandomForest No. ' num2str(i)]);
        treeFSEL            = TreeBagger(100,trainingdata(:,numVar),isACLR,...
                                         'oobvarimp','on','weights',weights);
        oobVarSignif(:,i)   = treeFSEL.OOBPermutedVarDeltaError;
        [YfitFSEL,SfitFSEL] = oobPredict(treeFSEL);
        [fpr,tpr,~,AUC(i)]  = perfcurve(treeFSEL.Y,SfitFSEL(:,2),'1');
        fprFSEL(:,i)        = interp1(linspace(0,1,length(fpr)),fpr,normaxis);
        tprFSEL(:,i)        = interp1(linspace(0,1,length(tpr)),tpr,normaxis);
        statsFSEL{i}        = classperf(treeFSEL.Y,YfitFSEL);
    end
    save([outputDatabase 'PCAanalysis3DOF_inter.mat'],'-append','statsFSEL','tprFSEL',...
         'fprFSEL','oobVarSignif','AUC');
else
    aux        = load([outputDatabase 'PCAanalysis3DOF_inter.mat'],'statsFSEL');
    statsFSEL  = aux.statsFSEL;
    clear aux;
end

hcorrela = corr(trainingdata,trainingdata(:,numVar));
inc      = 1;
for i=1:length(numVar)
    hcorrela(numVar(i),i) = NaN;
    idx                   = abs(hcorrela(:,i))>0.5;
    for j=1:length(idx)
        if idx(j)
            corrtable{inc,1}  = varlabels{numVar(i)};
            corrtable{inc,2}  = varlabels{j};
            corrtable{inc,3}  = hcorrela(j,i);
            corrtable{inc,4}  = mean(trainingdata(~isACLR,j));
            corrtable{inc,5}  =  std(trainingdata(~isACLR,j));
            corrtable{inc,6}  = mean(trainingdata( isACLR,j));
            corrtable{inc,7}  =  std(trainingdata( isACLR,j));
            [corrtable{inc,8},corrtable{inc,9}]  = ttest2(trainingdata(~isACLR,j),trainingdata(isACLR,j));
            mesvals = mes(trainingdata(~isACLR,j),trainingdata(isACLR,j),'hedgesg');
            corrtable{inc,10} = mesvals.hedgesg;
            inc = inc+1;
        end
    end
end

[coeff,scores,latent,~,explained] = pca(trainingdata(:,numVar));
pointCTRL = scores(~isACLR,1:2);
pointACLR = scores( isACLR,1:2);
hullCTRL  = convhull(pointCTRL);
hullACLR = convhull(pointACLR);
%plot3(pointCTRL(:,1),pointCTRL(:,2),pointCTRL(:,3),'bo');
plot(pointCTRL(:,1),pointCTRL(:,2),'bo');
hold on;
%plot3(pointACLR(:,1),pointACLR(:,2),pointACLR(:,3),'ro');
plot(pointACLR(:,1),pointACLR(:,2),'ro');
%trisurf(hullCTRL,pointCTRL(:,1),pointCTRL(:,2),pointCTRL(:,3)); shading interp;
%trisurf(hullACLR,pointACLR(:,1),pointACLR(:,2),pointACLR(:,3));
plot(pointCTRL(hullCTRL,1),pointCTRL(hullCTRL,2),'b-');
plot(pointACLR(hullACLR,1),pointACLR(hullACLR,2),'r-');
hold off;


performance      = zeros(7,1000);
performance(1,:) = AUC;
numIter          = 1000;
for i=1:numIter
    performance(2,i) = statsFSEL{i}.Sensitivity;
    performance(3,i) = statsFSEL{i}.Specificity;
    performance(4,i) = statsFSEL{i}.CorrectRate;
    performance(5,i) = statsFSEL{i}.PositiveLikelihood;
    performance(6,i) = statsFSEL{i}.PositivePredictiveValue;
    performance(7,i) = statsFSEL{i}.NegativePredictiveValue;
end
tabfin = [mean(performance,2) std(performance,[],2)];


disp('ACK');
    
    
% ==============================================================================
% Subfunctions
% ==============================================================================
function plotGroups(time,meanCtrl,stdCtrl,meanTest,stdTest,labels,flag)
% 
upperCtrl = meanCtrl + stdCtrl;
lowerCtrl = meanCtrl - stdCtrl;
upperTest = meanTest + stdTest;
lowerTest = meanTest - stdTest;

if flag
    upperCtrl = upperCtrl.*(upperCtrl>=0);
    lowerCtrl = lowerCtrl.*(lowerCtrl>=0);
    upperTest = upperTest.*(upperTest>=0);
    lowerTest = lowerTest.*(lowerTest>=0);
end

fill([time fliplr(time)],[upperCtrl fliplr(lowerCtrl)],'b','EdgeColor','none');
hold on;
fill([time fliplr(time)],[upperTest fliplr(lowerTest)],'r','EdgeColor','none');
alpha(0.1);
plot(time,meanCtrl,'b',time,meanTest,'r','LineWidth',2);
hold off;
hmax = ceil( 100*max([max(upperCtrl) max(upperTest)]))/100;
hmin = floor(100*min([min(lowerCtrl) min(lowerTest)]))/100;
axis([0 100 hmin hmax]);
%legend(labels.legends);
title(labels.title); %xlabel('Landing Phase (%)');
ylabel(labels.ylabel);
% ------------------------------------------------------------------------------
function weights = calculateWeights(indexes)
% 
subjectIdx  = unique(indexes);
weights     = zeros(length(indexes),1);
for i=1:length(subjectIdx)
    idx     = indexes==subjectIdx(i);
    weights = weights + (idx./sum(idx));
end
% ------------------------------------------------------------------------------
function [means,percentiles] = parallelanalysis(nvars,ncases)

ndatsets  = 1000;  % number of ndatsets
percent   = 95;    % desired percentile
for nds = 1:ndatsets
    evals(:,nds) = eig(corrcoef(randn(ncases,nvars)));
end
evals       = sort(evals,1,'descend');
means       = (mean(evals,2));   % mean eigenvalues for each position.
evals       = sort(evals,2);     % sorting the eigenvalues for each position.
percentiles = (evals(:,round((percent*ndatsets)/100)));  % percentiles.
