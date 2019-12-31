function [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% ���ؾ���ѵ���ķ���������׼ȷ�ȡ����´������´����� Classification Learner App ��ѵ
% ���ķ���ģ�͡�������ʹ�ø����ɵĴ�������������Զ�ѵ��ͬһģ�ͣ���ͨ�����˽�����Գ��򻯷�
% ʽѵ��ģ�͡�
%
%  ����:
%      trainingData: һ���������������������뵼�� App �е���ͬ�ľ���
%
%  ���:
%      trainedClassifier: һ������ѵ���ķ������Ľṹ�塣�ýṹ���о��и��ֹ�����ѵ����
%       ��������Ϣ���ֶΡ�
%
%      trainedClassifier.predictFcn: һ���������ݽ���Ԥ��ĺ�����
%
%      validationAccuracy: һ������׼ȷ�Ȱٷֱȵ�˫����ֵ���� App �У�"��ʷ��¼" ��
%       ����ʾÿ��ģ�͵Ĵ�����׼ȷ�ȷ�����
%
% ʹ�øô��������������ѵ��ģ�͡�Ҫ����ѵ������������ʹ��ԭʼ���ݻ���������Ϊ�������
% trainingData �������е��øú�����
%
% ���磬Ҫ����ѵ������ԭʼ���ݼ� T ѵ���ķ�������������:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% Ҫʹ�÷��ص� "trainedClassifier" �������� T2 ����Ԥ�⣬��ʹ��
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 �����ǽ���������ѵ����Ԥ������еľ����й���ϸ��Ϣ��������:
%   trainedClassifier.HowToPredict

% �� MATLAB �� 2019-11-05 16:25:34 �Զ�����


% ��ȡԤ���������Ӧ
% ���´��뽫���ݴ���Ϊ���ʵ���״��ѵ��ģ�͡�
%
% ������ת��Ϊ��
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_7;
isCategoricalPredictor = [false, false, false, false, false, false];

% �� PCA Ӧ����Ԥ���������
% ������ֵԤ��������� PCA��PCA ����Է���Ԥ����������κδ���
isCategoricalPredictorBeforePCA = isCategoricalPredictor;
numericPredictors = predictors(:, ~isCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% �� PCA �б��뽫 'inf' ֵ��Ϊȱʧ���ݡ�
numericPredictors(isinf(numericPredictors)) = NaN;
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors);
% �����㹻�ĳɷ�����������ķ�������
explainedVarianceToKeepAsFraction = 95/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = [array2table(pcaScores(:,1:numComponentsToKeep)), predictors(:, isCategoricalPredictor)];
isCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(isCategoricalPredictor))];

% ѵ��������
% ���´���ָ�����з�����ѡ�ѵ����������
classificationTree = fitctree(...
    predictors, ...
    response, ...
    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', 100, ...
    'Surrogate', 'off', ...
    'ClassNames', [1; 2]);

% ʹ��Ԥ�⺯����������ṹ��
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
treePredictFcn = @(x) predict(classificationTree, x);
trainedClassifier.predictFcn = @(x) treePredictFcn(pcaTransformationFcn(predictorExtractionFcn(x)));

% �����ṹ��������ֶ�
trainedClassifier.PCACenters = pcaCenters;
trainedClassifier.PCACoefficients = pcaCoefficients;
trainedClassifier.ClassificationTree = classificationTree;
trainedClassifier.About = '�˽ṹ���Ǵ� Classification Learner R2019b ������ѵ��ģ�͡�';
trainedClassifier.HowToPredict = sprintf('Ҫ����Ԥ������о��� X ����Ԥ�⣬��ʹ��: \n yfit = c.predictFcn(X) \n�� ''c'' �滻Ϊ��Ϊ�˽ṹ��ı��������ƣ����� ''trainedModel''��\n \nX ����������� 6 ���У���Ϊ��ģ����ʹ�� 6 ��Ԥ���������ѵ���ġ�\nX �����������ѵ�����ݾ�����ȫ��ͬ��˳��͸�ʽ��\nԤ������С���Ҫ������Ӧ�л�δ���� App ���κ��С�\n \n�й���ϸ��Ϣ������� <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>��');

% ��ȡԤ���������Ӧ
% ���´��뽫���ݴ���Ϊ���ʵ���״��ѵ��ģ�͡�
%
% ������ת��Ϊ��
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_7;
isCategoricalPredictor = [false, false, false, false, false, false];

% ������������֤
cvp = cvpartition(response, 'Holdout', 0.2);
trainingPredictors = predictors(cvp.training, :);
trainingResponse = response(cvp.training, :);
trainingIsCategoricalPredictor = isCategoricalPredictor;

% �� PCA Ӧ����Ԥ���������
% ������ֵԤ��������� PCA��PCA ����Է���Ԥ����������κδ���
isCategoricalPredictorBeforePCA = trainingIsCategoricalPredictor;
numericPredictors = trainingPredictors(:, ~trainingIsCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% �� PCA �б��뽫 'inf' ֵ��Ϊȱʧ���ݡ�
numericPredictors(isinf(numericPredictors)) = NaN;
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors);
% �����㹻�ĳɷ�����������ķ�������
explainedVarianceToKeepAsFraction = 95/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
trainingPredictors = [array2table(pcaScores(:,1:numComponentsToKeep)), trainingPredictors(:, trainingIsCategoricalPredictor)];
trainingIsCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(trainingIsCategoricalPredictor))];

% ѵ��������
% ���´���ָ�����з�����ѡ�ѵ����������
classificationTree = fitctree(...
    trainingPredictors, ...
    trainingResponse, ...
    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', 100, ...
    'Surrogate', 'off', ...
    'ClassNames', [1; 2]);

% ʹ��Ԥ�⺯����������ṹ��
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
treePredictFcn = @(x) predict(classificationTree, x);
validationPredictFcn = @(x) treePredictFcn(pcaTransformationFcn(x));

% �����ṹ��������ֶ�


% ������֤Ԥ��
validationPredictors = predictors(cvp.test, :);
validationResponse = response(cvp.test, :);
[validationPredictions, validationScores] = validationPredictFcn(validationPredictors);

% ������֤׼ȷ��
correctPredictions = (validationPredictions == validationResponse);
isMissing = isnan(validationResponse);
correctPredictions = correctPredictions(~isMissing);
validationAccuracy = sum(correctPredictions)/length(correctPredictions);
