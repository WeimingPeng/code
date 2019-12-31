function [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% 返回经过训练的分类器及其准确度。以下代码重新创建在 Classification Learner App 中训
% 练的分类模型。您可以使用该生成的代码基于新数据自动训练同一模型，或通过它了解如何以程序化方
% 式训练模型。
%
%  输入:
%      trainingData: 一个所含列数和数据类型与导入 App 中的相同的矩阵。
%
%  输出:
%      trainedClassifier: 一个包含训练的分类器的结构体。该结构体中具有各种关于所训练分
%       类器的信息的字段。
%
%      trainedClassifier.predictFcn: 一个对新数据进行预测的函数。
%
%      validationAccuracy: 一个包含准确度百分比的双精度值。在 App 中，"历史记录" 列
%       表显示每个模型的此总体准确度分数。
%
% 使用该代码基于新数据来训练模型。要重新训练分类器，请使用原始数据或新数据作为输入参数
% trainingData 从命令行调用该函数。
%
% 例如，要重新训练基于原始数据集 T 训练的分类器，请输入:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% 要使用返回的 "trainedClassifier" 对新数据 T2 进行预测，请使用
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 必须是仅包含用于训练的预测变量列的矩阵。有关详细信息，请输入:
%   trainedClassifier.HowToPredict

% 由 MATLAB 于 2019-11-05 16:25:34 自动生成


% 提取预测变量和响应
% 以下代码将数据处理为合适的形状以训练模型。
%
% 将输入转换为表
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_7;
isCategoricalPredictor = [false, false, false, false, false, false];

% 将 PCA 应用于预测变量矩阵。
% 仅对数值预测变量运行 PCA。PCA 不会对分类预测变量进行任何处理。
isCategoricalPredictorBeforePCA = isCategoricalPredictor;
numericPredictors = predictors(:, ~isCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% 在 PCA 中必须将 'inf' 值视为缺失数据。
numericPredictors(isinf(numericPredictors)) = NaN;
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors);
% 保留足够的成分来解释所需的方差量。
explainedVarianceToKeepAsFraction = 95/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = [array2table(pcaScores(:,1:numComponentsToKeep)), predictors(:, isCategoricalPredictor)];
isCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(isCategoricalPredictor))];

% 训练分类器
% 以下代码指定所有分类器选项并训练分类器。
classificationTree = fitctree(...
    predictors, ...
    response, ...
    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', 100, ...
    'Surrogate', 'off', ...
    'ClassNames', [1; 2]);

% 使用预测函数创建结果结构体
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
treePredictFcn = @(x) predict(classificationTree, x);
trainedClassifier.predictFcn = @(x) treePredictFcn(pcaTransformationFcn(predictorExtractionFcn(x)));

% 向结果结构体中添加字段
trainedClassifier.PCACenters = pcaCenters;
trainedClassifier.PCACoefficients = pcaCoefficients;
trainedClassifier.ClassificationTree = classificationTree;
trainedClassifier.About = '此结构体是从 Classification Learner R2019b 导出的训练模型。';
trainedClassifier.HowToPredict = sprintf('要对新预测变量列矩阵 X 进行预测，请使用: \n yfit = c.predictFcn(X) \n将 ''c'' 替换为作为此结构体的变量的名称，例如 ''trainedModel''。\n \nX 必须包含正好 6 个列，因为此模型是使用 6 个预测变量进行训练的。\nX 必须仅包含与训练数据具有完全相同的顺序和格式的\n预测变量列。不要包含响应列或未导入 App 的任何列。\n \n有关详细信息，请参阅 <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>。');

% 提取预测变量和响应
% 以下代码将数据处理为合适的形状以训练模型。
%
% 将输入转换为表
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_7;
isCategoricalPredictor = [false, false, false, false, false, false];

% 设置留出法验证
cvp = cvpartition(response, 'Holdout', 0.2);
trainingPredictors = predictors(cvp.training, :);
trainingResponse = response(cvp.training, :);
trainingIsCategoricalPredictor = isCategoricalPredictor;

% 将 PCA 应用于预测变量矩阵。
% 仅对数值预测变量运行 PCA。PCA 不会对分类预测变量进行任何处理。
isCategoricalPredictorBeforePCA = trainingIsCategoricalPredictor;
numericPredictors = trainingPredictors(:, ~trainingIsCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% 在 PCA 中必须将 'inf' 值视为缺失数据。
numericPredictors(isinf(numericPredictors)) = NaN;
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors);
% 保留足够的成分来解释所需的方差量。
explainedVarianceToKeepAsFraction = 95/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
trainingPredictors = [array2table(pcaScores(:,1:numComponentsToKeep)), trainingPredictors(:, trainingIsCategoricalPredictor)];
trainingIsCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(trainingIsCategoricalPredictor))];

% 训练分类器
% 以下代码指定所有分类器选项并训练分类器。
classificationTree = fitctree(...
    trainingPredictors, ...
    trainingResponse, ...
    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', 100, ...
    'Surrogate', 'off', ...
    'ClassNames', [1; 2]);

% 使用预测函数创建结果结构体
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
treePredictFcn = @(x) predict(classificationTree, x);
validationPredictFcn = @(x) treePredictFcn(pcaTransformationFcn(x));

% 向结果结构体中添加字段


% 计算验证预测
validationPredictors = predictors(cvp.test, :);
validationResponse = response(cvp.test, :);
[validationPredictions, validationScores] = validationPredictFcn(validationPredictors);

% 计算验证准确度
correctPredictions = (validationPredictions == validationResponse);
isMissing = isnan(validationResponse);
correctPredictions = correctPredictions(~isMissing);
validationAccuracy = sum(correctPredictions)/length(correctPredictions);
