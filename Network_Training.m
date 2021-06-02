% ***** NETWORK DESIGN AND TRAINING *****

% ***** import label data *****
% change filepath to the folder where the training labels are located
label_data = fileDatastore(fullfile(...
    ''), ...
    'ReadFcn',@load,'FileExtensions','.mat');
%preview(label_data)

% ***** Import example data *****
% change filepath to the folder where the training examples are located
input_data = fileDatastore(fullfile(...
    ''),...
    'ReadFcn',@load,'FileExtensions','.mat');
%preview(input_data)


% ***** Data Transofrmation (struct -> cell) *****
input_datat = transform(input_data,@(data) rearrange_datastore_feature(data));
label_datat = transform(label_data,@(data) rearrange_datastore_label(data));

% ***** Create combined datastore *****
data_raw = combine(input_datat, label_datat);
%preview(trainData)

% ***** Shuffle datastore *****
data_shuffled = shuffle(data_raw); 

% ***** Subdivide datastore into training, validation and test sets *****
% set percentages as 70% training, 25% validation, 5% test.

% 500 data examples now
% REMEMBER TO CHANGE THIS!
examples_tot = n_samples*noise_iterations;
train_number = examples_tot*0.7;
valid_number = examples_tot*0.25;
test_number = examples_tot*0.05;
% calculate subset indices
train_indices = 1:train_number;
valid_indices = (train_number + 1):(train_number + valid_number);
test_indices = (train_number + valid_number + 1):examples_tot;

% create sub datastores
train_data = subset(data_shuffled, train_indices);
valid_data = subset(data_shuffled, valid_indices);
test_data = subset(data_shuffled, test_indices);

% % ***** Save Datastores *****
% save('train_data.mat', 'train_data');
% save('valid_data.mat', 'valid_data');
% save('test_data.mat', 'test_data');

% ***** Network Architecture *****

% 4-D input data convnet  (19 freq bins)
layers = [
    image3dInputLayer([11 6 19 3],"Name","image3dinput")
    convolution3dLayer([3 3 3],32,"Name","conv3d_1","Padding",[1 1 1;1 1 1])
    reluLayer("Name","relu_1")
    batchNormalizationLayer("Name","batchnorm_1")
    convolution3dLayer([3 3 3],64,"Name","conv3d_2","Padding",[1 1 1;1 1 1],"Stride",[2 1 3])
    reluLayer("Name","relu_2")
    batchNormalizationLayer("Name","batchnorm_2")
    convolution3dLayer([2 2 2],128,"Name","conv3d_3")
    reluLayer("Name","relu_3")
    batchNormalizationLayer("Name","batchnorm_3")
    convolution3dLayer([2 2 2],256,"Name","conv3d_4","Stride",[2 2 2])
    reluLayer("Name","relu_4")
    batchNormalizationLayer("Name","batchnorm_4")
    fullyConnectedLayer(96,"Name","fc")
    sigmoidLayer("Name","sigmoid")
    regressionLayer("Name","regressionoutput")];
% %plot(layerGraph(layers));

% % 4D Grouped convolution - group size 2.
% lgraph = layerGraph();
% 
% tempLayers = image3dInputLayer([11 6 19 3],"Name","image3dinput");
% lgraph = addLayers(lgraph,tempLayers);
% 
% tempLayers = [
%     convolution3dLayer([3 3 3],32,"Name","conv3d_1_2","Padding",[1 1 1;1 1 1])
%     reluLayer("Name","relu_1_2")
%     batchNormalizationLayer("Name","batchnorm_1_2")
%     convolution3dLayer([3 3 3],64,"Name","conv3d_2_2","Padding",[1 1 1;1 1 1],"Stride",[2 1 3])
%     reluLayer("Name","relu_2_2")
%     batchNormalizationLayer("Name","batchnorm_2_2")
%     convolution3dLayer([2 2 2],128,"Name","conv3d_3_2")
%     reluLayer("Name","relu_3_2")
%     batchNormalizationLayer("Name","batchnorm_3_2")
%     convolution3dLayer([2 2 2],256,"Name","conv3d_4_2","Stride",[2 2 2])
%     reluLayer("Name","relu_4_2")
%     batchNormalizationLayer("Name","batchnorm_4_2")];
% lgraph = addLayers(lgraph,tempLayers);
% 
% tempLayers = [
%     convolution3dLayer([3 3 3],32,"Name","conv3d_1_1","Padding",[1 1 1;1 1 1])
%     reluLayer("Name","relu_1_1")
%     batchNormalizationLayer("Name","batchnorm_1_1")
%     convolution3dLayer([3 3 3],64,"Name","conv3d_2_1","Padding",[1 1 1;1 1 1],"Stride",[2 1 3])
%     reluLayer("Name","relu_2_1")
%     batchNormalizationLayer("Name","batchnorm_2_1")
%     convolution3dLayer([2 2 2],128,"Name","conv3d_3_1")
%     reluLayer("Name","relu_3_1")
%     batchNormalizationLayer("Name","batchnorm_3_1")
%     convolution3dLayer([2 2 2],256,"Name","conv3d_4_1","Stride",[2 2 2])
%     reluLayer("Name","relu_4_1")
%     batchNormalizationLayer("Name","batchnorm_4_1")];
% lgraph = addLayers(lgraph,tempLayers);
% 
% tempLayers = [
%     depthConcatenationLayer(2,"Name","depthcat")
%     fullyConnectedLayer(96,"Name","fc")
%     sigmoidLayer("Name","sigmoid")
%     regressionLayer("Name","regressionoutput")];
% lgraph = addLayers(lgraph,tempLayers);
% 
% % clean up helper variable
% clear tempLayers;
% lgraph = connectLayers(lgraph,"image3dinput","conv3d_1_2");
% lgraph = connectLayers(lgraph,"image3dinput","conv3d_1_1");
% lgraph = connectLayers(lgraph,"batchnorm_4_1","depthcat/in1");
% lgraph = connectLayers(lgraph,"batchnorm_4_2","depthcat/in2");
% layers = lgraph;

% 4D input data convnet 2
% layers = [
%     image3dInputLayer([11 6 19 3],"Name","image3dinput")
%     convolution3dLayer([3 3 3],64,"Name","conv3d_1","Padding",[1 1 1;1 1 1])
%     reluLayer("Name","relu_1")
%     batchNormalizationLayer("Name","batchnorm_1")
%     convolution3dLayer([3 3 3],128,"Name","conv3d_2","Padding",[1 1 1;1 1 1],"Stride",[2 1 3])
%     reluLayer("Name","relu_2")
%     batchNormalizationLayer("Name","batchnorm_2")
%     convolution3dLayer([2 2 2],256,"Name","conv3d_3")
%     reluLayer("Name","relu_3")
%     batchNormalizationLayer("Name","batchnorm_3")
%     convolution3dLayer([2 2 2],512,"Name","conv3d_4","Stride",[2 2 2])
%     reluLayer("Name","relu_4")
%     batchNormalizationLayer("Name","batchnorm_4")
%     fullyConnectedLayer(96,"Name","fc_1")
%     sigmoidLayer("Name","sigmoid")
%     regressionLayer("Name","regressionoutput")];
% 
% %plot(lgraph);

% % 4-D U-NET
% layers = [
%     image3dInputLayer([11 6 19 3],"Name","image3dinput")
%     convolution3dLayer([3 3 3],32,"Name","conv3d_1","Padding",[1 1 1;1 1 1])
%     reluLayer("Name","relu_1")
%     batchNormalizationLayer("Name","batchnorm_1")
%     convolution3dLayer([3 3 3],64,"Name","conv3d_2","Padding",[1 1 1;1 1 1],"Stride",[2 1 3])
%     reluLayer("Name","relu_2")
%     batchNormalizationLayer("Name","batchnorm_2")
%     convolution3dLayer([2 2 2],128,"Name","conv3d_3")
%     reluLayer("Name","relu_3")
%     batchNormalizationLayer("Name","batchnorm_3")
%     convolution3dLayer([2 2 2],256,"Name","conv3d_4","Stride",[2 2 2])
%     reluLayer("Name","relu_4")
%     batchNormalizationLayer("Name","batchnorm_4")
%     convolution3dLayer([1 1 1],64,"Name","conv3d_5","DilationFactor",[2 2 2],"Padding",[1 1 1;1 1 1])
%     leakyReluLayer(0.2,"Name","leakyrelu_1")
%     batchNormalizationLayer("Name","batchnorm_5")
%     convolution3dLayer([1 1 1],32,"Name","conv3d_6","DilationFactor",[2 2 2],"Padding",[1 1 1;1 1 1])
%     leakyReluLayer(0.2,"Name","leakyrelu_2")
%     batchNormalizationLayer("Name","batchnorm_6")
%     convolution3dLayer([1 2 2],16,"Name","conv3d_7","DilationFactor",[2 2 2],"Padding",[1 1 1;1 1 1])
%     leakyReluLayer(0.2,"Name","leakyrelu_3")
%     batchNormalizationLayer("Name","batchnorm_7")
%     fullyConnectedLayer(96,"Name","fc_1")
%     sigmoidLayer("Name","sigmoid")
%     regressionLayer("Name","regressionoutput")];
% % plot 
% %plot(lgraph)

% ***** Training Options *****

miniBatchSize = 128;
options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'MaxEpochs', 2000, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',2, ...
    'ValidationData', valid_data, ...
    'Shuffle','every-epoch', ...
    'Verbose',true, ...
    'Plots','training-progress', ...
    'LearnRateSchedule', 'none', ...
    'LearnRateDropPeriod', 10, ...
    'LearnRateDropFactor', 0.2, ...
    'L2Regularization', 0.0001, ...
    'InitialLearnRate', 0.001, ... 
'ValidationFrequency', 300);

% ***** Network Training *****

% train new network
net = trainNetwork(train_data, layers, options);


% % resume training using pretrained net
% insert network filepath
% net2 = coder.loadDeepLearningNetwork( ...
%     '');
% net = trainNetwork(train_data, net2.Layers, options);

% ***** Network Testing *****
% predict labels form test data
test_pred_labels = predict(net, test_data);
% read test data labels
test_read_labels = zeros(96,test_number);
for i = 1:test_number
    test_read = read(test_data);
    % test_read_labels is a matrix having the i-th label array on the i-th
    % column
    test_read_labels(:,i) = test_read{1,2};
end
% transpose as to match with labels_pred_test
test_read_labels = test_read_labels.'; 

% compute mean square error
% SE for each label
test_error = zeros(1, test_number);
for i = 1:test_number
    test_error(i) = sum(abs(test_pred_labels(i,:) - test_read_labels(i,:)).^2);
end
% MSE for all labels
test_MSE = sum(test_error, 'all')/test_number;

% ***** Save Network *****
save('Net_Final.mat', ...
    'net');


% ***** Fetaures and Labels transormation functions (struct -> cell) *****
function image = rearrange_datastore_feature(data)
image = data.FFT_Data_Noise;
image = {image};
end


function image = rearrange_datastore_label(data)
image = data.random_array;
image = {image};
end
