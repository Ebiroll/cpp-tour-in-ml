#include "NeuralNetwork.hpp"
#include <cmath>
#include <iostream>

NeuralNetwork::NeuralNetwork(const std::vector<int> &architecture,
                             double learningRate,
                             Activation activation)
    : mLearningRate(learningRate), mActivation(activation), mArchitecture(architecture)
{
    // Initialize neurons and errors
    for (size_t i = 0; i < mArchitecture.size(); ++i)
    {
        int layerSize = mArchitecture[i] + (i < mArchitecture.size() - 1 ? 1 : 0); // Add bias neuron
        mNeurons.emplace_back(Eigen::VectorXd::Zero(layerSize));
        mErrors.emplace_back(Eigen::VectorXd::Zero(layerSize));
        if (i < mArchitecture.size() - 1)
        {
            mNeurons[i](layerSize - 1) = 1.0; // Set bias neuron activation to 1
        }
    }
    // Initialize weights
    initializeWeights();
}

void NeuralNetwork::initializeWeights()
{
    for (size_t i = 0; i < mArchitecture.size() - 1; ++i)
    {
        int rows = mNeurons[i].size();
        int cols = mNeurons[i + 1].size() - (i + 1 < mArchitecture.size() - 1 ? 1 : 0); // bias
        Eigen::MatrixXd weightMatrix = Eigen::MatrixXd::Random(rows, cols);
        mWeights.push_back(weightMatrix);
    }
}

double NeuralNetwork::activation(double x)
{
    if (mActivation == SIGMOID)
    {
        return 1.0 / (1.0 + std::exp(-x));
    }
    else
    { // TANH
        return std::tanh(x);
    }
}

double NeuralNetwork::activationDerivative(double x)
{
    if (mActivation == SIGMOID)
    {
        double sigmoid = activation(x);
        return sigmoid * (1 - sigmoid);
    }
    else
    { // TANH
        double tanhVal = std::tanh(x);
        return 1 - tanhVal * tanhVal;
    }
}

void NeuralNetwork::forward(const Eigen::VectorXd &input)
{
    mNeurons[0].head(mArchitecture[0]) = input;

    for (size_t i = 1; i < mArchitecture.size(); ++i)
    {
        Eigen::VectorXd netInput = mWeights[i - 1].transpose() * mNeurons[i - 1];

        for (int j = 0; j < netInput.size(); ++j)
        {
            mNeurons[i](j) = activation(netInput(j));
        }

        // Set bias neuron if not the output layer
        if (i < mArchitecture.size() - 1)
        {
            mNeurons[i](mNeurons[i].size() - 1) = 1.0; // Set bias neuron
        }
    }
}

void NeuralNetwork::backward(const Eigen::VectorXd &target)
{
    mLastEpochLoss = 0.0;
    // Output layer errors
    Eigen::VectorXd netInput = mWeights.back().transpose() * mNeurons[mNeurons.size() - 2];
    for (int i = 0; i < mArchitecture.back(); ++i)
    {
        double output = mNeurons.back()(i);
        double error = output - target(i);
        mErrors.back()(i) = error * activationDerivative(netInput(i));
        // Calculate mean squared error for this output
        mLastEpochLoss += error * error;
    }
    mLastEpochLoss = mLastEpochLoss / mArchitecture.back();

    // Hidden layer errors
    for (int i = mArchitecture.size() - 2; i > 0; --i)
    {
        Eigen::VectorXd netInput = mWeights[i - 1].transpose() * mNeurons[i - 1];
        Eigen::VectorXd delta = mWeights[i] * mErrors[i + 1].head(mArchitecture[i + 1]);
        for (int j = 0; j < mArchitecture[i]; ++j)
        {
            mErrors[i](j) = delta(j) * activationDerivative(netInput(j));
        }
    }

    // Update weights
    for (size_t i = 0; i < mWeights.size(); ++i)
    {
        Eigen::MatrixXd deltaWeight = mNeurons[i] * mErrors[i + 1].head(mArchitecture[i + 1]).transpose();
        mWeights[i] -= mLearningRate * deltaWeight;
    }
}

double NeuralNetwork::last_epoch_loss() const
{
    return mLastEpochLoss;
}

void NeuralNetwork::visualize_weights() const
{
    for (size_t i = 0; i < mWeights.size(); ++i)
    {
        std::cout << "Layer " << i << " weights:\n"
                  << mWeights[i] << std::endl;
    }
}

void NeuralNetwork::train(const Eigen::VectorXd &input, const Eigen::VectorXd &target)
{
    forward(input);
    backward(target);
}

Eigen::VectorXd NeuralNetwork::predict(const Eigen::VectorXd &input)
{
    forward(input);
    return mNeurons.back().head(mArchitecture.back());
}