#ifndef NEURAL_NETWORK_HPP
#define NEURAL_NETWORK_HPP

#include <Eigen/Dense>
#include <vector>

class NeuralNetwork {
public:
    enum Activation { TANH, SIGMOID };

    NeuralNetwork(const std::vector<int>& architecture,
                  double learningRate = 0.01,
                  Activation activation = TANH);

    void train(const Eigen::VectorXd& input, const Eigen::VectorXd& target);
    Eigen::VectorXd predict(const Eigen::VectorXd& input);

private:
    double mLearningRate;
    Activation mActivation;
    std::vector<int> mArchitecture;
    std::vector<Eigen::VectorXd> mNeurons;
    std::vector<Eigen::VectorXd> mErrors;
    std::vector<Eigen::MatrixXd> mWeights;

    void initializeWeights();
    void forward(const Eigen::VectorXd& input);
    void backward(const Eigen::VectorXd& target);
    double activation(double x);
    double activationDerivative(double x);
};

#endif // NEURAL_NETWORK_HPP