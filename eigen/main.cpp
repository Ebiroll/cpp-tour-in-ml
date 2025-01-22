#include "NeuralNetwork.hpp"
#include <iostream>

int main()
{
    // Define network architecture: 2 inputs, 2 hidden neurons, 1 output
    std::vector<int> architecture = {2, 2, 1};
    NeuralNetwork nn(architecture, 0.1, NeuralNetwork::SIGMOID);

    // Training data for XOR problem
    std::vector<Eigen::VectorXd> inputs = {
        (Eigen::VectorXd(2) << 0, 0).finished(),
        (Eigen::VectorXd(2) << 0, 1).finished(),
        (Eigen::VectorXd(2) << 1, 0).finished(),
        (Eigen::VectorXd(2) << 1, 1).finished()};
    std::vector<Eigen::VectorXd> targets = {
        (Eigen::VectorXd(1) << 0).finished(),
        (Eigen::VectorXd(1) << 1).finished(),
        (Eigen::VectorXd(1) << 1).finished(),
        (Eigen::VectorXd(1) << 0).finished()};

    // Train the network
    for (int epoch = 0; epoch < 10000; ++epoch)
    {
        for (size_t i = 0; i < inputs.size(); ++i)
        {
            nn.train(inputs[i], targets[i]);
        }
    }

    // Test the network
    for (size_t i = 0; i < inputs.size(); ++i)
    {
        Eigen::VectorXd output = nn.predict(inputs[i]);
        // Round to nearest integer
        output = output.unaryExpr([](double elem)
                                  { return (elem > 0.5) ? 1.0 : 0.0; });
        std::cout << "Input: " << inputs[i].transpose() << " Output: " << output(0) << std::endl;
    }

    return 0;
}