/*
    Simple neural network using Eigen math library
*/

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

double Sigmoid(const double z)
{
    return 1.0 / (1.0 + exp(-z));
}

double SigmoidDerivative(const double z)
{
    return z * (1.0 - z);
}

class NeuralNetwork
{
public:
    NeuralNetwork(const Eigen::MatrixXd &x, const Eigen::VectorXd &y);
    void Feedforward();
    void BackPropagation(double epsilon);
    Eigen::VectorXd predicted_output;
    void visualize_weights() const;

private:
    double GetLossValue() const;
    // Make int mutable to avoid const_cast
    mutable uint32_t loss_call_times_ = 0;
    uint32_t num_samples_;
    uint32_t num_features_;
    Eigen::MatrixXd input_;
    Eigen::VectorXd real_output_;
    const uint32_t neuron_count_1_ = 4;
    const uint32_t neuron_count_2_ = 1;
    Eigen::MatrixXd weights_1_;
    Eigen::MatrixXd weights_2_;
    Eigen::MatrixXd output_1_;
};

NeuralNetwork::NeuralNetwork(const Eigen::MatrixXd &x, const Eigen::VectorXd &y)
{
    input_ = x;
    real_output_ = y;
    num_samples_ = (uint32_t)input_.rows();
    num_features_ = (uint32_t)input_.cols();
    assert(num_samples_ == real_output_.rows());
    predicted_output = Eigen::VectorXd::Zero(num_samples_);
    weights_1_ = Eigen::MatrixXd::Random(num_features_, neuron_count_1_);
    weights_2_ = Eigen::MatrixXd::Random(neuron_count_1_, neuron_count_2_);
    std::cout << "layer 1 weights:\n " << weights_1_ << std::endl;
    std::cout << "layer 2 weights:\n " << weights_2_ << std::endl;
}

void NeuralNetwork::Feedforward()
{
    // Layer 1
    Eigen::MatrixXd z_1 = input_ * weights_1_;
    // Squishify values in interval 0-1 using Sigmoid function
    output_1_ = z_1.unaryExpr(std::ref(Sigmoid));

    // Layer 2
    Eigen::VectorXd z_2 = output_1_ * weights_2_;
    predicted_output = z_2.unaryExpr(std::ref(Sigmoid));
}

double NeuralNetwork::GetLossValue() const
{
    Eigen::VectorXd loss_vector = real_output_ - predicted_output;
    loss_vector = loss_vector.array().square();
    double loss_value = loss_vector.sum();
    if (loss_call_times_++ % 1000 == 0)
    {
        std::cout << "Loss value: " << loss_value << std::endl;
    }
    return loss_value;
}

void NeuralNetwork::visualize_weights() const
{
    std::cout << "layer 1 weights:\n " << weights_1_ << std::endl;
    std::cout << "layer 2 weights:\n " << weights_2_ << std::endl;
}
void NeuralNetwork::BackPropagation(double epsilon)
{
    // Calculate the loss and gradients
    GetLossValue();
    Eigen::VectorXd delta_output = real_output_ - predicted_output;
    Eigen::VectorXd predicted_output_der = predicted_output.unaryExpr(std::ref(SigmoidDerivative));

    // Calculate weight updates for the output layer
    Eigen::VectorXd delta_w_2 = output_1_.transpose() * (2 * delta_output.cwiseProduct(predicted_output_der));

    // Calculate weight updates for the hidden layer
    Eigen::MatrixXd delta_w_1 = input_.transpose() * (2 * delta_output.cwiseProduct(predicted_output_der) * weights_2_.transpose()).cwiseProduct(output_1_.unaryExpr(std::ref(SigmoidDerivative)));

    // Update weights with learning rate (epsilon)
    weights_1_ += epsilon * delta_w_1;
    weights_2_ += epsilon * delta_w_2;
}

int main()
{
    const uint32_t num_samples = 4;
    const uint32_t num_features = 2;
    Eigen::MatrixXd input(num_samples, num_features);
    input << 0.0, 0.0,
        0.0, 1.0,
        1.0, 0.0,
        1.0, 1.0;
    Eigen::VectorXd output(num_samples);
    output << 0, 1, 1, 0;

    NeuralNetwork net(input, output);

    double epsilon = 0.05;
    for (int i = 0; i < 10000; i++)
    {
        net.Feedforward();
        net.BackPropagation(epsilon);
    }
    // Use eigen rounding to get binary output
    net.predicted_output = net.predicted_output.unaryExpr([](double x)
                                                          { return std::round(x); });

    std::cout << "Predicted output:\n"
              << net.predicted_output << std::endl;

    // Input dummy before pinting weights
    std::cout << "Press any key to visualize weights" << std::endl;
    std::cin.get();
    net.visualize_weights();

    return 0;
}
