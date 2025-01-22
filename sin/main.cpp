#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

#include "tensorflow/cc/client/client_session.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/cc/ops/resource_variable_ops.h"
#include "tensorflow/cc/framework/scope.h"
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/cc/framework/gradients.h"
#include "tensorflow/core/public/session_options.h"


// gnuplot-iostream for plotting
#include "gnuplot-iostream.h"

static const int SAMPLES = 1000;
static const double PI = 3.14159265358979323846;

// Shuffle helper
template<typename T>
void shuffle_in_unison(std::vector<T>& a, std::vector<T>& b) {
    static std::random_device rd;
    static std::mt19937 g(rd());
    for (size_t i = a.size() - 1; i > 0; --i) {
        std::uniform_int_distribution<size_t> distrib(0, i);
        size_t j = distrib(g);
        std::swap(a[i], a[j]);
        std::swap(b[i], b[j]);
    }
}

// Template plot functio

void plotDoubles(Gnuplot& gp, const std::vector<std::pair<double, double>> data, const std::string& title) {
    gp << "plot '-' with points title '" << title << "'\n";
    gp.send1d(data);
}

void plotXY(Gnuplot& gp, const std::vector<float> x, const std::vector<float> y, const std::string& title) {
    std::vector<std::pair<double, double>> data;
    for (int i = 0; i < x.size(); ++i) {
        data.push_back({static_cast<double>(x[i]), (double)y[i]});
    }
    plotDoubles(gp, data, title);
}

int main() {
    Gnuplot gp;

    // 1) Generate data
    std::vector<float> x_values(SAMPLES), y_values(SAMPLES);
    
    {
        std::mt19937 gen(1337);
        std::uniform_real_distribution<float> dist(0.0f, 2.0f * static_cast<float>(PI));
        for (int i = 0; i < SAMPLES; ++i) {
            float x = dist(gen);
            x_values[i] = x;
            y_values[i] = std::sin(x);
        }
    }

    shuffle_in_unison(x_values, y_values);


    std::vector<std::pair<double, double>> sin_data;
    for (int i = 0; i < SAMPLES; ++i) {
        sin_data.push_back({x_values[i], y_values[i]});
    }

    // plotDoubles(gp,sin_data, std::string("sin(x)"));



    // Splits
    int TRAIN_SPLIT = static_cast<int>(0.6f * SAMPLES);
    int TEST_SPLIT = static_cast<int>(0.2f * SAMPLES) + TRAIN_SPLIT;

    std::vector<float> x_train(x_values.begin(), x_values.begin() + TRAIN_SPLIT);
    std::vector<float> y_train(y_values.begin(), y_values.begin() + TRAIN_SPLIT);
    std::vector<float> x_test(x_values.begin() + TRAIN_SPLIT, x_values.begin() + TEST_SPLIT);
    std::vector<float> y_test(y_values.begin() + TRAIN_SPLIT, y_values.begin() + TEST_SPLIT);

    // Set x_ordered to y_train length and hten set data from 0,2pi
    std::vector<float> x_ordered(y_train.size());
    std::iota(x_ordered.begin(), x_ordered.end(), 0.0f);

    plotXY(gp, x_ordered, y_train, "Training data");

    using namespace tensorflow;
    using namespace tensorflow::ops;

    // 2) Build graph
    const int HIDDEN_SIZE = 16;
    Scope root = Scope::NewRootScope();

    // Placeholders
    auto X = Placeholder(root.WithOpName("X"), DT_FLOAT);
    auto Y = Placeholder(root.WithOpName("Y"), DT_FLOAT);

    // Variables using ResourceVariables
    auto W1_handle = VarHandleOp(root.WithOpName("W1_handle"), DT_FLOAT, TensorShape({1, HIDDEN_SIZE}));
    auto b1_handle = VarHandleOp(root.WithOpName("b1_handle"), DT_FLOAT, TensorShape({1, HIDDEN_SIZE}));
    auto W2_handle = VarHandleOp(root.WithOpName("W2_handle"), DT_FLOAT, TensorShape({HIDDEN_SIZE, HIDDEN_SIZE}));
    auto b2_handle = VarHandleOp(root.WithOpName("b2_handle"), DT_FLOAT, TensorShape({HIDDEN_SIZE, HIDDEN_SIZE}));
    auto W3_handle = VarHandleOp(root.WithOpName("W3_handle"), DT_FLOAT, TensorShape({HIDDEN_SIZE, 1}));
    auto b3_handle = VarHandleOp(root.WithOpName("b3_handle"), DT_FLOAT, TensorShape({1, 1}));

    // Initialize variables with random values
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> urand(-0.1f, 0.1f);

    auto create_random_tensor = [&](const std::vector<int64_t>& shape) {
        Tensor t(DT_FLOAT, TensorShape(shape));
        auto flat = t.flat<float>();
        for (int i = 0; i < flat.size(); ++i) {
            flat(i) = urand(rng);
        }
        return t;
    };

    auto W1_init = AssignVariableOp(root, W1_handle, create_random_tensor({1, HIDDEN_SIZE}));
    auto b1_init = AssignVariableOp(root, b1_handle, create_random_tensor({1, HIDDEN_SIZE}));
    auto W2_init = AssignVariableOp(root, W2_handle, create_random_tensor({HIDDEN_SIZE, HIDDEN_SIZE}));
    auto b2_init = AssignVariableOp(root, b2_handle, create_random_tensor({HIDDEN_SIZE, HIDDEN_SIZE}));
    auto W3_init = AssignVariableOp(root, W3_handle, create_random_tensor({HIDDEN_SIZE, 1}));
    auto b3_init = AssignVariableOp(root, b3_handle, create_random_tensor({1, 1}));

    // Read variables for forward pass
    auto W1 = ReadVariableOp(root.WithOpName("W1"), W1_handle, DT_FLOAT);
    auto b1 = ReadVariableOp(root.WithOpName("b1"), b1_handle, DT_FLOAT);
    auto W2 = ReadVariableOp(root.WithOpName("W2"), W2_handle, DT_FLOAT);
    auto b2 = ReadVariableOp(root.WithOpName("b2"), b2_handle, DT_FLOAT);
    auto W3 = ReadVariableOp(root.WithOpName("W3"), W3_handle, DT_FLOAT);
    auto b3 = ReadVariableOp(root.WithOpName("b3"), b3_handle, DT_FLOAT);

    // Forward pass
    auto hidden = Relu(root, Add(root, MatMul(root, X, W1), b1));
    auto output = Add(root, MatMul(root, hidden, W2), b2);

    // Loss
    //Tensor reduction_indices_tensor(DT_INT32, TensorShape({1}));
    //reduction_indices_tensor.vec<int32>()(0) = 0;
    //auto reduction_indices = Const(root, reduction_indices_tensor);
    auto diff = Sub(root, output, Y);
    auto sq = Square(root, diff);
    auto loss = Mean(root.WithOpName("loss"), sq, 0 /*reduction_indices*/);

    // Compute gradients
    std::vector<Output> grad_outputs;
    TF_CHECK_OK(AddSymbolicGradients(root, {loss}, {Output(W1), Output(b1), Output(W2), Output(b2),Output(W3),Output(b3)}, &grad_outputs));
    
    // Optimizer
    float learning_rate = 0.01f;
    auto apply_W1 = ApplyGradientDescent(root.WithOpName("apply_W1"), W1_handle, learning_rate, grad_outputs[0]);
    auto apply_b1 = ApplyGradientDescent(root.WithOpName("apply_b1"), b1_handle, learning_rate, grad_outputs[1]);
    auto apply_W2 = ApplyGradientDescent(root.WithOpName("apply_W2"), W2_handle, learning_rate, grad_outputs[2]);
    auto apply_b2 = ApplyGradientDescent(root.WithOpName("apply_b2"), b2_handle, learning_rate, grad_outputs[3]);
    auto apply_W3 = ApplyGradientDescent(root.WithOpName("apply_W3"), W3_handle, learning_rate, grad_outputs[4]);
    auto apply_b3 = ApplyGradientDescent(root.WithOpName("apply_b3"), b3_handle, learning_rate, grad_outputs[5]);

    // Create session and initialize variables
    ClientSession session(root);
    std::vector<Operation> init_ops = {W1_init, b1_init, W2_init, b2_init,W3_init,b3_init};
    TF_CHECK_OK(session.Run({}, {}, init_ops, nullptr));


    // Training loop
    const int EPOCHS = 500;
    const int BATCH_SIZE = 16;
    int train_size = static_cast<int>(x_train.size());

    for (int epoch = 0; epoch < EPOCHS; ++epoch) {
        float epoch_loss = 0.0f;
        int batches = 0;

        for (int start = 0; start < train_size; start += BATCH_SIZE) {
            int end = std::min(start + BATCH_SIZE, train_size);
            int cur_batch_size = end - start;

            // Prepare batch data
            Tensor x_batch(DT_FLOAT, TensorShape({cur_batch_size, 1}));
            Tensor y_batch(DT_FLOAT, TensorShape({cur_batch_size, 1}));
            
            auto x_flat = x_batch.flat<float>();
            auto y_flat = y_batch.flat<float>();
            
            for (int i = 0; i < cur_batch_size; ++i) {
                x_flat(i) = x_train[start + i];
                y_flat(i) = y_train[start + i];
            }

            // Run training step
            std::vector<std::pair<string, Tensor>> inputs = {{"X", x_batch}, {"Y", y_batch}};
            std::vector<Output> fetch_outputs = {loss};
            //std::vector<Operation> train_ops{apply_W1, apply_b1, apply_W2, apply_b2,apply_W3,apply_b3};
            std::vector<Operation> train_ops{{grad_outputs[0],grad_outputs[1],grad_outputs[2],grad_outputs[3],grad_outputs[4],grad_outputs[5]}};

            //train_ops.push_back(apply_W1);
            //train_ops.push_back(apply_b1);
            //train_ops.push_back(apply_W2);
            //train_ops.push_back(apply_b2);
            std::vector<Tensor> out_tensors;

            tensorflow::RunOptions run_options;
            //TF_CHECK_OK(session.Run(run_options, inputs, fetch_outputs, train_ops, &out_tensors, nullptr));

            // candidate: ‘tsl::Status tensorflow::ClientSession::Run(const tensorflow::RunOptions&, 
            //const FeedType&, 
            // const std::vector<tensorflow::Output>&, 
            // const std::vector<tensorflow::Operation>&, s
            // std::vector<tensorflow::Tensor>*, tensorflow::RunMetadata*) const’
            tensorflow::RunMetadata run_metadata;
            // ‘tsl::Status tensorflow::ClientSession::Run(const tensorflow::RunOptions&, 
            //const FeedType&, const std::vector<tensorflow::Output>&,
            //  const std::vector<tensorflow::Operation>&, std::vector<tensorflow::Tensor>*, tensorflow::RunMetadata*) const’
            TF_CHECK_OK(session.Run(run_options,fetch_outputs, inputs, train_ops, &out_tensors, &run_metadata));

            float batch_loss = out_tensors[0].scalar<float>()();
            epoch_loss += batch_loss;
            ++batches;
        }

        epoch_loss /= batches;
        if ((epoch + 1) % 50 == 0) {
            std::cout << "Epoch " << (epoch + 1) << " - Loss: " << epoch_loss << std::endl;
        }
    }

    // Evaluate on test set
    int test_size = static_cast<int>(x_test.size());
    Tensor x_test_t(DT_FLOAT, TensorShape({test_size, 1}));
    Tensor y_test_t(DT_FLOAT, TensorShape({test_size, 1}));
    
    auto x_test_flat = x_test_t.flat<float>();
    auto y_test_flat = y_test_t.flat<float>();
    
    for (int i = 0; i < test_size; ++i) {
        x_test_flat(i) = x_test[i];
        y_test_flat(i) = y_test[i];
    }

    std::vector<std::pair<string, Tensor>> test_inputs = {{"X", x_test_t}, {"Y", y_test_t}};
    std::vector<Output> test_fetch = {loss};
    std::vector<Tensor> test_outputs;
    
    // TF_CHECK_OK(session.Run(test_inputs, test_fetch, {}, &test_outputs));
    
    float test_loss = test_outputs[0].scalar<float>()();
    std::cout << "Test Loss = " << test_loss << std::endl;

    // Generate predictions for plotting
    std::vector<std::pair<double, double>> original_data;
    std::vector<std::pair<double, double>> predicted_data;
    
    int N_PLOT = 200;
    double step = 2.0 * PI / (N_PLOT - 1);
    
    Tensor x_plot(DT_FLOAT, TensorShape({N_PLOT, 1}));
    auto x_plot_flat = x_plot.flat<float>();
    
    for (int i = 0; i < N_PLOT; ++i) {
        double x = i * step;
        original_data.push_back({x, std::sin(x)});
        x_plot_flat(i) = static_cast<float>(x);
    }

    std::vector<std::pair<string, Tensor>> plot_inputs = {{"X", x_plot}};
    std::vector<Output> plot_fetch = {output};
    std::vector<Tensor> pred_outputs;
    
    tensorflow::RunOptions eval_options;
    // TF_CHECK_OK(session.Run(eval_options,plot_inputs, plot_fetch, {}, &pred_outputs));
    
    auto pred_flat = pred_outputs[0].flat<float>();
    for (int i = 0; i < N_PLOT; ++i) {
        predicted_data.push_back({original_data[i].first, pred_flat(i)});
    }

    // Plot results
    
    gp << "set title 'Sine function approximation'\n";
    gp << "plot '-' with points title 'Original sin(x)', "
          "'-' with lines title 'NN Prediction'\n";
    gp.send1d(original_data);
    gp.send1d(predicted_data);

    return 0;
}

#if 0
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

#include "tensorflow/cc/client/client_session.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/cc/framework/scope.h"
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/core/platform/env.h"

// gnuplot-iostream for plotting
#include "gnuplot-iostream.h"

static const int SAMPLES = 1000;
static const double PI = 3.14159265358979323846;

// Shuffle helper
template<typename T>
void shuffle_in_unison(std::vector<T>& a, std::vector<T>& b) {
    static std::random_device rd;
    static std::mt19937 g(rd());
    for (size_t i = a.size() - 1; i > 0; --i) {
        std::uniform_int_distribution<size_t> distrib(0, i);
        size_t j = distrib(g);
        std::swap(a[i], a[j]);
        std::swap(b[i], b[j]);
    }
}

int main() {
    // 1) Generate data
    std::vector<float> x_values(SAMPLES), y_values(SAMPLES);
    
    {
        std::mt19937 gen(1337);
        std::uniform_real_distribution<float> dist(0.0f, 2.0f * static_cast<float>(PI));
        for (int i = 0; i < SAMPLES; ++i) {
            float x = dist(gen);
            x_values[i] = x;
            y_values[i] = std::sin(x);
        }
    }
    shuffle_in_unison(x_values, y_values);

    // Splits
    int TRAIN_SPLIT = static_cast<int>(0.6f * SAMPLES);
    int TEST_SPLIT = static_cast<int>(0.2f * SAMPLES) + TRAIN_SPLIT;

    std::vector<float> x_train(x_values.begin(), x_values.begin() + TRAIN_SPLIT);
    std::vector<float> y_train(y_values.begin(), y_values.begin() + TRAIN_SPLIT);
    std::vector<float> x_test(x_values.begin() + TRAIN_SPLIT, x_values.begin() + TEST_SPLIT);
    std::vector<float> y_test(y_values.begin() + TRAIN_SPLIT, y_values.begin() + TEST_SPLIT);

    using namespace tensorflow;
    using namespace tensorflow::ops;

    // 2) Build graph
    const int HIDDEN_SIZE = 16;
    
    Scope root = Scope::NewRootScope();
    
    // Placeholders
    auto X = Placeholder(root.WithOpName("X"), DT_FLOAT);
    auto Y = Placeholder(root.WithOpName("Y"), DT_FLOAT);


// Layer 1: 16 neurons with ReLU
auto W1 = VarHandleOp(root, DT_FLOAT, TensorShape({1, HIDDEN_SIZE}));
auto b1 = VarHandleOp(root, DT_FLOAT, TensorShape({HIDDEN_SIZE}));
auto layer1 = Relu(root, Add(root, MatMul(root, X, ReadVariableOp(root, W1, DT_FLOAT)), ReadVariableOp(root, b1, DT_FLOAT)));

// Layer 2: 16 neurons with ReLU
auto W2 = VarHandleOp(root, DT_FLOAT, TensorShape({HIDDEN_SIZE, HIDDEN_SIZE}));
auto b2 = VarHandleOp(root, DT_FLOAT, TensorShape({HIDDEN_SIZE}));
auto layer2 = Relu(root, Add(root, MatMul(root, layer1, ReadVariableOp(root, W2, DT_FLOAT)), ReadVariableOp(root, b2, DT_FLOAT)));

// Output layer
auto W3 = VarHandleOp(root, DT_FLOAT, TensorShape({HIDDEN_SIZE, 1}));
auto b3 = VarHandleOp(root, DT_FLOAT, TensorShape({1}));
auto output = Add(root, MatMul(root, layer2, ReadVariableOp(root, W3, DT_FLOAT)), ReadVariableOp(root, b3, DT_FLOAT));


    

    // Initialize variables
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> urand(-0.1f, 0.1f);

    auto create_random_tensor = [&](const std::vector<int64_t>& shape) {
        Tensor t(DT_FLOAT, TensorShape(shape));
        auto flat = t.flat<float>();
        for (int i = 0; i < flat.size(); ++i) {
            flat(i) = urand(rng);
        }
        return t;
    };

    // Create initialization ops with proper initialization
    auto W1_init = Assign(root, W1, create_random_tensor({1, HIDDEN_SIZE}));
    auto b1_init = Assign(root, b1, create_random_tensor({1, HIDDEN_SIZE}));
    auto W2_init = Assign(root, W2, create_random_tensor({HIDDEN_SIZE, 1}));
    auto b2_init = Assign(root, b2, create_random_tensor({1, 1}));

    // Forward pass
    auto hidden = Relu(root, Add(root, MatMul(root, X, W1), b1));
    auto output = Add(root, MatMul(root, hidden, W2), b2);
    
    // Loss - fixed Mean operation with reduction dimensions
    auto diff = Sub(root, output, Y);
    auto sq = Square(root, diff);
    std::vector<int> reduction_dims = {0};  // Reduce along batch dimension
    // OLAS auto loss = Mean(root.WithOpName("loss"), sq, reduction_dims);

    // Create a session and initialize variables
    ClientSession session(root);
    
    // Fixed initialization with proper vector
    std::vector<Output> init_ops = {W1_init, b1_init, W2_init, b2_init};
    std::vector<Tensor> outputs;
    //TF_CHECK_OK(session.Run({}, {}, init_ops, &outputs));

    // Training loop
    std::vector<std::pair<double, double>> epoch_loss_data;

    const int EPOCHS = 500;
    const int BATCH_SIZE = 16;
    
    int train_size = static_cast<int>(x_train.size());

    for (int epoch = 0; epoch < EPOCHS; ++epoch) {
        float epoch_loss = 0.0f;
        int batches = 0;

        for (int start = 0; start < train_size; start += BATCH_SIZE) {
            int end = std::min(start + BATCH_SIZE, train_size);
            int cur_batch_size = end - start;

            // Prepare batch data
            Tensor x_batch(DT_FLOAT, TensorShape({cur_batch_size, 1}));
            Tensor y_batch(DT_FLOAT, TensorShape({cur_batch_size, 1}));
            
            auto x_flat = x_batch.flat<float>();
            auto y_flat = y_batch.flat<float>();
            
            for (int i = 0; i < cur_batch_size; ++i) {
                x_flat(i) = x_train[start + i];
                y_flat(i) = y_train[start + i];
            }

            // Fixed Run call with proper types
            std::vector<std::pair<string, Tensor>> inputs = {
                {X.node()->name(), x_batch},
                {Y.node()->name(), y_batch}
            };
            //std::vector<Output> fetch_outputs = {loss, output};
            //std::vector<Tensor> out_tensors;
            
            //TF_CHECK_OK(session.Run(inputs, fetch_outputs, &out_tensors));
            
            //float batch_loss = out_tensors[0].scalar<float>()();
            //epoch_loss += batch_loss;
            //++batches;
        }

        epoch_loss /= batches;
        if ((epoch + 1) % 50 == 0) {
            epoch_loss_data.push_back({epoch + 1, epoch_loss});
            std::cout << "Epoch " << (epoch + 1) << " - Loss: " << epoch_loss << std::endl;
        }
    }
    Gnuplot learning_plot;
    learning_plot << "set title 'Loss vs Epoch'\n";
    gp.send1d(epoch_loss_data);


    // Evaluate on test set
    int test_size = static_cast<int>(x_test.size());
    Tensor x_test_t(DT_FLOAT, TensorShape({test_size, 1}));
    Tensor y_test_t(DT_FLOAT, TensorShape({test_size, 1}));
    
    auto x_test_flat = x_test_t.flat<float>();
    auto y_test_flat = y_test_t.flat<float>();
    
    for (int i = 0; i < test_size; ++i) {
        x_test_flat(i) = x_test[i];
        y_test_flat(i) = y_test[i];
    }

    std::vector<std::pair<string, Tensor>> test_inputs = {
        {X.node()->name(), x_test_t},
        {Y.node()->name(), y_test_t}
    };
    //std::vector<Output> test_fetch = {loss};
    std::vector<Tensor> test_outputs;
    
    //TF_CHECK_OK(session.Run(test_inputs, test_fetch, &test_outputs));
    
    float test_loss = test_outputs[0].scalar<float>()();
    std::cout << "Test Loss = " << test_loss << std::endl;

    // Generate predictions for plotting
    std::vector<std::pair<double, double>> original_data;
    std::vector<std::pair<double, double>> predicted_data;
    
    int N_PLOT = 200;
    double step = 2.0 * PI / (N_PLOT - 1);
    
    Tensor x_plot(DT_FLOAT, TensorShape({N_PLOT, 1}));
    auto x_plot_flat = x_plot.flat<float>();
    
    for (int i = 0; i < N_PLOT; ++i) {
        double x = i * step;
        original_data.push_back({x, std::sin(x)});
        x_plot_flat(i) = static_cast<float>(x);
    }

    std::vector<std::pair<string, Tensor>> plot_inputs = {{X.node()->name(), x_plot}};
    std::vector<Output> plot_fetch = {output};
    std::vector<Tensor> pred_outputs;
    
    //TF_CHECK_OK(session.Run(plot_inputs, plot_fetch, &pred_outputs));
    
    auto pred_flat = pred_outputs[0].flat<float>();
    for (int i = 0; i < N_PLOT; ++i) {
        predicted_data.push_back({original_data[i].first, pred_flat(i)});
    }

    // Plot results
    Gnuplot gp;
    gp << "set title 'Sine function approximation'\n";
    gp << "plot '-' with points title 'Original sin(x)', "
          "'-' with lines title 'NN Prediction'\n";
    gp.send1d(original_data);
    gp.send1d(predicted_data);

    return 0;
}
#endif