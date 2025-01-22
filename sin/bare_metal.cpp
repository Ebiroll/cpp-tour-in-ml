#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <functional>
#include <random>

using namespace std;

#define SAMPLES 200
#define PI 3.141592653589793238463
#define N 10
#define epsilon 0.005
#define epoch 600


inline double approx_function(double x) {
	return sin(x);
	// Normal distribution
	//return (exp(-(x-PI) * (x-PI) / 2) / sqrt(2 * PI));
}

// Automatic differentiation variable
class Variable {
public:
    double data;
    double grad;
    vector<shared_ptr<Variable>> parents;
    function<void()> backward_fn;

    Variable(double data) : data(data), grad(0) {}

    void backward() {
        grad = 1.0;
        _backward();
    }

    void _backward() {
        if(backward_fn) backward_fn();
        for(auto& parent : parents) {
            parent->_backward();
        }
    }

    void zero_grad() {
        grad = 0;
        for(auto& parent : parents) {
            parent->zero_grad();
        }
    }
};

// Operation overloads
shared_ptr<Variable> operator+(shared_ptr<Variable> a, shared_ptr<Variable> b) {
    auto result = make_shared<Variable>(a->data + b->data);
    result->parents = {a, b};
    result->backward_fn = [=]() {
        a->grad += result->grad;
        b->grad += result->grad;
    };
    return result;
}

shared_ptr<Variable> operator-(shared_ptr<Variable> a, shared_ptr<Variable> b) {
    auto result = make_shared<Variable>(a->data - b->data);
    result->parents = {a, b};
    result->backward_fn = [=]() {
        a->grad += result->grad;
        b->grad -= result->grad;
    };
    return result;
}

shared_ptr<Variable> operator*(shared_ptr<Variable> a, shared_ptr<Variable> b) {
    auto result = make_shared<Variable>(a->data * b->data);
    result->parents = {a, b};
    result->backward_fn = [=]() {
        a->grad += b->data * result->grad;
        b->grad += a->data * result->grad;
    };
    return result;
}

shared_ptr<Variable> sigmoid(shared_ptr<Variable> x) {
    auto result = make_shared<Variable>(1.0 / (1.0 + exp(-x->data)));
    result->parents = {x};
    result->backward_fn = [=]() {
        x->grad += result->data * (1 - result->data) * result->grad;
    };
    return result;
}

// Network parameters
vector<shared_ptr<Variable>> W(N);
vector<shared_ptr<Variable>> V(N);
vector<shared_ptr<Variable>> c(N);

shared_ptr<Variable> b = make_shared<Variable>(0.0);

// Neural network function
shared_ptr<Variable> f_theta(shared_ptr<Variable> x) {
    auto result = b;
    for(int i = 0; i < N; i++) {
        auto wx = W[i] * x;
        auto c_wx = c[i] + wx;
        auto sig = sigmoid(c_wx);
        result = result + V[i] * sig;
    }
    return result;
}


// Training function
void train(shared_ptr<Variable> x, shared_ptr<Variable> y) {
    auto pred = f_theta(x);
    auto error = pred - y;  // Now using overloaded operator-
    auto loss = error * error;
    
    loss->backward();
    
    // Update parameters
    for(int i = 0; i < N; i++) {
        W[i]->data -= epsilon * W[i]->grad;
        V[i]->data -= epsilon * V[i]->grad;
        c[i]->data -= epsilon * c[i]->grad;
    }
    b->data -= epsilon * b->grad;
   
    // Reset gradients
    for(int i = 0; i < N; i++) {
        W[i]->zero_grad();
        V[i]->zero_grad();
        c[i]->zero_grad();
    }

    b->zero_grad();
    x->zero_grad();
    y->zero_grad();
}

// Shuffle data for vector of pairs
template<typename T>
void shuffle_data(std::vector<std::pair<T, T>>& vec) {
    static std::random_device rd;
    static std::mt19937 g(rd());
    for (size_t i = vec.size() - 1; i > 0; --i) {
        std::uniform_int_distribution<size_t> distrib(0, i);
        size_t j = distrib(g);
        std::swap(vec[i], vec[j]);
    }
}


int main() {
    // Initialize parameters
    srand(time(NULL));
    for(int i = 0; i < N; i++) {
        W[i] = make_shared<Variable>(2.0 * rand()/RAND_MAX - 1.0);
        V[i] = make_shared<Variable>(2.0 * rand()/RAND_MAX - 1.0);
        c[i] = make_shared<Variable>(2.0 * rand()/RAND_MAX - 1.0);
    }

    // Create training data
    vector<pair<shared_ptr<Variable>, shared_ptr<Variable>>> startSet;
    for(int i = 0; i < SAMPLES; i++) {
        double x_val = i * 2* PI / SAMPLES;
        double y_val = approx_function(x_val) + 0.05 * (1.0 * rand()/RAND_MAX - 1.0);
        startSet.emplace_back(
            make_shared<Variable>(x_val),
            make_shared<Variable>(y_val)
        );
    }

	// Shuffle training data at start of each epoch
    shuffle_data(startSet);

    // Splits
    int TRAIN_SPLIT = static_cast<int>(0.8f * SAMPLES);
    int TEST_SPLIT = static_cast<int>(0.2f * SAMPLES) + TRAIN_SPLIT;

	vector<pair<shared_ptr<Variable>, shared_ptr<Variable>> > trainSet(startSet.begin(), startSet.begin() + TRAIN_SPLIT);
	vector<pair<shared_ptr<Variable>, shared_ptr<Variable>> > testSet(startSet.begin() + TRAIN_SPLIT, startSet.begin() + TEST_SPLIT);

    vector<double> loss_history;
    vector<double> validation_history;


    // Training loop
    for(int j = 0; j < epoch; j++) {
	    shuffle_data(trainSet);
        double total_loss = 0.0;
		double validation_loss = 0.0;

        
        for(auto& [x, y] : trainSet) {
            train(x, y);
            total_loss += pow(f_theta(x)->data - y->data, 2);
        }

		for(auto& [x, y] : testSet) {
			validation_loss += pow(f_theta(x)->data - y->data, 2);
		}

        
		validation_history.push_back(validation_loss / SAMPLES);
        loss_history.push_back(total_loss / SAMPLES);
        
        if(j % 100 == 0) {
            cout << "Epoch: " << j << " Loss: " << loss_history.back()  << " " << validation_history.back() << endl;
        }
    }

    // Generate plot data
    vector<double> x_plot, y_true, y_pred;
    for(int i = 0; i < 1000; i++) {
        double x_val = i * 2* PI / 1000;
        auto x_var = make_shared<Variable>(x_val);
        auto y_var = f_theta(x_var);
        x_plot.push_back(x_val);
        y_true.push_back(approx_function(x_val));
        y_pred.push_back(y_var->data);
    }

    // Plotting code (requires gnuplot)
    FILE *gp = popen("gnuplot -persist", "w");
    if (!gp) {
        cerr << "Error opening gnuplot!" << endl;
        return 1;
    }

    fprintf(gp, "set terminal wxt size 1000,500\n");
    fprintf(gp, "set multiplot layout 1,3\n");
    
    // Plot 1: Function approximation
    fprintf(gp, "set title 'Function Approximation'\n");
    fprintf(gp, "plot '-' with lines title 'True', '-' with lines title 'Predicted'\n");
    for(size_t i = 0; i < x_plot.size(); ++i) {
        fprintf(gp, "%f %f\n", x_plot[i], y_true[i]);
    }
    fprintf(gp, "e\n");
    for(size_t i = 0; i < x_plot.size(); ++i) {
        fprintf(gp, "%f %f\n", x_plot[i], y_pred[i]);
    }
    fprintf(gp, "e\n");
	// Plot training set
		fprintf(gp, "set title 'Training Set'\n");
		fprintf(gp, "plot '-' with points title 'Training Set'\n");
		for (size_t i = 0; i < trainSet.size(); ++i) {
			fprintf(gp, "%f %f\n", trainSet[i].first->data, trainSet[i].second->data);
		}	
	fprintf(gp, "e\n"); 

    // Plot 2: Training loss
    fprintf(gp, "set title 'Validation Loss'\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "plot '-' with lines title 'Loss'\n");
    for(size_t i = 0; i < validation_history.size(); ++i) {
        fprintf(gp, "%zu %f\n", i, loss_history[i]);
    }
    fprintf(gp, "e\n");

    fprintf(gp, "unset multiplot\n");
    fflush(gp);
    
    cout << "Training complete. Close plots to exit...\n";
    getchar();
    pclose(gp);

    return 0;
}