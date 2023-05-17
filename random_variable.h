#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <cstdlib>
#include <map>
#include <cmath>
#include <vector>

class RandomVariable
{
public:
    RandomVariable(size_t experiments, size_t n, double p);

    double get_expected_value();
    double get_sample_expected_value();
    double get_variance();
    double get_sample_variance();
    double get_sample_median();
    double get_range();

    std::vector<double> get_probabilities();
    std::vector<double> get_left_borders();
    const std::vector<size_t>& get_random_variables();
    void simulate_experiment();
    std::map<size_t, size_t> get_frequencies(const std::vector<size_t>& rvs);

private:
    std::vector<size_t> count_random_variables();

private:
    size_t experiments;
    size_t n;
    double p;
    std::vector<size_t> random_variables;
    std::vector<double> probabilities;
};

#endif // RANDOMVARIABLE_H
