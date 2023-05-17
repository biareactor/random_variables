#include "random_variable.h"
#include <ctime>
#include <iostream>
#include <algorithm>

RandomVariable::RandomVariable(size_t _experiments, size_t _n, double _p)
    : experiments(_experiments), n(_n), p(_p)
{
    probabilities = get_probabilities();
    for (const auto& p : probabilities)
        std::cout << p << ' ';
    std::cout << '\n';
    //random_variables = count_random_variables();
    simulate_experiment();
    std::sort(random_variables.begin(), random_variables.end());
}

static std::vector<std::vector<size_t>> get_pascal_triangle(size_t n)
{
    std::vector<std::vector<size_t>> B(n+1, std::vector<size_t>(n+1));      // Создаем массив B из n+1 строки
    for(size_t i=0;i<=n;++i) // Заполняем i-ю строку массива
    {
       B[i][0]=1;         // На концах строки стоят единицы
       B[i][i]=1;
       for(size_t j=1;j<i;++j)
       {   // Заполняем оставшиеся элементы i-й строки
           B[i][j]=B[i-1][j-1]+B[i-1][j];
       }
    }
    return B;
}

std::vector<double> RandomVariable::get_left_borders()
{
    std::vector<double> left_borders;

    std::vector<std::vector<size_t>> C = get_pascal_triangle(n);
    double left_border = 0;
    for (size_t i = 0; i < n+1; ++i)
    {
        left_border += C[n][i]*std::pow(p, i)*std::pow(1-p,n-i);
        left_borders.push_back(left_border);
    }
    return left_borders;
}

std::vector<double> RandomVariable::get_probabilities()
{
    std::vector<double> probabilities;
    std::vector<std::vector<size_t>> C = get_pascal_triangle(n);

    for (size_t i = 0; i < n+1; ++i)
        probabilities.push_back(C[n][i]*std::pow(p, i)*std::pow(1-p,n-i));

    return probabilities;
}

double RandomVariable::get_expected_value()
{
    double expected_value = 0;
    for (size_t i = 0; i < probabilities.size(); ++i)
        expected_value += i*probabilities[i];
    return expected_value;
}

double RandomVariable::get_sample_expected_value()
{
    double sample_expected_value = 0;
    for (size_t i = 0; i < random_variables.size(); ++i)
        sample_expected_value += random_variables[i];
    return sample_expected_value / random_variables.size();
}

double RandomVariable::get_variance()
{
    double expected_value = get_expected_value();
    double variance = 0;
    for (size_t i = 0; i < probabilities.size(); ++i)
        variance += std::pow(i - expected_value, 2) * probabilities[i];
    return variance;
}

double RandomVariable::get_sample_variance()
{
    double sample_expected_value = get_sample_expected_value();
    double sample_variance = 0;
    for (size_t i = 0; i < random_variables.size(); ++i)
        sample_variance += std::pow(random_variables[i] - sample_expected_value, 2);
    return sample_variance / random_variables.size();
}

double RandomVariable::get_sample_median()
{
    size_t sz = random_variables.size();
    double res = 0;
    if (sz % 2 == 1)
        res = random_variables[sz/2];
    else
        res = (random_variables[sz/2-1] + random_variables[sz/2]) / 2.0;

    return res;
}

double RandomVariable::get_range()
{
    double res = random_variables.back() - random_variables.front();
    return res;
}

std::vector<size_t> RandomVariable::count_random_variables()
{
    std::vector<size_t> random_variables;
    std::srand(time(0));

    std::vector<double> left_borders = get_left_borders();
    left_borders.pop_back();
    std::cout << "left borders\n";
    for (const auto& p : left_borders)
        std::cout << p << ' ';
    std::cout << '\n';

    for (size_t i = 0; i < experiments; ++i)
    {
        const double u = ((double) std::rand()) / RAND_MAX;
        auto it = std::lower_bound(left_borders.begin(), left_borders.end(), u);
        const size_t random_value = it - left_borders.begin();
        random_variables.push_back(random_value);
    }

    return random_variables;
}

void RandomVariable::simulate_experiment()
{
    std::vector<size_t> rv;
    for (size_t i = 0; i < experiments; ++i)
    {
        size_t count = 0;
        for (size_t j = 0; j < n; ++j)
        {
            const double u = ((double) std::rand()) / RAND_MAX;
            count += u <= p ? 1 : 0;
        }
        rv.push_back(count);
    }
    random_variables = rv;
}

const std::vector<size_t>& RandomVariable::get_random_variables()
{
    return random_variables;
}

std::map<size_t, size_t> RandomVariable::get_frequencies(const std::vector<size_t>& rvs)
{
    std::map<size_t, size_t> res;
    for (const auto rv : rvs)
        res[rv]++;

    return res;
}
