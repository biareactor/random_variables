#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include <iostream>
#include "random_variable.h"
#include <QChart>
#include <QLineSeries>
#include <QString>
#include <QtMath>
#include <numeric>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

static void set_axes(QAbstractAxis* axisX, QAbstractAxis* axisY, QString label1, QString label2, size_t fontsize)
{
    auto font = QFont();
    font.setPixelSize(24);

    axisX->setTitleText(label1);
    axisX->setTitleFont(font);

    axisY->setTitleText(label2);
    axisY->setTitleFont(font);
}

double Gamma(double a)
{
    if (std::abs(a - 1) < 1e-10)
        return 1;
    if (std::abs(a - 0.5) < 1e-10)
        return qSqrt(M_PI);
    return ((a - 1) * Gamma(a - 1));
}

double fx(double x, double r)
{
    if (x > 0)
        return qPow(2, -r/2) * qPow(Gamma(r/2), -1) * qPow(x, r/2 - 1) * qExp(-x/2);
    return 0;
}

double F_R0(double R0, double k, double n, double N)
{
    double SUM = 0;
    for (size_t i = 1; i < n+1; i++)
        SUM += R0 * (fx(R0 * ((i-1)/n), k - 1) + fx(R0 * (i/n),k - 1))/(2*N);
    return 1 - SUM;
}

void MainWindow::on_calculate_clicked()
{
    while (ui->table->columnCount() > 0)
    {
        ui->table->removeColumn(0);
    }

    while (ui->table2->rowCount() > 0)
    {
        ui->table2->removeRow(0);
    }

    while (ui->table3->columnCount() > 0)
    {
        ui->table3->removeColumn(0);
    }

    const size_t experiments = ui->experiments->text().toUInt();
    const size_t n = ui->n->text().toUInt();
    const double p = ui->p->text().toDouble();

    RandomVariable rv(experiments, n, p);

    auto random_variables = rv.get_random_variables();
//    auto random_variables = rv.simulate_experiment();
    auto freqs = rv.get_frequencies(random_variables);

    for (size_t i = 0; i < n+1; i++)
    {
        ui->table->insertColumn(i);

        ui->table->setItem(0, i, new QTableWidgetItem(QString::number(i)));
        ui->table->setItem(1, i, new QTableWidgetItem(QString::number(freqs[i])));
        ui->table->setItem(2, i, new QTableWidgetItem(QString::number(1.0*freqs[i]/experiments)));
    }

    ui->table2->insertRow(0);
    ui->table2->setItem(0, 0, new QTableWidgetItem(QString::number(rv.get_expected_value())));
    ui->table2->setItem(0, 1, new QTableWidgetItem(QString::number(rv.get_sample_expected_value())));
    ui->table2->setItem(0, 2, new QTableWidgetItem(QString::number(std::fabs(rv.get_expected_value() - rv.get_sample_expected_value()))));
    ui->table2->setItem(0, 3, new QTableWidgetItem(QString::number(rv.get_variance())));
    ui->table2->setItem(0, 4, new QTableWidgetItem(QString::number(rv.get_sample_variance())));
    ui->table2->setItem(0, 5, new QTableWidgetItem(QString::number(std::fabs(rv.get_variance() - rv.get_sample_variance()))));
    ui->table2->setItem(0, 6, new QTableWidgetItem(QString::number(rv.get_sample_median())));
    ui->table2->setItem(0, 7, new QTableWidgetItem(QString::number(rv.get_range())));

    double D = -1;

    auto* chart = new QChart();
    auto* series = new QLineSeries();
    series->setColor(QColorConstants::Blue);

    auto lb = rv.get_left_borders();

    for (int i = 0; i < n+1; i++)
    {
        series->append(i, lb[i]);
        series->append(i+0.999, lb[i]);

        chart->addSeries(series);
        series = new QLineSeries();
        series->setColor(QColorConstants::Green);

        size_t c = 0;
        size_t c_prev = 0;
        for (const auto rv : random_variables)
        {
            if (rv < i) c++;
            if (rv < i-1) c_prev++;
        }

        series->append(i-0.999, 1.0*c/experiments);
        series->append(i, 1.0*c/experiments);

        chart->addSeries(series);
        series = new QLineSeries();
        series->setColor(QColorConstants::Blue);

        double x = std::max(1.0*c/experiments - lb[i+1], lb[i] - 1.0*c_prev/experiments);//std::fabs(lb[i] - 1.0*c/experiments);
        if (x > D) D = x;
    }

    ui->D->setText(QString::number(D));
    chart->createDefaultAxes();

    auto axisX = chart->axes(Qt::Horizontal).back();
    auto axisY = chart->axes(Qt::Vertical).back();
    set_axes(axisX, axisY, "x", "Fη(x)", 24);

    chart->legend()->hide();
    ui->F->setChart(chart);

    /**********************************************/

    auto probs = rv.get_probabilities();
    double max_freq_prob = -1;

    std::vector<double> yi, pi, ni;

    for (size_t i = 0; i < n+1; i++)
    {
        ui->table3->insertColumn(i);

        yi.push_back(i);

        ui->table3->setItem(0, i, new QTableWidgetItem(QString::number(i)));
        ui->table3->setItem(1, i, new QTableWidgetItem(QString::number(probs[i])));
        ui->table3->setItem(2, i, new QTableWidgetItem(QString::number(1.0*freqs[i]/experiments)));

        if (std::fabs(1.0*freqs[i]/experiments - probs[i]) > max_freq_prob)
            max_freq_prob = std::fabs(1.0*freqs[i]/experiments - probs[i]);
    }

    ui->max_freq_prob->setText(QString::number(max_freq_prob));

    /************************************************/

    auto intervals_txt = ui->intervals->text().split(" ");
    std::vector<double> intervals;
    for (const auto& i : intervals_txt)
        intervals.push_back(i.toDouble());

    const size_t k = intervals.size();
    /*
    q = [0] * k
    ni = [0] * k

    it = 0
    for i in range(1, k + 2):
        while it < len(df) and df['yi'][it] < vi[i]:
            q[i-1] += df['Pi'][it]
            ni[i-1] += df['ni'][it]
            it += 1
    S = 0
    N = 0
    for i in range(k):
        S += q[i]
        N += ni[i]
        print(q[i], end = " ")
    print()
    print(S, N)*/

    std::vector<double> qi;
    std::vector<double> ni_;

    size_t j = 0;
    double q = 0;
    double n_ = 0;

    for (size_t i = 0; i < yi.size(); )
    {
        if (j == intervals.size() || yi[i] < intervals[j])
        {
            q += probs[i];
            n_ += freqs[i];
            i++;
        }
        else
        {
            qi.push_back(q);
            ni_.push_back(n_);

            q = 0;
            n_ = 0;
            j++;
        }
    }

    qi.push_back(q);
    ni_.push_back(n_);

    double S = std::accumulate(qi.begin(), qi.end(), 0.0);
    double N = std::accumulate(ni_.begin(), ni_.end(), 0.0);

    ui->S->setText(QString::number(S));
    ui->N->setText(QString::number(N));

    const double alpha = ui->alpha->text().toDouble();

    double R0 = 0;
    for (size_t i = 0; i < qi.size(); i++)
        R0 += ((ni_[i]-experiments*qi[i])*(ni_[i]-experiments*qi[i]))/(experiments*qi[i]);

    double Fres = F_R0(R0, k, n, experiments);
    ui->Fres->setText(QString::number(Fres));

    if(Fres > alpha) ui->hyp->setText("Гипотеза принята");
    else ui->hyp->setText("Гипотеза отклонена");
}
