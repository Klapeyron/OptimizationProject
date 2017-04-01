#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <string>
#include "exprtk.hpp"
#include <QMainWindow>
#include <QObject>
#include <ui_mainwindow.h>
#include "random.hpp"

using symbol_table_t = exprtk::symbol_table<double>;
using expression_t   = exprtk::expression<double>;
using parser_t       = exprtk::parser<double>;

namespace Ui {
  class MainWindow;
}

class Function
{
    symbol_table_t symbol_table_;
    parser_t parser_;
    expression_t expression_;
    std::map<std::string, double> uniqueSymbols_;
    RandomGenerator generator;
public:
    Function(std::string const& expression);

    bool setValue(std::string const& symbol, double const& newValue);
    double calculateExpression();
};

class Algorithm :public QObject
{
    Q_OBJECT

    Ui::MainWindow* ui;
public:
    Algorithm(Ui::MainWindow* ui) :ui(ui) {}
    virtual ~Algorithm() {};
public slots:
    void startCalculations();
    void putSymbolsToTable();
};

#endif // ALGORITHM_H
