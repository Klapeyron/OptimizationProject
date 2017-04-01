#include "exprtk.hpp"
#include "algorithm.h"

#include <regex>
#include <map>

Function::Function(std::string const& expression)
{
    std::regex symbols("x[0-9]+");

    auto begin = std::sregex_iterator(expression.begin(), expression.end(), symbols);
    auto end = std::sregex_iterator();

    for (std::sregex_iterator i = begin; i != end; ++i)
    {
        uniqueSymbols_[(*i).str()] = 0.0;
    }

    for(auto & symbolPair : uniqueSymbols_)
    {
        symbol_table_.add_variable(symbolPair.first, uniqueSymbols_[symbolPair.first]);
    }

    expression_.register_symbol_table(symbol_table_);
    parser_.compile(expression, expression_);
}

bool Function::setValue(std::string const& symbol, double const& newValue)
{
    auto it = uniqueSymbols_.find(symbol);
    if(it == uniqueSymbols_.end())
        return false;
    it->second = newValue;
    return true;
}

double Function::calculateExpression()
{
    expression_.value();
}

void Algorithm::putSymbolsToTable()
{
    ui->tableWidget->setRowCount(6);
    ui->tableWidget->setColumnCount(4);
}

void Algorithm::startCalculations()
{
    Function firstFunction(std::string(ui->function1->text().toUtf8().constData()));
    Function secondFunction(std::string(ui->function2->text().toUtf8().constData()));

    firstFunction.setValue("x2", 14);

    ui->function2->setText(QString::fromStdString(std::to_string(firstFunction.calculateExpression())) +
                           QString::fromStdString(std::to_string(secondFunction.calculateExpression())));
}
