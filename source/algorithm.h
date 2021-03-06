#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <string>
#include <QMainWindow>
#include <QObject>
#include <ui_mainwindow.h>
#include <memory>
#include <set>

#include "random.hpp"
#include "muParser.h"

namespace Ui {
  class MainWindow;
}

class Function
{
    mu::Parser parser_;
    std::map<std::string, double> uniqueSymbols_;
public:
    Function(std::string const& expression);
    bool setValue(std::string const& symbol, double const& newValue);
    std::set<std::string> getSymbols();
    double calculateExpression();
};

class Point
{
public:
    Point(std::shared_ptr<Function> firstFunction,
          std::shared_ptr<Function> secondFunction,
          std::map<std::string, double> const& symbols);
    const std::map<std::string, double>& getSymbols() const;
    double getXValue() const;
    double getYValue() const;
    std::shared_ptr<Function> getFirstFunction() const;
    std::shared_ptr<Function> getSecondFunction() const;
private:
    std::shared_ptr<Function> firstFunction;
    std::shared_ptr<Function> secondFunction;

    double xValue, yValue;
    std::map<std::string, double> symbols;
};

struct Constraint
{
    Constraint(std::string min = "0.0", std::string max = "0.0") :minString(min), maxString(max) {}
    double min = 0.0, max = 0.0;
    std::string minString, maxString;

    bool calculateConstraints();
};

class Algorithm :public QObject
{
    Q_OBJECT

    Ui::MainWindow* ui;

    std::map<std::string, Constraint> constraints;
    RandomGenerator generator;
    std::shared_ptr<Function> firstFunction, secondFunction;

    void printPoints(const std::vector<Point>& points);
    void tabularizePoints(std::vector<Point>& points);
    std::vector<Point> chooseRandomPointsFromSet(const std::vector<Point>& points,
                                                 const unsigned& amountOfPoints);
public:
    Algorithm(Ui::MainWindow* ui) :ui(ui) {}
    virtual ~Algorithm() {}

    std::vector<Point> generatePoints(std::size_t numberOfPoints);

    std::pair<Point, Point> crossover(Point const& firstPoint, Point const& secondPoint);
    Point mutate(Point const& firstPoint);
public slots:
    void startCalculations();
    void putSymbolsToTable();
    void updateConstraints(int row, int column);
    void clearCell(int row, int column);
};

#endif // ALGORITHM_H
