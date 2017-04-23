#include "exprtk.hpp"
#include "algorithm.h"

#include <regex>
#include <map>
#include <set>

Point::Point(Function& firstFunction, Function& secondFunction, std::map<std::string, double> const& initialSymbols)
    :firstFunction(firstFunction),
     secondFunction(secondFunction)
{
    auto applySymbols = [&](Function& function)
    {
        auto functionSymbols = function.getSymbols();

        for(auto const& symbol : functionSymbols)
        {
            try
            {
                auto initialValue = initialSymbols.at(symbol);
                function.setValue(symbol, initialValue);
                symbols[symbol] = initialValue;
            }
            catch(std::out_of_range & exception)
            {
                symbols[symbol] = 0.0;
            }
        }
    };

    applySymbols(firstFunction);
    applySymbols(secondFunction);

    xValue = firstFunction.calculateExpression();
    yValue = secondFunction.calculateExpression();
}

double Point::getXValue() const
{
    return xValue;
}

double Point::getYValue() const
{
    return yValue;
}

const std::map<std::string, double>& Point::getSymbols() const
{
    return symbols;
}

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

std::set<std::string> Function::getSymbols()
{
    std::set<std::string> symbols;

    for(auto const& symbol : uniqueSymbols_)
    {
        symbols.insert(symbol.first);
    }
    return symbols;
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
    return expression_.value();
}

void Algorithm::clearCell(int row, int column)
{
    // ui->tableWidget->item(row, column)->setText(QString(""));
}

void Algorithm::updateConstraints(int row, int column)
{
    auto symbol = std::string(ui->tableWidget->item(row, 0)->text().toUtf8().constData());
    if(column == 1) // min value changed
    {
        constraints[symbol].min = std::stod(std::string(ui->tableWidget->item(row, 1)->text().toUtf8().constData()));
    }
    else if(column == 2) // max value changed
    {
        constraints[symbol].max = std::stod(std::string(ui->tableWidget->item(row, 2)->text().toUtf8().constData()));
    }
}

void Algorithm::putSymbolsToTable()
{
    auto firstFunctionSymbols = Function(std::string(ui->function1->text().toUtf8().constData())).getSymbols();
    auto secondFunctionSymbols = Function(std::string(ui->function2->text().toUtf8().constData())).getSymbols();

    std::set<std::string> symbols;

    std::merge(firstFunctionSymbols.begin(), firstFunctionSymbols.end(),
               secondFunctionSymbols.begin(), secondFunctionSymbols.end(),
               std::inserter(symbols, symbols.end()));

    ui->tableWidget->clear();

    ui->tableWidget->setRowCount(symbols.size());
    ui->tableWidget->setColumnCount(4);

    ui->tableWidget->setHorizontalHeaderLabels(QStringList{"Symbol", "Min value", "Max value", "Result"});
    // ui->tableWidget->setVerticalHeaderLabels(QStringList{"x1", "x2", "x3", "x4"});

    auto i = 0u;
    for(auto const& symbol : symbols)
    {
        auto* item = new QTableWidgetItem(QString::fromStdString(symbol));
        ui->tableWidget->setItem(i, 0, item);
        item->setFlags(item->flags() ^ Qt::ItemIsEditable);

        try
        {
            auto* minItem = new QTableWidgetItem(QString::fromStdString(std::to_string(constraints.at(symbol).min)));
            ui->tableWidget->setItem(i, 1, minItem);
        } catch(std::out_of_range&) {}

        try
        {
            auto* maxItem = new QTableWidgetItem(QString::fromStdString(std::to_string(constraints.at(symbol).max)));
            ui->tableWidget->setItem(i, 2, maxItem);
        } catch(std::out_of_range&) {}

        i++;
    }
}

bool isDominated(std::vector<Point> setOfPoints, Point comparePoint)
{
    return std::find_if(setOfPoints.begin(), setOfPoints.end(), [&](Point const& pointInCollection)
    {
        return comparePoint.getXValue() > pointInCollection.getXValue() and comparePoint.getYValue() > pointInCollection.getYValue();
    }) != setOfPoints.end();
}

double norm(Point firstPoint, Point secondPoint)
{
    auto xLength = firstPoint.getXValue() - secondPoint.getXValue();
    auto yLength = firstPoint.getYValue() - secondPoint.getYValue();

    return std::sqrt(std::pow(xLength, 2) + std::pow(yLength, 2));
}

std::size_t numberOfPointsInCircleWithNicheRadius(std::vector<Point> compareSet, Point point, double nicheRadius)
{
    std::size_t numberOfPoints = 0;
    for(auto const& comparePoint : compareSet)
    {
        if (norm(comparePoint, point) < nicheRadius)
        {
            numberOfPoints++;
        }
    }
    return numberOfPoints;
}

std::vector<Point> Algorithm::generatePoints(std::size_t numberOfPoints)
{
    std::vector<Point> points;
    for(std::size_t i = 0; i < numberOfPoints; ++i)
    {
        auto firstFunctionSymbols = firstFunction->getSymbols();
        auto secondFunctionSymbols = secondFunction->getSymbols();

        std::set<std::string> functionSymbols;
        std::merge(firstFunctionSymbols.begin(), firstFunctionSymbols.end(),
                   secondFunctionSymbols.begin(), secondFunctionSymbols.end(),
                   std::inserter(functionSymbols, functionSymbols.end()));

        std::map<std::string, double> generatedSymbols;

        for(auto const& symbol : functionSymbols)
        {
            auto min = constraints[symbol].min;
            auto max = constraints[symbol].max;

            auto generatedSymbolValue = generator.generateDouble(min, max);

            generatedSymbols[symbol] = generatedSymbolValue;
        }
        auto generatedPoint = Point(*firstFunction, *secondFunction, std::move(generatedSymbols));
        points.push_back(std::move(generatedPoint));
    }

    return points;
}

Point Algorithm::crossover(Point const& firstPoint, Point const& secondPoint)
{
    std::map<std::string, double> newSymbols;

    auto pointSymbolsMap = firstPoint.getSymbols();

    auto symbolIt = pointSymbolsMap.begin();
    for(unsigned i = 0; i != pointSymbolsMap.size()/2; i++)
    {
        symbolIt++;
        const auto  symbolName = (*symbolIt).first;

        auto newValue = secondPoint.getSymbols().at(symbolName);

        newSymbols[symbolName] = newValue;
    }

    return Point(*firstFunction, *secondFunction, newSymbols);
}

Point Algorithm::mutate(Point const& point)
{
    auto symbolsOfBasedPoint = point.getSymbols();

    auto symbolIndex = generator.generateInt(0, symbolsOfBasedPoint.size() - 1);

    auto symbolIt = symbolsOfBasedPoint.begin();
    std::advance(symbolIt, symbolIndex);

    auto constraint = constraints[(*symbolIt).first];
    auto multiplicand = generator.generateDouble(-1, 1);
    auto newValue = (*symbolIt).second + (*symbolIt).second * multiplicand;

    bool underConstaints = newValue >= constraint.min && newValue <= constraint.max;
    if(underConstaints)
    {
        (*symbolIt).second = newValue;
    }
    else
    {
        (*symbolIt).second = generator.generateDouble(constraint.min, constraint.max);
    }

    return Point(*firstFunction, *secondFunction, std::move(symbolsOfBasedPoint));
}

void Algorithm::startCalculations()
{
    firstFunction = std::make_unique<Function>(std::string(ui->function1->text().toUtf8().constData()));
    secondFunction = std::make_unique<Function>(std::string(ui->function2->text().toUtf8().constData()));

    // ui->function2->setText(QString::fromStdString(std::to_string(firstFunction.calculateExpression())) +
    //                       QString::fromStdString(std::to_string(secondFunction.calculateExpression())));

    auto pm = ui->mutationProbability->text().toDouble();
    auto pc = ui->crossingProbability->text().toDouble();
    auto T = ui->maxGenertation->text().toUInt();
    auto N = ui->populationSize->text().toUInt();
    auto sigma = ui->sigma->text().toDouble();

    auto tdom = 7u;

    std::vector<Point> p0 = generatePoints(N);

    for (unsigned t = 0; t < T; ++t)
    {
        std::vector<Point> temporarySet;

        for(std::size_t i = 0; i < N; i++)
        {
            auto firstPoint = p0[generator.generateInt(0, p0.size() - 1)];
            auto secondPoint = p0[generator.generateInt(0, p0.size() - 1)];

            std::vector<Point> compareSet = generatePoints(tdom);

            if(not isDominated(compareSet, firstPoint) and isDominated(compareSet, secondPoint))
            {
                temporarySet.push_back(firstPoint);
            }
            else if(isDominated(compareSet, firstPoint) and not isDominated(compareSet, secondPoint))
            {
                temporarySet.push_back(secondPoint);
            }
            else
            {
                auto nFirstPoint = numberOfPointsInCircleWithNicheRadius(compareSet, firstPoint, sigma);
                auto nSecondPoint = numberOfPointsInCircleWithNicheRadius(compareSet, secondPoint, sigma);

                if (nFirstPoint < nSecondPoint)
                {
                    temporarySet.push_back(firstPoint);
                }
                else
                {
                    temporarySet.push_back(secondPoint);
                }
            }
        }

        // TODO: VEGA 3, 4, 5
        std::vector<Point> crossoverSet;

        for(std::size_t i = 0; i < std::floor(N/2); i++)
        {
            auto firstPos = generator.generateInt(0, temporarySet.size() - 1);
            auto firstPoint = temporarySet[firstPos];
            // temporarySet.erase(temporarySet.begin() + firstPos);

            auto secondPos = generator.generateInt(0, temporarySet.size() - 1);
            auto secondPoint = temporarySet[secondPos];
            // temporarySet.erase(temporarySet.begin() + secondPos);

            auto decisionVariable = generator.generateDouble(0,100);
            bool shouldCross = decisionVariable < pc;

            if(shouldCross)
            {
                auto newPoint1 = crossover(firstPoint, secondPoint);
                auto newPoint2 = crossover(secondPoint, firstPoint);
                crossoverSet.push_back(std::move(newPoint1));
                crossoverSet.push_back(std::move(newPoint2));
            }
            else
            {
                crossoverSet.push_back(std::move(firstPoint));
                crossoverSet.push_back(std::move(secondPoint));
            }
        }

        std::vector<Point> mutationSet;

        for (auto const& crossoverPoint : crossoverSet)
        {
            auto decisionVariable = generator.generateDouble(0,100);
            bool shouldMutate = decisionVariable < pm;

            if(shouldMutate)
            {
                auto newPoint = mutate(crossoverPoint);
                mutationSet.push_back(std::move(newPoint));
            }
            else
            {
                mutationSet.push_back(crossoverPoint);
            }
        }

        p0.clear();
        p0.swap(mutationSet);
    }
}
