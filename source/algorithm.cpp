#include "exprtk.hpp"
#include "algorithm.h"

#include <regex>
#include <map>
#include <set>
#include <thread>
#include <chrono>

Point::Point(std::shared_ptr<Function> firstFunction, std::shared_ptr<Function> secondFunction, std::map<std::string, double> const& initialSymbols)
    :firstFunction(firstFunction),
     secondFunction(secondFunction)
{
    auto applySymbols = [&](std::shared_ptr<Function> function)
    {
        auto functionSymbols = function->getSymbols();

        for(auto const& symbol : functionSymbols)
        {
            try
            {
                auto initialValue = initialSymbols.at(symbol);
                function->setValue(symbol, initialValue);
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

    xValue = firstFunction->calculateExpression();
    yValue = secondFunction->calculateExpression();
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

bool Constraint::calculateConstraints()
{
    try
    {
        min = std::stod(minString);
        max = std::stod(maxString);
        if(max < min)
            return false;
        return true;
    } catch(std::invalid_argument&)
    {
        return false;
    }
}

void Algorithm::clearCell(int row, int column)
{
    auto* newClearItem = new QTableWidgetItem(QString(""));
    ui->tableWidget->setItem(row, column, newClearItem);
}

void Algorithm::updateConstraints(int row, int column)
{
    auto symbol = std::string(ui->tableWidget->item(row, 0)->text().toUtf8().constData());
    if(column == 1) // min value changed
    {
        constraints[symbol].minString = std::string(ui->tableWidget->item(row, 1)->text().toUtf8().constData());
    }
    else if(column == 2) // max value changed
    {
        constraints[symbol].maxString = std::string(ui->tableWidget->item(row, 2)->text().toUtf8().constData());
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
    ui->tableWidget->setColumnCount(3);

    ui->tableWidget->setHorizontalHeaderLabels(QStringList{"Symbol", "Min value", "Max value"});

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
        auto generatedPoint = Point(firstFunction, secondFunction, std::move(generatedSymbols));
        points.push_back(std::move(generatedPoint));
    }

    return points;
}

std::pair<Point, Point> Algorithm::crossover(Point const& firstPoint, Point const& secondPoint)
{   
    std::map<std::string, double> newFirstSymbols, newSecondSymbols;

    auto allSymbols = firstPoint.getSymbols();

    for (auto const& symbolPair : allSymbols)
    {
        auto & symbol = symbolPair.first;
        auto decisionVariable = generator.generateInt(0, 1);
        bool selectFirst = decisionVariable == 0;
        if(selectFirst)
        {
            newFirstSymbols[symbol] = firstPoint.getSymbols().at(symbol);
            newSecondSymbols[symbol] = secondPoint.getSymbols().at(symbol);
        }
        else
        {
            newFirstSymbols[symbol] = secondPoint.getSymbols().at(symbol);
            newSecondSymbols[symbol] = firstPoint.getSymbols().at(symbol);
        }
    }

    return std::make_pair<Point, Point>(Point(firstFunction, secondFunction, newFirstSymbols),
                                        Point(firstFunction, secondFunction, newSecondSymbols));
}

std::vector<Point> Algorithm::chooseRandomPointsFromSet(const std::vector<Point>& points,
                                                        const unsigned& amountOfPoints)
{
    std::vector<unsigned> indexVector;
    unsigned tempIndex;

    for(auto i = 0u; i != amountOfPoints; i++)
    {
        do
        {
            tempIndex = generator.generateInt(0, points.size() - 1);
        } while(std::find(indexVector.begin(), indexVector.end(), tempIndex) != indexVector.end());
        indexVector.push_back(tempIndex);
    }

    std::vector<Point> result;
    for (auto const& index : indexVector)
    {
        result.push_back(points[index]);
    }

    return result;
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

    auto underConstaints = [&](auto newValue)
    {
        return newValue >= constraint.min and newValue <= constraint.max;
    };

    while(not underConstaints(newValue))
    {
        multiplicand = generator.generateDouble(-1, 1);
        newValue = (*symbolIt).second + (*symbolIt).second * multiplicand;
    }

    (*symbolIt).second = newValue;

    auto newPoint = Point(firstFunction, secondFunction, std::move(symbolsOfBasedPoint));

    if(not isDominated({point}, newPoint))
    {
        return newPoint;
    }
    else
    {
        return point;
    }
}

void Algorithm::printPoints(const std::vector<Point>& points)
{
    ui->figurePlot->addGraph();

    QVector<double> xValue(points.size()), yValue(points.size());

    for (auto const& point : points)
    {
        xValue.push_back(point.getXValue());
        yValue.push_back(point.getYValue());
    }

    auto minmaxf1 = std::minmax_element(points.begin(), points.end(),
                                  [](const auto& point1, const auto& point2)
    {
        return point1.getXValue() < point2.getXValue();
    });
    auto minmaxf2 = std::minmax_element(points.begin(), points.end(),
                                  [](const auto& point1, const auto& point2)
    {
        return point1.getYValue() < point2.getYValue();
    });

    ui->figurePlot->graph(0)->setData(xValue, yValue);
    ui->figurePlot->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->figurePlot->graph(0)->setScatterStyle(QCPScatterStyle::ssDot);
    ui->figurePlot->graph(0)->setPen(QPen(QBrush(Qt::blue), 4));

    ui->figurePlot->xAxis->setLabel("f1");
    ui->figurePlot->yAxis->setLabel("f2");

    ui->figurePlot->xAxis->setRange(minmaxf1.first->getXValue(), minmaxf1.second->getXValue());
    ui->figurePlot->yAxis->setRange(minmaxf2.first->getYValue(), minmaxf2.second->getYValue());

    ui->figurePlot->replot();
}

void Algorithm::tabularizePoints(std::vector<Point>& points)
{
    std::sort(points.begin(), points.end(),
              [](Point const& left, Point const& right)
    {
        return left.getXValue() < right.getXValue();
    });

    auto firstFunctionSymbols = firstFunction->getSymbols();
    auto secondFunctionSymbols = secondFunction->getSymbols();

    std::set<std::string> symbols;

    std::merge(firstFunctionSymbols.begin(), firstFunctionSymbols.end(),
               secondFunctionSymbols.begin(), secondFunctionSymbols.end(),
               std::inserter(symbols, symbols.end()));

    ui->resultsTable->clear();

    ui->resultsTable->setRowCount(points.size());
    ui->resultsTable->setColumnCount(symbols.size()+2);

    auto tableColumnsLabels = QStringList{"f1", "f2"};
    for (auto const& symbol : symbols)
    {
        tableColumnsLabels.push_back(QString::fromStdString(symbol));
    }
    ui->resultsTable->setHorizontalHeaderLabels(tableColumnsLabels);

    auto i = 0u;
    for(auto const& point : points)
    {
        try
        {
            auto* f1Item = new QTableWidgetItem(QString::fromStdString(std::to_string(point.getXValue())));
            ui->resultsTable->setItem(i, 0, f1Item);
            f1Item->setFlags(f1Item->flags() ^ Qt::ItemIsEditable);

            auto* f2Item = new QTableWidgetItem(QString::fromStdString(std::to_string(point.getYValue())));
            ui->resultsTable->setItem(i, 1, f2Item);
            f2Item->setFlags(f2Item->flags() ^ Qt::ItemIsEditable);
        } catch(std::out_of_range&) {}

        auto j = 2u;
        for (auto const& symbol : symbols)
        {
            try
            {
                auto* symbolItem = new QTableWidgetItem(QString::fromStdString(std::to_string(point.getSymbols().at(symbol))));
                ui->resultsTable->setItem(i, j, symbolItem);
                symbolItem->setFlags(symbolItem->flags() ^ Qt::ItemIsEditable);
            }
            catch (std::out_of_range&) {}
            j++;
        }
        i++;
    }
}

void Algorithm::startCalculations()
{
    for (auto& constraint : constraints)
    {
        if(not constraint.second.calculateConstraints())
        {
            // TODO: do something if invalid input in constraints
            QMessageBox msgBox;
            msgBox.setWindowTitle("Invalid constraints");
            msgBox.setText("Invalid constraints for " + QString::fromStdString(constraint.first));
            msgBox.exec();
            return;
        }
    }

    firstFunction = std::make_shared<Function>(std::string(ui->function1->text().toUtf8().constData()));
    secondFunction = std::make_shared<Function>(std::string(ui->function2->text().toUtf8().constData()));

    auto pm = ui->mutationProbability->text().toDouble();
    auto pc = ui->crossingProbability->text().toDouble();
    auto T = ui->maxGenertation->text().toUInt();
    auto N = ui->populationSize->text().toUInt();
    auto sigma = ui->sigma->text().toDouble();
    auto tdom = ui->tdom->text().toUInt();

    std::vector<Point> p0 = generatePoints(N);

    for (unsigned t = 0; t < T; ++t)
    {
        printPoints(p0);
        std::vector<Point> temporarySet;

        for(std::size_t i = 0; i < N; i++)
        {
            auto firstPoint = p0[generator.generateInt(0, p0.size() - 1)];
            auto secondPoint = p0[generator.generateInt(0, p0.size() - 1)];

            std::vector<Point> compareSet = chooseRandomPointsFromSet(p0, tdom);

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

//        std::sort(temporarySet.begin(), temporarySet.end(),
//                  [](Point const& left, Point const& right)
//        {
//            return left.getXValue() < right.getXValue();
//        });
        std::vector<Point> crossoverSet;

        for(std::size_t i = 0; i < std::floor(N/2); i++)
        {
            auto firstPos = generator.generateInt(0, temporarySet.size() - 1);
            auto firstPoint = temporarySet[firstPos];

            unsigned secondPos;
            do
            {
                secondPos = generator.generateInt(0, temporarySet.size() - 1);
            } while(firstPos == secondPos);
            auto secondPoint = temporarySet[secondPos];

            temporarySet.erase(temporarySet.begin() + secondPos);
            temporarySet.erase(temporarySet.begin() + firstPos);

            auto decisionVariable = generator.generateDouble(0, 100);
            bool shouldCross = decisionVariable < pc;

            if(shouldCross)
            {
                auto pointsPair = crossover(firstPoint, secondPoint);

                crossoverSet.push_back(std::move(pointsPair.first));
                crossoverSet.push_back(std::move(pointsPair.second));
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
            auto decisionVariable = generator.generateDouble(0, 100);
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
//        std::this_thread::sleep_for(std::chrono::milliseconds(300));
        ui->progressBar->setValue((t * 100/T) + 1);
    }
    auto pointsToPrint = std::vector<Point> {};
    for (auto const& point : p0)
    {
        bool dominated = isDominated(p0, point);
        if (not dominated)
        {
            pointsToPrint.push_back(point);
        }
    }
    printPoints(pointsToPrint);
    ui->progressBar->setValue(100);
    tabularizePoints(pointsToPrint);
}
