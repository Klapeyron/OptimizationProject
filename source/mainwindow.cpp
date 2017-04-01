#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    algorithm(ui)
{
    ui->setupUi(this);
    connect(ui->startButton, SIGNAL(clicked(bool)), &algorithm, SLOT(startCalculations()));
    connect(ui->function1, SIGNAL(editingFinished()), &algorithm, SLOT(putSymbolsToTable()));
    connect(ui->function2, SIGNAL(editingFinished()), &algorithm, SLOT(putSymbolsToTable()));
}

MainWindow::~MainWindow()
{
    delete ui;
}
