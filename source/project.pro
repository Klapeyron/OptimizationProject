#-------------------------------------------------
#
# Project created by QtCreator 2017-03-31T21:03:53
#
#-------------------------------------------------

QT       += core gui \
         printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = project
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += main.cpp\
        mainwindow.cpp \
        qcustomplot.cpp \
    algorithm.cpp \
    random.cpp \
    muParser.cpp \
    muParserBase.cpp \
    muParserBytecode.cpp \
    muParserCallback.cpp \
    muParserDLL.cpp \
    muParserError.cpp \
    muParserInt.cpp \
    muParserTest.cpp \
    muParserTokenReader.cpp

HEADERS  += \
        qcustomplot.h \
    algorithm.h \
    random.hpp \
    muParser.h \
    muParserBase.h \
    muParserBytecode.h \
    muParserCallback.h \
    muParserDef.h \
    muParserDLL.h \
    muParserError.h \
    muParserFixes.h \
    muParserInt.h \
    muParserStack.h \
    muParserTemplateMagic.h \
    muParserTest.h \
    muParserToken.h \
    muParserTokenReader.h

FORMS    += mainwindow.ui

CONFIG += c++14
