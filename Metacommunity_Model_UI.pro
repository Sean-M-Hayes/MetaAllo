#-------------------------------------------------
#
# Project created by QtCreator 2015-03-24T20:48:37
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Metacommunity_Model_UI
TEMPLATE = app

QMAKE_CXXFLAGS += -std=gnu++0x -fopenmp
QMAKE_LFLAGS += -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp
QMAKE_CFLAGS_RELEASE += -fopenmp

CONFIG += static

SOURCES += main.cpp\
        mainwindow.cpp \
    allometric_metacommunity_class.cpp \
    opkdmain-lsoda.f \
    opkda1-lsoda.f \
    opkda2-lsoda.f \
    table_validator_delegate.cpp \
    food_web_methods.cpp \
    parameter_generator_functions.cpp \
    measure_time_series.cpp \
    lsoda_link.cpp \
    input_reading_functions.cpp \
    xml_methods.cpp \
    spatial_structure_methods.cpp \
    write_simulation_data_functions.cpp \
    find_cycles.cpp \
    read_csv.cpp \
    find_extrema.cpp \
    find_patterns.cpp \
    find_sequences.cpp \
    k_means.cpp

HEADERS  += mainwindow.h \
    allometric_metacommunity.h \
    lsoda_link.h \
    table_validator_delegate.h \
    food_web_methods.h \
    parameter_generator_functions.h \
    measure_time_series.h \
    input_reading_functions.h \
    forward_declarations.h \
    xml_methods.h \
    spatial_structure_methods.h \
    write_simulation_data_functions.h \
    find_cycles.h \
    read_csv.h \
    find_extrema.h \
    find_patterns.h \
    find_sequences.h \
    k_means.h

FORMS    += mainwindow.ui

INCLUDEPATH += /home/syn_mal/include/

LIBS += /usr/lib/libboost_thread.a \
        /usr/lib/libboost_random.a \
        /usr/lib/libboost_system.a \
        -lgfortran \
        -lquadmath \
        -fopenmp
