#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include "AudioFile.h"
#include <iostream>
#include <math.h>
#include <string>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;
using namespace std;
using namespace matplot;
int const X = 1500;
int const Y = 1500;

vector<double> filtr(const vector<double>& sygnal, double czestotliwosc_graniczna, double probkowanie, int rzadFiltru) {
    vector<double> przefiltrowanySygnal(sygnal.size());
    vector<double> wspolczynnikiFiltru(rzadFiltru);

    // Generowanie wspó³czynników filtra (prostok¹tny FIR)
    for (int i = 0; i < rzadFiltru; ++i) {
        if (i == rzadFiltru / 2) {
            wspolczynnikiFiltru[i] = 2 * czestotliwosc_graniczna / probkowanie;
        }
        else {
            wspolczynnikiFiltru[i] = sin(2 * pi * czestotliwosc_graniczna / probkowanie * (i - rzadFiltru / 2)) / (pi * (i - rzadFiltru / 2));
        }
    }

    // Normalizacja wspó³czynników filtra
    double sumaWspolczynnikow = accumulate(wspolczynnikiFiltru.begin(), wspolczynnikiFiltru.end(), 0.0);
    for (auto& wspolczynnik : wspolczynnikiFiltru) {
        wspolczynnik /= sumaWspolczynnikow;
    }

    // Zastosowanie filtra do sygna³u
    for (size_t i = rzadFiltru / 2; i < sygnal.size() - rzadFiltru / 2; ++i) {
        double przefiltrowanaWartosc = 0.0;
        for (int j = 0; j < rzadFiltru; ++j) {
            przefiltrowanaWartosc += sygnal[i - rzadFiltru / 2 + j] * wspolczynnikiFiltru[j];
        }
        przefiltrowanySygnal[i] = przefiltrowanaWartosc;
    }

    return przefiltrowanySygnal;
}
double prostokat(double x) {
    return (sin(x) >= 0) ? 1 : -1;
}
double pila(double x) {
    return 2 * (x / (2 * pi) - floor(0.5 + x / (2 * pi)));
}
void rysowanie_wykresow_podstawowych(int czestotliwosc, string typ) {
    vector<double> x(X);

    if (typ == "sin") {
        vector<double> y_sin(Y);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = i * 2.0 * pi / (x.size() - 1); // generowanie wartoœci x od 0 do 2*pi
            y_sin[i] = sin(czestotliwosc * x[i]); // obliczanie wartoœci y jako sin(x)
        }
        figure();
        plot(x, y_sin, "-"); // rysowanie wykresu sinusa
        title("Sinus");
    }
    else if (typ == "cos") {
        vector<double> y_cos(Y);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = i * 2.0 * pi / (x.size() - 1);
            y_cos[i] = cos(czestotliwosc * x[i]);
        }
        figure();
        plot(x, y_cos, "-");
        title("Cosinus");
    }
    else if (typ == "piloksztaltny") {
        vector<double> y_piloksztaltny(Y);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = i * 2.0 * pi / (x.size() - 1);
            y_piloksztaltny[i] = pila(czestotliwosc * x[i]);
        }
        figure();
        plot(x, y_piloksztaltny, "-");
        title("Pi³okszta³tny");
    }
    else if (typ == "prostokatny") {
        vector<double> y_prostokatny(Y);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = i * 2.0 * pi / (x.size() - 1);
            y_prostokatny[i] = prostokat(czestotliwosc * x[i]);
        }
        figure();
        plot(x, y_prostokatny, "-");
        title("Prostokatny (oryginalny)");
    }
    show();
}
void rysowanie_wykresow_filtr(double czestotliwosc, double czestotliwosc_graniczna, double probkowanie, int rzadFiltru, string typ) {
    vector<double> x(X);

    if (typ == "sin") {
        vector<double> y_sin(Y);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = i * 2.0 * pi / (x.size() - 1);
            y_sin[i] = sin(czestotliwosc * x[i]);
        }
        vector<double> filtered_y_sin = filtr(y_sin, czestotliwosc_graniczna, probkowanie, rzadFiltru);
        figure();
        plot(x, filtered_y_sin, "-r");
        title("Sinus (przefiltrowany)");
    }
    else if (typ == "cos") {
        vector<double> y_cos(Y);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = i * 2.0 * pi / (x.size() - 1);
            y_cos[i] = cos(czestotliwosc * x[i]);
        }
        vector<double> filtered_y_cos = filtr(y_cos, czestotliwosc_graniczna, probkowanie, rzadFiltru);
        figure();
        plot(x, filtered_y_cos, "-r");
        title("Cosinus (przefiltrowany)");
    }
    else if (typ == "piloksztaltny") {
        vector<double> y_piloksztaltny(Y);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = i * 2.0 * pi / (x.size() - 1);
            y_piloksztaltny[i] = pila(czestotliwosc * x[i]);
        }
        vector<double> filtered_y_piloksztaltny = filtr(y_piloksztaltny, czestotliwosc_graniczna, probkowanie, rzadFiltru);
        figure();
        plot(x, filtered_y_piloksztaltny, "-r");
        title("Piloksztaltny (przefiltrowany)");
    }
    else if (typ == "prostokatny") {
        vector<double> y_prostokatny(Y);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = i * 2.0 * pi / (x.size() - 1);
            y_prostokatny[i] = prostokat(czestotliwosc * x[i]);
        }
        vector<double> filtered_y_prostokatny = filtr(y_prostokatny, czestotliwosc_graniczna, probkowanie, rzadFiltru);
        figure();
        plot(x, filtered_y_prostokatny, "-r");
        title("Prostokatny (przefiltrowany)");
    }
    show();
}
void genAudio(double czestotliwosc_graniczna, int rzadFiltru) {
    AudioFile<double> audioFile; // Utwórz obiekt AudioFile

    cout << "Proba zaladowania pliku dzwiekowego..." << endl;
    bool zaladowano = audioFile.load("C:/Users/oneep/Source/Repos/scikit_build_example/audio.wav"); // SprawdŸ, czy plik znajduje siê w tym samym katalogu co AudioFile.h

    if (!zaladowano)
    {
        cout << "Nie udalo siê zaladowac pliku dzwiekowego." << endl;

    }
    else {
        cout << "Plik dzwiêkowy zaladowany pomyslnie." << endl;
        int probkowanie = audioFile.getSampleRate();
        int ile_probek = audioFile.getNumSamplesPerChannel();

        vector<double> t(ile_probek); // wektor na osi czasu
        vector<double> y(ile_probek); // wektor na osi y dla próbek dŸwiêkowych

        for (size_t i = 0; i < t.size(); ++i) {
            t[i] = (double)i / probkowanie; // generowanie wartoœci czasu
            y[i] = audioFile.samples[0][i]; // pobieranie próbek dŸwiêkowych
        }

        // Filtracja sygna³u dŸwiêkowego
        vector<double> filtered_y = filtr(y, czestotliwosc_graniczna, probkowanie, rzadFiltru);

        figure();
        plot(t, y, "-"); // rysowanie wykresu próbek dŸwiêkowych
        title("Wykres próbek dŸwiêkowych (oryginalny)");
        show();

        figure();
        plot(t, filtered_y, "-r"); // rysowanie wykresu przefiltrowanych próbek dŸwiêkowych
        title("Wykres próbek dŸwiêkowych (przefiltrowany)");
        show();
    }
}



PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");


    m.def("rysowanie_wykresow_podstawowych", &rysowanie_wykresow_podstawowych, "rysowanie_wykresow_podstawowych");
    m.def("rysowanie_wykresow_filtr", &rysowanie_wykresow_filtr, "rysowanie_wykresow_filtr");
    m.def("genAudio", &genAudio, "genAudio");


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}