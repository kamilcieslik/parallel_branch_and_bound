#include <iostream>
#include "TimeMeasurement.h"
#include "TSPLIB_Parser.h"
#include "TravellingSalesmanProblem.h"
#include "AlgorithmTest.h"

void displayMenu(const std::string &info) //Menu dla problemu komiwojażera.
{
    std::cout << std::endl;
    std::cout << info << std::endl;
    std::cout << "1. Wczytaj z pliku." << std::endl;
    std::cout << "2. Wygeneruj miasta losowo." << std::endl;
    std::cout << "3. Wyświetl informacje o danych wejściowych problemu." << std::endl;
    std::cout << "4. Algorytm - zachłanny." << std::endl;
    std::cout << "5. Algorytm - metoda podziału i ograniczeń." << std::endl;
    std::cout << "0. Powrót do menu." << std::endl;
    std::cout << "Podaj opcje: ";
}

void menu_travelling_salesman_problem() //Obsługa problemu komiwojażera.
{
    TravellingSalesmanProblem s;
    std::string path;
    int option;
    do {
        displayMenu("*** Problem komiwojażera ***");
        std::cin >> option;
        std::cout << std::endl;
        switch (option) {
            case 1: //Tworzenie zbioru miast z pliku TSPLIB.
                std::cout << "Podaj sciezke pliku z danymi TSPLIB: ";
                std::cin >> path;
                try {
                    TSPLIB_Parser parser(path);
                    s.LoadArrayOfMatrixOfCities(parser.GetArrayOfMatrixCities(), parser.GetDimension(),
                                                parser.GetFileName(), parser.GetGraphType());
                } catch (std::logic_error &e) {
                    std::cout << e.what() << std::endl;
                } catch (std::invalid_argument &e_2) {
                    std::cout << e_2.what() << std::endl;
                }
                break;

            case 2: //Generowanie miast pseudolosowo.
                try {
                    s.GenerateRandomCities();
                }
                catch (std::invalid_argument &e) {
                    std::cout << e.what() << std::endl;
                }
                break;

            case 3: //Wyświetlanie zbioru miast.
                try {
                    s.PrintCitiesForTheTravellingSalesman(true);
                }
                catch (std::logic_error &e) {
                    std::cout << e.what() << std::endl;
                }
                break;

            case 4: //Algorytm 1. - zachłanny.
                try {
                    s.GreedyAlgorithm();
                    s.PrintCitiesForTheTravellingSalesman(false);
                    std::cout << std::endl;
                    s.PrintSolution();
                }
                catch (std::logic_error &e) {
                    std::cout << e.what() << std::endl;
                }
                break;

            case 5: //Algorytm 2. - metoda podziału i ograniczeń.
                try {
                    s.BranchAndBoundAlgorithm();
                    s.PrintCitiesForTheTravellingSalesman(false);
                    std::cout << std::endl;
                    s.PrintSolution();
                }
                catch (std::logic_error &e) {
                    std::cout << e.what() << std::endl;
                }
                break;

            default:
                break;
        }

    } while (option != 0);
}

void menu_tests() //Obsługa testów czasowych.
{
    AlgorithmTest test;
    int option;
    int numberOfRepetitions;
    int maxDistanceBetweenCity;
    do {
        std::cout << std::endl;
        std::cout << "*** Testy czasowe ***" << std::endl;
        std::cout << "1. Testy czasowe dla algorytmu zachłannego problemu komiwojażera." << std::endl;
        std::cout << "2. Testy czasowe dla algorytmu podziału i ograniczeń problemu komiwojażera." << std::endl;
        std::cout << "Podaj opcje: ";
        std::cin >> option;
        std::cout << std::endl;
        switch (option) {
            case 1: //Testy czasowe dla algorytmu zachłannego problemu komiwojażera.
                std::cout << "Podaj ilość instancji każdego zestawu danych w celu uśrednienia wyniku: ";
                std::cin >> numberOfRepetitions;
                std::cout << "Podaj maksymalną wartość odległości pomiędzy miastami: ";
                std::cin >> maxDistanceBetweenCity;
                test.TravellingSalesmanProblem_Test_Greedy(numberOfRepetitions, maxDistanceBetweenCity);
                break;

            case 2: //Testy czasowe dla algorytmu podziału i ograniczeń problemu komiwojażera.
                std::cout << "Podaj ilość instancji każdego zestawu danych w celu uśrednienia wyniku: ";
                std::cin >> numberOfRepetitions;
                std::cout << "Podaj maksymalną wartość odległości pomiędzy miastami: ";
                std::cin >> maxDistanceBetweenCity;
                test.TravellingSalesmanProblem_Test_BranchAndBound(numberOfRepetitions, maxDistanceBetweenCity);
                break;

            default:
                break;
        }

    } while (option != 0);
}

int main() {
    int option;
    do {
        std::cout << std::endl;
        std::cout << "==== MENU GŁÓWNE ===" << std::endl;
        std::cout << "1. Problem komiwojażera." << std::endl;
        std::cout << "2. Testy czasowe." << std::endl;
        std::cout << "0. Wyjście." << std::endl;
        std::cout << "Podaj opcje: ";
        std::cin >> option;
        std::cout << std::endl;

        switch (option) {
            case 1:
                menu_travelling_salesman_problem();
                break;
            case 2:
                menu_tests();
                break;
            default:
                break;
        }

    } while (option != 0);

    return 0;
}
