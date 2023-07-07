#include "rk4.h"
#include <cassert>

template<class T>
RungeKutta4<T>::RungeKutta4(DataStruct<T>& _Un) : Un(_Un)
{
    nSteps = 2; // Reduce the number of steps to 2.
    currentStep = 0;

    coeffsA = new T[2];
    coeffsB = new T[2];
    coeffsA[0] = 0.;
    coeffsA[1] = 1.;
    coeffsB[0] = 0.5; // Updated coefficients for the new method.
    coeffsB[1] = 0.5;

    // Only one DataStruct required to store intermediate values "f".
    fi = new DataStruct<T>[1];
    fi[0].setSize(Un.getSize());

    Ui.setSize(Un.getSize());
};

template<class T>
RungeKutta4<T>::~RungeKutta4()
{
    delete[] fi;
    delete[] coeffsA;
    delete[] coeffsB;
};

template<class T>
int RungeKutta4<T>::getNumSteps()
{
    return nSteps;
};

template<class T>
void RungeKutta4<T>::initRK()
{
    currentStep = 0;
};

template<class T>
void RungeKutta4<T>::stepUi(T dt)
{
    assert(currentStep < nSteps);

    if (currentStep == 0)
    {
        T* dataUi = Ui.getData();
        const T* dataU = Un.getData();

        for (int n = 0; n < Ui.getSize(); n++)
        {
            dataUi[n] = dataU[n];
        }
    }
    else
    {
        T* datafi = fi[0].getData(); // Use the single DataStruct for "f".
        T* dataUi = Ui.getData();
        const T* dataU = Un.getData();

        for (int n = 0; n < Ui.getSize(); n++)
        {
            dataUi[n] = dataU[n] + coeffsA[currentStep] * dt * datafi[n];
        }
    }
};

template<class T>
void RungeKutta4<T>::finalizeRK(const T dt)
{
    T* dataUn = Un.getData();
    T* dataUi = Ui.getData();

    // set Ui to 0
    for (int n = 0; n < Ui.getSize(); n++)
    {
        dataUi[n] = 0.;
    }

    for (int s = 0; s < nSteps; s++)
    {
        const T* dataFi = fi[0].getData(); // Use the single DataStruct for "f".
        const T b = coeffsB[s];

        for (int n = 0; n < Ui.getSize(); n++)
        {
            dataUi[n] += b * dataFi[n];
        }
    }

    const T oneDiv6 = 1. / 6.;
    for (int n = 0; n < Ui.getSize(); n++)
    {
        dataUn[n] += dt * oneDiv6 * dataUi[n];
    }
};

template<class T>
void RungeKutta4<T>::setFi(DataStruct<T>& _F)
{
    T* dataFi = fi[0].getData(); // Use the single DataStruct for "f".
    const T* dataF = _F.getData();

    for (int n = 0; n < Ui.getSize(); n++)
    {
        dataFi[n] = dataF[n];
    }

    currentStep++;
};

template<class T>
DataStruct<T>* RungeKutta4<T>::currentU()
{
    return &Ui;
};

// Instantiate the template classes
template class RungeKutta4<float>;
template class RungeKutta4<double>;
