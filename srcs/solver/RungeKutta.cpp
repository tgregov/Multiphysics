#include "rungeKutta.hpp"

void RK1(double t, Field& field, PartialField& partialField, const Matrix& matrix,
         const Mesh& mesh, const SolverParams& solverParams, Field& temp, UsedF usedF)
{
    usedF(t, field, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
        field.u[unk] += field.DeltaU[unk]*solverParams.timeStep;
}


void RK2(double t, Field& field, PartialField& partialField, const Matrix& matrix,
         const Mesh& mesh, const SolverParams& solverParams, Field& temp, UsedF usedF)
{

    double h = solverParams.timeStep;

    usedF(t, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k1[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k1[unk]/2;
    }

    usedF(t + h/2, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k2[unk] = temp.DeltaU[unk]*h;
        field.u[unk] += field.k2[unk];
    }
}


void RK3(double t, Field& field, PartialField& partialField, const Matrix& matrix,
         const Mesh& mesh, const SolverParams& solverParams, Field& temp, UsedF usedF)
{
    double h = solverParams.timeStep;

    usedF(t, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k1[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k1[unk]/2;
    }

    usedF(t + h/2, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k2[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] - field.k1[unk] + 2*field.k2[unk];
    }

    usedF(t + h, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < 3 ; ++unk)
    {
        field.k3[unk] = temp.DeltaU[unk]*h;
        field.u[unk] += (field.k1[unk] + 4*field.k2[unk] + field.k3[unk])/6;
    }
}


void RK4(double t, Field& field, PartialField& partialField, const Matrix& matrix,
         const Mesh& mesh, const SolverParams& solverParams, Field& temp, UsedF usedF)
{
    double h = solverParams.timeStep;

    usedF(t, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k1[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k1[unk]/2;
    }

    usedF(t + h/2, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k2[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k2[unk]/2;
    }

    usedF(t + h/2, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k3[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k3[unk];
    }

    usedF(t + h, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k4[unk] = temp.DeltaU[unk]*h;
        field.u[unk] += (field.k1[unk] + 2*field.k2[unk] + 2*field.k3[unk] + field.k4[unk])/6;
    }
}
