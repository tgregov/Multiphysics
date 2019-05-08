#include "RungeKutta.hpp"
#include "../mpi/sendReceive.hpp"
#include <iostream>
#include <mpi.h>

void RK1(double t, Field& field, PartialField& partialField,
         CompleteField& compField, const Matrix& matrix, const DomainDiv& domainDiv,
         unsigned int rank, const Mesh& mesh, const SolverParams& solverParams,
         Field& temp, CompleteField& tempCompField, UsedF usedF)
{
    usedF(t, field, compField, matrix, mesh, solverParams, domainDiv, rank);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
        field.u[unk] += field.DeltaU[unk]*solverParams.timeStep;

    solverParams.flux(field, partialField, solverParams, false);
    MPI_Barrier(MPI_COMM_WORLD);
    exchangeFlux(field, compField, domainDiv, rank, solverParams, mesh);
    exchangeUnk(field, compField, domainDiv, rank, solverParams, mesh);
    MPI_Barrier(MPI_COMM_WORLD);
}


void RK2(double t, Field& field, PartialField& partialField,
         CompleteField& compField, const Matrix& matrix, const DomainDiv& domainDiv,
         unsigned int rank, const Mesh& mesh, const SolverParams& solverParams,
         Field& temp, CompleteField& tempCompField, UsedF usedF)
{

    double h = solverParams.timeStep;

    usedF(t, temp, tempCompField, matrix, mesh, solverParams, domainDiv, rank);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k1[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k1[unk]/2;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    exchangeFlux(temp, tempCompField, domainDiv, rank, solverParams, mesh);
    exchangeUnk(temp, tempCompField, domainDiv, rank, solverParams, mesh);
    MPI_Barrier(MPI_COMM_WORLD);

    usedF(t + h/2, temp, tempCompField, matrix, mesh, solverParams, domainDiv, rank);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k2[unk] = temp.DeltaU[unk]*h;
        field.u[unk] += field.k2[unk];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    exchangeFlux(field, compField, domainDiv, rank, solverParams, mesh);
    exchangeUnk(field, compField, domainDiv, rank, solverParams, mesh);
    MPI_Barrier(MPI_COMM_WORLD);
}


void RK3(double t, Field& field, PartialField& partialField,
         CompleteField& compField, const Matrix& matrix, const DomainDiv& domainDiv,
         unsigned int rank, const Mesh& mesh, const SolverParams& solverParams,
         Field& temp, CompleteField& tempCompField, UsedF usedF)
{
    double h = solverParams.timeStep;

    usedF(t, temp, compField, matrix, mesh, solverParams, domainDiv, rank);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k1[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k1[unk]/2;
    }

    usedF(t + h/2, temp, compField, matrix, mesh, solverParams, domainDiv, rank);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k2[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] - field.k1[unk] + 2*field.k2[unk];
    }

    usedF(t + h, temp, compField, matrix, mesh, solverParams, domainDiv, rank);
    for(unsigned short unk = 0 ; unk < 3 ; ++unk)
    {
        field.k3[unk] = temp.DeltaU[unk]*h;
        field.u[unk] += (field.k1[unk] + 4*field.k2[unk] + field.k3[unk])/6;
    }
}


void RK4(double t, Field& field, PartialField& partialField,
         CompleteField& compField, const Matrix& matrix, const DomainDiv& domainDiv,
         unsigned int rank, const Mesh& mesh, const SolverParams& solverParams,
         Field& temp, CompleteField& tempCompField, UsedF usedF)
{
    double h = solverParams.timeStep;

    usedF(t, temp, compField, matrix, mesh, solverParams, domainDiv, rank);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k1[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k1[unk]/2;
    }

    usedF(t + h/2, temp, compField, matrix, mesh, solverParams, domainDiv, rank);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k2[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k2[unk]/2;
    }

    usedF(t + h/2, temp, compField, matrix, mesh, solverParams, domainDiv, rank);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k3[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k3[unk];
    }

    usedF(t + h, temp, compField, matrix, mesh, solverParams, domainDiv, rank);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k4[unk] = temp.DeltaU[unk]*h;
        field.u[unk] += (field.k1[unk] + 2*field.k2[unk] + 2*field.k3[unk] + field.k4[unk])/6;
    }
}
