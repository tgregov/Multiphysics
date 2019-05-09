#include <iostream> //Without this missing #define makes
                    //the compiler crashes in mpi.h
#include <mpi.h>
#include "sendReceive.hpp"

void exchangeFlux(const Field& field, CompleteField& compField,
                  const DomainDiv& domainDiv, unsigned int rank,
                  const SolverParams& solverParams, const Mesh& mesh)
{
    for(unsigned short dim = 0 ; dim < mesh.dim ; ++dim)
    {
        for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
        {
            MPI_Allgatherv(&field.flux[dim][unk][0], domainDiv.node[rank], MPI_DOUBLE,
                            &compField.flux[dim][unk][0], domainDiv.node.data(),
                            domainDiv.nodePrec.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        }
    }
}

void exchangeUnk(const Field& field, CompleteField& compField,
                  const DomainDiv& domainDiv, unsigned int rank,
                  const SolverParams& solverParams, const Mesh& mesh)
{
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        MPI_Allgatherv(&field.u[unk][0], domainDiv.node[rank], MPI_DOUBLE,
                    &compField.u[unk][0], domainDiv.node.data(),
                    domainDiv.nodePrec.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    }
}
