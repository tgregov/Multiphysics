#ifndef shallow_writer_hpp_included
#define shallow_writer_hpp_included

#include <string>
#include "../../solver/field.hpp"
#include "../../params/Params.hpp"


/**
 * \brief Write data for shallow waters. You can write H, u, v, 0.*(u²+v²)
 * or the velocity field (boolean in whatToWrite).
 * \param uDisplay Vector (per element) of vector (per nodes) to write the data.
 * \param elementNumNodes Vector containg the number of nodes per element.
 * \param elementTags Vector containing the tag of all elements.
 * \param modelName Name of the model.
 * \param nbreStep Current time step.
 * \param t Current simulation physical time.
 * \param field Structure that contains all the main variables.
 * \param fluxCoeffs Coefficient of the physical flux.
 * \param whatToWrite Vector containing boolean describing which unknown to write.
 * \param viewTags Vector containing the rag of the different writing data's.
 */
void writeShallow(std::vector<std::vector<double>>& uDisplay,
                  const std::vector<unsigned int>& elementNumNodes,
                  const std::vector<int>& elementTags, const std::string& modelName,
                  unsigned int nbreStep, double t, const Field& field,
                  const std::vector<double>& fluxCoeffs,
                  const std::vector<bool>& whatToWrite, std::vector<int>& viewTags);

#endif /* shallow_writer_hpp_included */
