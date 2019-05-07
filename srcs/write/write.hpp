#ifndef write_hpp_included
#define write_hpp_included

#include "../solver/field.hpp"
#include "../params/Params.hpp"

void writeShallow(std::vector<std::vector<double>>& uDisplay,
                  const std::vector<unsigned int>& elementNumNodes,
                  const std::vector<int>& elementTags, const std::string& modelName,
                  unsigned int nbreStep, double t, const Field& field,
                  const std::vector<bool>& whatToWrite, std::vector<int>& viewTags);

void writeTransport(std::vector<std::vector<double>>& uDisplay,
                  const std::vector<unsigned int>& elementNumNodes,
                  const std::vector<int>& elementTags, const std::string& modelName,
                  unsigned int nbreStep, double t, const Field& field,
                  const std::vector<bool>& whatToWrite, std::vector<int>& viewTags);

void writeEnd(const std::vector<int>& viewTags, const std::vector<bool>& whatToWrite);

#endif /* write_hpp_included */
