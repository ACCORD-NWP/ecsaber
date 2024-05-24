/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_VARIABLES_H_
#define OOPS_BASE_VARIABLES_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/VariablesBase.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Class to set an container of variable names and manipulate it.
///        One option is either construct the variable object with meta data
///        or add metadata to the variable object. The metadata keys typically are
///        variable names but don't need to be consistent with the vector of variable
///        names.
///
///        Most operators / methods do not affect/interact with the meta data.
///        The exceptions are:
///         += JediVariables;  where it updates the metadata and appends extra metadata
///                        from the right hand side JediVariables object.
///         == JediVariables;  where it compares the metadata for metadata keys that
///                        are internally consistent with the variable names.
///         .addMetaData(  that will either update a value or add a value within
///                        some metadata for a metadata (variable name) key.
class JediVariables : public VariablesBase {
 public:
  static const std::string classname() {return "oops::JediVariables";}

  JediVariables() = default;
  JediVariables(const eckit::Configuration &, const std::string &);
  explicit JediVariables(const std::vector<std::string> &);
  JediVariables(const eckit::Configuration &, const std::vector<std::string> & vars);

  JediVariables & operator+=(const JediVariables &);
  JediVariables & operator-=(const JediVariables &);
  JediVariables & operator-=(const std::string &);

  bool operator==(const JediVariables &) const;
  bool operator!=(const JediVariables &) const;
  bool operator<=(const JediVariables &) const;

  /// make this JediVariables an intersection between this JediVariables and other variables
  void intersection(const JediVariables & other);

  const eckit::Configuration & variablesMetaData() const {return varMetaData_;}

  void sort();

  // Metadata
  bool hasMetaData(const std::string & varname,
                   const std::string & keyname) const;
  template<typename T>
  void addMetaData(const std::string & varname,
                   const std::string & keyname,
                   const T & keyvalue);
  template<typename T>
  T getMetaData(const std::string & varname,
                const std::string & keyname) const;
  // TODO(Later): might be replaced with getMetaData<int>(varname, "levels")
  int getLevels(const std::string &) const;

 private:
  void print(std::ostream &) const override;
  void setConf();

  template<typename T>
  void getVariableSubKeyValue(const std::string & varname,
                              const std::string & keyname,
                              const eckit::Configuration & conf,
                              T & keyvalue) const;

  template<typename T>
  void setVariableSubKeyValue(const std::string & varname,
                              const std::string & keyname,
                              const T & keyvalue,
                              eckit::LocalConfiguration & lconf);

  eckit::LocalConfiguration varMetaData_;
};

// -----------------------------------------------------------------------------

template<typename T>
void JediVariables::addMetaData(const std::string & varname,
                            const std::string & keyname,
                            const T & keyvalue) {
  setVariableSubKeyValue(varname, keyname, keyvalue, varMetaData_);
}

// -----------------------------------------------------------------------------

template<typename T>
T JediVariables::getMetaData(const std::string & varname,
                         const std::string & keyname) const {
  T keyvalue;
  getVariableSubKeyValue(varname, keyname, varMetaData_, keyvalue);
  return keyvalue;
}

// -----------------------------------------------------------------------------

template<typename T>
void JediVariables::getVariableSubKeyValue(const std::string & varname,
                                       const std::string & keyname,
                                       const eckit::Configuration & variablesconf,
                                       T & keyvalue) const {
  ASSERT(!variablesconf.empty());
  ASSERT(variablesconf.has(varname));
  ASSERT(variablesconf.getSubConfiguration(varname).has(keyname));
  variablesconf.getSubConfiguration(varname).get(keyname, keyvalue);
}

// -----------------------------------------------------------------------------

template<typename T>
void JediVariables::setVariableSubKeyValue(const std::string & varname,
                                       const std::string & keyname,
                                       const T & keyvalue,
                                       eckit::LocalConfiguration & variableslconf) {
  eckit::LocalConfiguration variablelconf =
    variableslconf.has(varname) ? variableslconf.getSubConfiguration(varname) :
                                  eckit::LocalConfiguration();
  variablelconf.set(keyname, keyvalue);
  variableslconf.set(varname, variablelconf);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_VARIABLES_H_
