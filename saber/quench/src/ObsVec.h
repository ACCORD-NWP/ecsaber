/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "src/ObsSpace.h"

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace quench {

// -----------------------------------------------------------------------------
/// ObsVec class

class ObsVec : public util::Printable,
               private util::ObjectCounter<ObsVec> {
 public:
  static const std::string classname()
    {return "quench::ObsVec";}

  explicit ObsVec(const ObsSpace & obsSpace);
  ObsVec(const ObsVec &,
         const bool copy = true);
  ~ObsVec()
    {}

  ObsVec & operator= (const ObsVec &);
  ObsVec & operator*= (const double &);
  ObsVec & operator+= (const ObsVec &);
  ObsVec & operator-= (const ObsVec &);
  ObsVec & operator*= (const ObsVec &);
  ObsVec & operator/= (const ObsVec &);

  void zero();
  void axpy(const double &,
            const ObsVec &);
  void invert();
  void random();
  double dot_product_with(const ObsVec &) const;
  double rms() const;

  size_t size() const
    {return obsSpace_.size();}
  double & operator() (const size_t ii)
    {return data_[ii];}
  const double & operator() (const size_t ii) const
    {return data_[ii];}

  void read(const std::string & name)
    {obsSpace_.getdb(name, data_);}
  void save(const std::string & name) const
    {obsSpace_.putdb(name, data_);}

 private:
  void print(std::ostream &) const;

  const eckit::mpi::Comm & comm_;
  const ObsSpace & obsSpace_;
  std::vector<double> data_;
};

//-----------------------------------------------------------------------------

}  // namespace quench
