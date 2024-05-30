/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <string>
#include <vector>

 public:
  const eckit::mpi::Comm & timeComm() const
    {return eckit::mpi::self();}
  const std::vector<double> & vert_coord_avg(const std::string & var) const
    {return groups_[groupIndex_.at(var)].vert_coord_avg_;}
