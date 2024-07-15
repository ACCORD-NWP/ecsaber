/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

 public:
  Covariance(const Geometry &,
             const Variables &,
             const eckit::Configuration &,
             const State &)
    {}

  void linearize(const State &,
                 const Geometry &)
    {}

  void multiplySqrt(const IncrModCtlVec &,
                    Increment &) const
    {throw eckit::NotImplemented(Here());}
  void multiplySqrtTrans(const Increment &,
                         IncrModCtlVec &) const
    {throw eckit::NotImplemented(Here());}
