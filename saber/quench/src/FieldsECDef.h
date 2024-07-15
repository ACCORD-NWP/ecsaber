/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <algorithm>
#include <vector>

 public:
  double min(const Variables &) const;
  double max(const Variables &) const;
  void interpolate(const Locations &,
                   GeoVaLs &) const;
  void interpolateAD(const Locations &,
                     const GeoVaLs &);
  void forceWith(const Fields &,
                 const Variables &);
  void synchronizeFields();
  friend eckit::Stream & operator<<(eckit::Stream &,
                                    const Fields &);
  friend eckit::Stream & operator>>(eckit::Stream &,
                                    Fields &);

 private:
  std::vector<quench::Interpolation>::iterator setupObsInterpolation(const Locations &) const;
  void reduceDuplicatePoints();
