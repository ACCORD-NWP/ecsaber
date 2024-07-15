/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

 public:
  State(const Geometry & resol,
        const Model &,
        const eckit::Configuration & conf)
    : State(resol, conf) {}
  State(const Geometry & resol,
        const LinearModel &,
        const eckit::Configuration & conf)
    : State(resol, conf) {}
  State(const Geometry & resol,
        const Model &,
        const State & other)
    : State(resol, other) {}
  State(const Geometry & resol,
        const Model &,
        const State & other,
        const eckit::Configuration &)
    : State(resol, other) {}

  void interpolate(const Locations & locs,
                   GeoVaLs & gv) const
    {fields_->interpolate(locs, gv);}

  void forceWith(const State & other,
                 const Variables & vars)
    {fields_->forceWith(*(other.fields_), vars);}

  const atlas::FieldSet & fieldSet() const
    {return fields_->fieldSet();}
  atlas::FieldSet & fieldSet()
    {return fields_->fieldSet();}
  void synchronizeFields()
    {fields_->synchronizeFields();}

