 public:
  Increment(const Geometry &,
            const Variables &,
            const util::DateTime &,
            const util::DateTime &);


  void interpolateTL(const Locations & locs,
                     GeoVaLs & gv) const
    {fields_->interpolate(locs, gv);}
  void interpolateAD(const Locations & locs,
                     const GeoVaLs & gv)
    {fields_->interpolateAD(locs, gv);}

  double max(const Variables & var) const
    {return fields_->max(var);}
  double min(const Variables & var) const
    {return fields_->min(var);}
  util::DateTime & validTime()
    {return fields_->time();}

  friend eckit::Stream & operator<<(eckit::Stream &,
                                    const Increment &);
  friend eckit::Stream & operator>>(eckit::Stream &,
                                    Increment &);

  const atlas::FieldSet & fieldSet() const
    {return fields_->fieldSet();}
  atlas::FieldSet & fieldSet()
    {return fields_->fieldSet();}
  void synchronizeFields()
    {fields_->synchronizeFields();}
