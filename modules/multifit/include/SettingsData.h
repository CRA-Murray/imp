/**
 *  \file SettingsData.h
 *  \brief stored multifit settings data
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPMULTIFIT_SETTINGS_DATA_H
#define IMPMULTIFIT_SETTINGS_DATA_H

#include <IMP/algebra/Transformation3D.h>
#include "multifit_config.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

IMPMULTIFIT_BEGIN_NAMESPACE

//! Holds data about a component needed for optimization
class IMPMULTIFITEXPORT ComponentHeader {
  public:
    ComponentHeader() {
      name_="";
      filename_="";
      pdb_ap_fn_="";
      num_ap_=0;
      transformations_fn_="";
      cmm_ap_fn_="";
      reference_fn_="";
    }
    void set_name(const std::string &name) {name_=name;}
    inline std::string get_name() const {return name_;}
    void set_filename(const std::string &filename){filename_=filename;}
    inline std::string get_filename() const {return filename_;}
    void set_pdb_ap_fn(const std::string &pdb_ap_fn) {pdb_ap_fn_=pdb_ap_fn;}
    inline std::string get_pdb_ap_fn() const {return pdb_ap_fn_;}
    void set_num_ap(int num_ap) {num_ap_=num_ap;}
    inline int get_num_ap() const {return num_ap_;}
    void set_transformations_fn(std::string transformations_fn)
     { transformations_fn_=transformations_fn;}
    std::string get_transformations_fn() const {return transformations_fn_;}
    void set_cmm_ap_fn(const std::string &cmm_fn) {cmm_ap_fn_=cmm_fn;}
    std::string get_cmm_ap_fn() const {return cmm_ap_fn_;}
    void set_reference_fn(const std::string &ref_fn){reference_fn_=ref_fn;}
    std::string get_reference_fn() const {return reference_fn_;}
  protected:
    std::string name_;
    std::string filename_;
    std::string pdb_ap_fn_;
    int num_ap_;
    std::string transformations_fn_;
    std::string cmm_ap_fn_;
    std::string reference_fn_;
  };
//! Holds data about the assembly density needed for optimization
class IMPMULTIFITEXPORT AssemblyHeader {
  public:
    AssemblyHeader(){
      dens_fn_="";
      resolution_=0.;
      spacing_=0.;
      pdb_fine_ap_fn_="";
      pdb_coarse_ap_fn_="";
      cmm_ap_fn_="";
    }
    void set_dens_fn(const std::string &dens_fn) {dens_fn_=dens_fn;}
    std::string get_dens_fn() const {return dens_fn_;}
    void set_resolution(float res) {resolution_=res;}
    float get_resolution() const {return resolution_;}
    float get_spacing() const {return spacing_;}
    void set_spacing(float spacing) {spacing_=spacing;}
    algebra::Vector3D get_origin() const {return origin_;}
    void set_origin(algebra::Vector3D origin) {origin_=origin;}
    std::string get_pdb_fine_ap_fn () const {return pdb_fine_ap_fn_;}
    void set_pdb_fine_ap_fn (const std::string &new_fn) {
      pdb_fine_ap_fn_=new_fn;}
    std::string get_pdb_coarse_ap_fn () const {return pdb_coarse_ap_fn_;}
    void set_pdb_coarse_ap_fn (const std::string &new_fn) {
     pdb_coarse_ap_fn_=new_fn;}
    std::string get_cmm_ap_fn () const {return cmm_ap_fn_;}
    void set_cmm_ap_fn (const std::string &new_fn) {cmm_ap_fn_=new_fn;}
  protected:
    std::string dens_fn_;
    float resolution_;
    float spacing_;
    algebra::Vector3D origin_;
    std::string pdb_fine_ap_fn_;
    std::string pdb_coarse_ap_fn_;
    std::string cmm_ap_fn_;
  };


//! Holds header data for optimization
class IMPMULTIFITEXPORT SettingsData {
public:
public:
  SettingsData(){}
  void add_component_header(ComponentHeader ch) {
    comp_data_.push_back(ch);
  }
  void set_assembly_header(AssemblyHeader h) {
    dens_data_=h;
  }
  ComponentHeader get_component_header(int i) const {
    IMP_USAGE_CHECK(i<comp_data_.size(),"index out of range\n");
    return comp_data_[i];
  }
  int get_number_of_component_headers() const {
    return comp_data_.size();
  }
  AssemblyHeader get_assembly_header() const {
    return dens_data_;
  }
protected:
  std::vector<ComponentHeader> comp_data_;
  AssemblyHeader dens_data_;
};

IMPMULTIFITEXPORT SettingsData read_settings(const char *filename);
IMPMULTIFIT_END_NAMESPACE
#endif /* IMPMULTIFIT_SETTINGS_DATA_H */
