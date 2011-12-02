/**
 *  \file RMF/HDF5DataSetD.h
 *  \brief Handle read/write of Model data from/to files.
 *
 *  Copyright 2007-2011 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPLIBRMF_HDF_5DATA_SET_D_H
#define IMPLIBRMF_HDF_5DATA_SET_D_H

#include "RMF_config.h"
#include "hdf5_types.h"
#include "hdf5_handle.h"
#include "infrastructure_macros.h"
#include <boost/scoped_ptr.hpp>
#include <boost/filesystem/operations.hpp>
#include <algorithm>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/functional/hash.hpp>


namespace RMF {

  class HDF5Group;

  /** Store an index into a data set. Typedefs are provides
      for 1,2 and 3 dimension indexes, name like
      HDF5DataSetIndex2D.
   */
  template <int D>
  class HDF5DataSetIndexD
  {
    hsize_t d_[D];
    int compare(const HDF5DataSetIndexD<D> &o) const {
      for (unsigned int i=0; i< D; ++i) {
        if (d_[i] < o.d_[i]) return -1;
        else if (d_[i] > o.d_[i]) return 1;
      }
      return 0;
    }
  public:
    HDF5DataSetIndexD(const Ints &o) {
      IMP_RMF_USAGE_CHECK(o.size()==D, "Wrong number of values");
      std::copy(o.begin(), o.end(), d_);
    }
    HDF5DataSetIndexD() {
      std::fill(d_, d_+D, -1);
    }
    HDF5DataSetIndexD(unsigned int i) {
      IMP_RMF_USAGE_CHECK(D==1, "Constructor does not match dimension.");
      d_[0]=i;
    }
    HDF5DataSetIndexD(unsigned int i, unsigned int j) {
      IMP_RMF_USAGE_CHECK(D==2, "Constructor does not match dimension.");
      d_[0]=i;
      if (D>1) d_[1]=j;
    }
    HDF5DataSetIndexD(unsigned int i, unsigned int j, unsigned int k) {
      IMP_RMF_USAGE_CHECK(D==3, "Constructor does not match dimension.");
      d_[0]=i;
      // for clang
      if (D>1) d_[1]=j;
      if (D>2) d_[2]=k;
    }
#ifndef SWIG
    hsize_t& operator[](unsigned int i) {
      IMP_RMF_USAGE_CHECK(i < D, "Index out of range: "<< i);
      return d_[i];
    }
    hsize_t operator[](unsigned int i) const {
      IMP_RMF_USAGE_CHECK(i < D, "Index out of range: "<< i);
      return d_[i];
    }
    typedef const hsize_t * const_iterator;
    const_iterator begin() const {
      return d_;
    }
    const_iterator end() const {
      return d_+D;
    }
    typedef hsize_t * iterator;
    iterator begin() {
      return d_;
    }
    iterator end() {
      return d_+D;
    }
    hsize_t* get() const {
      return const_cast<hsize_t*>(d_);
    }
#endif
    int __getitem__(unsigned int i) const {
      if (i >= D) {
        IMP_RMF_THROW("Bad index " << i, std::runtime_error);
      }
      return operator[](i);
    }

    unsigned int get_dimension() const {return D;}
    void show(std::ostream &out) const {
      using std::operator<<;
      out << "(";
      for (unsigned int i=0; i< D; ++i) {
        if (i > 0) {
          out << ", ";
        }
        out << operator[](i);
      }
      out << ")";
    }
    IMP_RMF_COMPARISONS(HDF5DataSetIndexD);
    IMP_RMF_HASHABLE(HDF5DataSetIndexD,
                     size_t ret=0;
                     for (unsigned int i=0; i< D; ++i) {
                       boost::hash_combine(ret, static_cast<size_t>(d_[i]));
                     }
                     return ret;);
  };

#ifndef IMP_DOXYGEN
  typedef HDF5DataSetIndexD<1> HDF5DataSetIndex1D;
  typedef vector<HDF5DataSetIndex1D> HDF5DataSetIndex1Ds;
  typedef HDF5DataSetIndexD<2> HDF5DataSetIndex2D;
  typedef vector<HDF5DataSetIndex2D> HDF5DataSetIndex2Ds;
  typedef HDF5DataSetIndexD<3> HDF5DataSetIndex3D;
  typedef vector<HDF5DataSetIndex3D> HDF5DataSetIndex3Ds;
#endif

/** Data sets can be compressed using one of several algorithms.
 */
  enum Compression {GZIP_COMPRESSION, SLIB_COMPRESSION, NO_COMPRESSION};

  /** Wrap an HDF5 data set. Typedefs and python types are provided for
      data sets in 1,2, and 3 dimensions with all the
      \ref rmf_types "supported types". They are named as
      RMF::HDF5IndexDataSet2D (or RMF.HDF5IndexDataSet2).
   See
   \external{http://www.hdfgroup.org/HDF5/doc/UG/UG_frame10Datasets.html,
  the HDF5 manual} for more information.
  */
  template <class TypeTraits, unsigned int D>
  class HDF5DataSetD {
    struct Data {
      HDF5Handle h_;
      HDF5Handle ids_;
      HDF5Handle rds_;
      HDF5Handle sel_;
      hsize_t ones_[D];
    };

    boost::shared_ptr<Data> data_;

    int compare(const HDF5DataSetD<TypeTraits, D> &o) const {
      // not great, but...
      if (data_ && !o.data_) return -1;
      else if (o.data_ && !data_) return 1;
      else if (!o.data_ && !data_) return 0;
      else if (get_name() < o.get_name()) return -1;
      else if (get_name() > o.get_name()) return 1;
      else return 0;
    }

    bool get_is_null_value(const HDF5DataSetIndexD<D> &ijk) const {
      return TypeTraits::get_is_null_value(get_value(ijk));
    }
    const HDF5Handle& get_row_data_space() const {
      return data_->rds_;
    }
    const HDF5Handle& get_data_space() const {
      return data_->sel_;
    }
    void initialize_handles() {
      data_->sel_.open(H5Dget_space(data_->h_.get_hid()), &H5Sclose);
      // must be second
      hsize_t ret[D];
      std::fill(ret, ret+D, -1);
      IMP_HDF5_CALL(H5Sget_simple_extent_dims(get_data_space(),
                                              ret, NULL));
      IMP_RMF_INTERNAL_CHECK(ret[D-1] <1000000,
                             "extents not returned properly");
      if (ret[D-1] > 0) {
        // some versions will spew an error on this
        // we will call this function again before rds_ is needed
        //std::cout << "inializing row to " << ret[data_->dim_-1] << std::endl;
        data_->rds_.open(H5Screate_simple(1, ret+D-1,
                                          NULL), &H5Sclose);
      } else {
        //std::cout << "clearing row data" << std::endl;
        data_->rds_.close();
      }
    }
    void initialize() {
      hsize_t one=1;
      data_->ids_.open(H5Screate_simple(1, &one, NULL), &H5Sclose);
      std::fill(data_->ones_, data_->ones_+D, 1);
      //pos_.reset(new hsize_t[dim_]);
      //sel_= new HDF5SharedHandle(H5Dget_space(h_->get_hid()), &H5Sclose);
      initialize_handles();
    }
    friend class HDF5Group;
    HDF5DataSetD(HDF5SharedHandle* parent, std::string name,
                 Compression comp= NO_COMPRESSION):
      data_(new Data()) {
      //std::cout << "Creating data set " << name << std::endl;
      IMP_RMF_USAGE_CHECK(!H5Lexists(parent->get_hid(),
                                     name.c_str(), H5P_DEFAULT),
                          "Data set " << name << " already exists");
      hsize_t dims[D]={0};
      hsize_t cdims[D]={64};
      if (D >2) {
        std::fill(cdims+1, cdims+D-1, 2);
      }
      if (D >1) {
        cdims[D-1]=1;
      }
      hsize_t maxs[D];
      std::fill(maxs, maxs+D, H5S_UNLIMITED);
      IMP_HDF5_HANDLE(ds, H5Screate_simple(D, dims, maxs), &H5Sclose);
      IMP_HDF5_HANDLE(plist, H5Pcreate(H5P_DATASET_CREATE), &H5Pclose);
      IMP_HDF5_CALL(H5Pset_chunk(plist, D, cdims));
      IMP_HDF5_CALL(H5Pset_fill_value(plist, TypeTraits::get_hdf5_fill_type(),
                                      &TypeTraits::get_fill_value()));
      IMP_HDF5_CALL(H5Pset_fill_time(plist, H5D_FILL_TIME_IFSET));
      IMP_HDF5_CALL(H5Pset_alloc_time(plist, H5D_ALLOC_TIME_LATE));
      /*IMP_HDF5_CALL(H5Pset_szip (plist, H5_SZIP_NN_OPTION_MASK,
        32));*/
      if (comp==GZIP_COMPRESSION) {
        IMP_HDF5_CALL(H5Pset_deflate(plist, 9));
      } else if (comp == SLIB_COMPRESSION) {
        IMP_HDF5_CALL(H5Pset_szip (plist, H5_SZIP_NN_OPTION_MASK,
                                   32));
      }
      //std::cout << "creating..." << name << std::endl;
      data_->h_.open(H5Dcreate2(parent->get_hid(),
                               name.c_str(),
                               TypeTraits::get_hdf5_disk_type(),
                               ds, H5P_DEFAULT, plist, H5P_DEFAULT),
                     &H5Dclose);
      initialize();
      //std::cout << "done..." << std::endl;
    }
    // bool is to break symmetry
    HDF5DataSetD(HDF5SharedHandle* parent,
                 std::string name, bool): data_(new Data()) {
      IMP_RMF_USAGE_CHECK(H5Lexists(parent->get_hid(),
                                    name.c_str(), H5P_DEFAULT),
                          "Data set " << name << " does not exist");
      data_->h_.open(H5Dopen2(parent->get_hid(),
                             name.c_str(), H5P_DEFAULT),
                     &H5Dclose);
      //IMP_HDF5_HANDLE(s, H5Dget_space(h_->get_hid()), H5Sclose);
      IMP_HDF5_HANDLE(sel, H5Dget_space(data_->h_.get_hid()), &H5Sclose);
      IMP_RMF_USAGE_CHECK(H5Sget_simple_extent_ndims(sel)==D,
                          "Dimensions don't match. Got "
                          << H5Sget_simple_extent_ndims(sel)
                          << " but expected " << D);
      initialize();
    }
    void check_index(const HDF5DataSetIndexD<D> &ijk) const {
      HDF5DataSetIndexD<D> sz= get_size();
      for (unsigned int i=0; i< D; ++i) {
        IMP_RMF_USAGE_CHECK(ijk[i] < sz[i], "Index is out of range: "
                            << ijk[i] << " >= " << sz[i]);
      }
    }
  public:
#if !defined(SWIG) && !defined(IMP_DOXYGEN)
    HDF5DataSetD(hid_t file, std::string name): data_(new Data()) {
      IMP_RMF_USAGE_CHECK(H5Lexists(file,
                                    name.c_str(), H5P_DEFAULT),
                          "Data set " << name << " does not exist");
      data_->h_.open(H5Dopen2(file,
                             name.c_str(), H5P_DEFAULT),
                     &H5Dclose);
      //IMP_HDF5_HANDLE(s, H5Dget_space(h_->get_hid()), H5Sclose);
      IMP_HDF5_HANDLE(sel, H5Dget_space(data_->h_.get_hid()), &H5Sclose);
      IMP_RMF_USAGE_CHECK(H5Sget_simple_extent_ndims(sel)==D,
                          "Dimensions don't match. Got "
                          << H5Sget_simple_extent_ndims(sel)
                          << " but expected " << D);
      initialize();
    }
#endif
    typedef HDF5DataSetIndexD<D> Index;
    std::string get_name() const {
      char buf[10000];
      IMP_HDF5_CALL(H5Iget_name(data_->h_.get_hid(), buf, 10000));
      return std::string(buf);
    }
    void show(std::ostream &out) const {
      using std::operator<<;
      out << "HDF5DataSet " << get_name();
    }
    HDF5DataSetD(){}
    HDF5DataSetIndexD<D> get_size() const {
      //IMP_HDF5_HANDLE(s, H5Dget_space(h_->get_hid()), H5Sclose);
      HDF5DataSetIndexD<D> ret;
      IMP_HDF5_CALL(H5Sget_simple_extent_dims(get_data_space(),
                                              ret.begin(), NULL));
      return ret;
    }
    hid_t get_handle() const {
      return data_->h_.get_hid();
    }

    void set_value(const HDF5DataSetIndexD<D> &ijk,
                   typename TypeTraits::Type value) {
      IMP_RMF_IF_CHECK {
        check_index(ijk);
      }
      //IMP_HDF5_HANDLE(sel, H5Dget_space(h_->get_hid()), &H5Sclose);
      IMP_HDF5_CALL(H5Sselect_hyperslab(get_data_space(),
                                        H5S_SELECT_SET, ijk.get(),
                                        data_->ones_, data_->ones_,
                                        NULL));
      TypeTraits::write_value_dataset(data_->h_.get_hid(),
                                      data_->ids_.get_hid(),
                                      get_data_space(), value);
    }
    typename TypeTraits::Type get_value(const HDF5DataSetIndexD<D> &ijk) const {
      IMP_RMF_IF_CHECK {
        check_index(ijk);
      }
      //IMP_HDF5_HANDLE(sel, H5Dget_space(h_->get_hid()), &H5Sclose);
      IMP_HDF5_CALL(H5Sselect_hyperslab(get_data_space(),
                                        H5S_SELECT_SET, ijk.get(),
                                        data_->ones_, data_->ones_,
                                        NULL));
      return TypeTraits::read_value_dataset(data_->h_.get_hid(),
                                            data_->ids_.get_hid(),
                                            get_data_space());
    }
#ifndef SWIG
    typedef HDF5DataSetIndexD<D-1> RowIndex;
    void set_row( const RowIndex& ijkr,
                  const typename TypeTraits::Types& value) {
      HDF5DataSetIndexD<D> ijk;
      std::copy(ijkr.begin(), ijkr.end(), ijk.begin());
      ijk[D-1]=0;
      IMP_RMF_IF_CHECK {
        check_index(ijk);
      }
      hsize_t size[D]; std::fill(size, size+D-1, 1);
      size[D-1]= get_size()[D-1]; // set last to size of row
      //IMP_HDF5_HANDLE(sel, H5Dget_space(h_->get_hid()), &H5Sclose);
      IMP_HDF5_CALL(H5Sselect_hyperslab(get_data_space(),
                                        H5S_SELECT_SET, ijk.get(),
                                        data_->ones_, &size[0],
                                        NULL));
      TypeTraits::write_values_dataset(data_->h_.get_hid(),
                                       get_row_data_space().get_hid(),
                                       get_data_space(), value);
    }
    typename TypeTraits::Types get_row( const RowIndex ijkr) const {
      HDF5DataSetIndexD<D> ijk;
      std::copy(ijkr.begin(), ijkr.end(), ijk.begin());
      ijk[D-1]=0;
      IMP_RMF_IF_CHECK {
        check_index(ijk);
      }
      hsize_t size[D]; std::fill(size, size+D-1, 1);
      size[D-1]= get_size()[D-1]; // set last to size of row
      //IMP_HDF5_HANDLE(sel, H5Dget_space(h_->get_hid()), &H5Sclose);
      IMP_HDF5_CALL(H5Sselect_hyperslab(get_data_space(),
                                        H5S_SELECT_SET, ijk.get(),
                                        data_->ones_, &size[0],
                                        NULL));
      return TypeTraits::read_values_dataset(data_->h_.get_hid(),
                                             get_row_data_space().get_hid(),
                                             get_data_space(),
                                             size[D-1]);
    }
#endif
    //! Write a rectangular block starting at ln of size size
    void set_block(const Index&lb, const Index &size,
                   const typename TypeTraits::Types& value) {
      IMP_RMF_IF_CHECK {
        check_index(lb);
         Index last=lb;
        // offset size by one and check...
        unsigned int total=1;
        for (unsigned int i=0; i< D; ++i) {
          total*= size[i];
          last[i]+=size[i]-1;
        }
        IMP_RMF_USAGE_CHECK(total==value.size(),
                            "Block has size " << total << " but found "
                            << value.size() << " values");
        check_index(last);
      }
      //IMP_HDF5_HANDLE(sel, H5Dget_space(h_->get_hid()), &H5Sclose);
      IMP_HDF5_CALL(H5Sselect_hyperslab(get_data_space(),
                                        H5S_SELECT_SET, lb.get(),
                                        data_->ones_, size.get(),
                                        NULL));
      hsize_t sz= value.size();
      IMP_HDF5_HANDLE(input, H5Screate_simple(1, &sz,
                                        NULL), &H5Sclose);
      TypeTraits::write_values_dataset(data_->h_.get_hid(),
                                       input,
                                       get_data_space(), value);
    }
    //! Write a rectangular block starting at ln of size size
    typename TypeTraits::Types get_block( const Index &lb,
                                          const Index &size) const {
      hsize_t total=1;
      for (unsigned int i=0; i< D; ++i) {
        total*= size[i];
      }
      IMP_RMF_IF_CHECK {
        check_index(lb);
      }
      //IMP_HDF5_HANDLE(sel, H5Dget_space(h_->get_hid()), &H5Sclose);
      IMP_HDF5_CALL(H5Sselect_hyperslab(get_data_space(),
                                        H5S_SELECT_SET, lb.get(),
                                        data_->ones_, size.get(),
                                        NULL));
      IMP_HDF5_HANDLE(input, H5Screate_simple(1, &total,
                                        NULL), &H5Sclose);
      return TypeTraits::read_values_dataset(data_->h_.get_hid(),
                                             input,
                                             get_data_space(),
                                             total);
    }

    void set_size(const HDF5DataSetIndexD<D> &ijk) {
      hsize_t nd[D]; std::copy(ijk.begin(), ijk.end(), nd);;
      IMP_HDF5_CALL(H5Dset_extent(data_->h_.get_hid(),
                                  &nd[0]));
      initialize_handles();
    }
#if !defined(IMP_DOXYGEN) && !defined(SWIG)
    // replace with safe bool
    operator bool() const {
      return data_;
    }
    bool operator!() const {
      return !data_;
    }
#endif
    IMP_RMF_COMPARISONS(HDF5DataSetD);
  };

#ifndef IMP_DOXYGEN
  /** An HDF5 data set for integers with dimension 2. */
  template <class Traits, unsigned int D>
  class HDF5DataSetDTraits: public StringTraits {
    static HDF5DataSetD<Traits, D> get_data_set(hid_t dsc,
                                                std::string name) {
      if (name.empty()) {
        return HDF5DataSetD<Traits, D>();
      } else {
        hid_t file=H5Iget_file_id(dsc);
        return HDF5DataSetD<Traits, D>(file, name);
      }
    }
  public:
    typedef HDF5DataSetD<Traits, D> Type;
    typedef vector<Type> Types;
    static void write_value_dataset(hid_t d, hid_t is,
                                    hid_t s,
                                    Type v) {
      return StringTraits::write_value_dataset(d, is, s, v.get_name());
    }
    static Type read_value_dataset(hid_t d, hid_t is,
                                   hid_t sp) {
      String i= StringTraits::read_value_dataset(d, is, sp);
      if (!i.empty()) return get_data_set(d, i);
      else return Type();
    }
    static void write_values_dataset(hid_t d, hid_t is,
                                     hid_t s,
                                     const Types& v) {
      Strings vi(v.size());
      for (unsigned int i=0; i< v.size(); ++i) {
        vi[i]= v[i].get_name();
      }
      StringTraits::write_values_dataset(d, is, s, vi);
    }
    static Types read_values_dataset(hid_t d, hid_t is,
                                     hid_t sp,
                                     unsigned int sz)  {
      Strings reti= StringTraits::read_values_dataset(d, is, sp, sz);
      Types ret(reti.size());
      for (unsigned int i=0; i< ret.size(); ++i) {
        ret[i]= get_data_set(d, reti[i]);
      }
      return ret;
    }
    static Types
      read_values_attribute(hid_t a, unsigned int size) {
      Strings is= StringTraits::read_values_attribute(a, size);
      Types ret(is.size());
      for (unsigned int i=0; i< ret.size(); ++i) {
        ret[i]=get_data_set(a, is[i]);
      }
      return ret;
    }
    static void write_values_attribute(hid_t a,
                                       const Types &v) {
      Strings is(v.size());
      for (unsigned int i=0; i< v.size(); ++i) {
        is[i]=v[i].get_name();
      }
      StringTraits::write_values_attribute(a, is);
    }
    static std::string get_name() {
      std::ostringstream oss;
      using std::operator<<;
      oss << "data set " << Traits::get_name() << " " << D << "d";
      return oss.str();
    }
    static unsigned int get_index() {
      return 10+Traits::get_index()+10*D;
    }
    static Type get_null_value() {
      return Type();
    }
    static bool get_is_null_value(Type i) {
      return i== get_null_value();
    }
  };
#endif

#ifndef IMP_DOXYGEN
  /** \ingroup hdf5 */
  typedef HDF5DataSetDTraits<IndexTraits, 2> IndexDataSet2DTraits;
  /** \ingroup hdf5 */
  typedef HDF5DataSetDTraits<FloatTraits, 2> FloatDataSet2DTraits;

#define IMP_RMF_DECLARE_DATA_SET(lcname, Ucname, PassValue, ReturnValue, \
                                 PassValues, ReturnValues)              \
  typedef HDF5DataSetD<Ucname##Traits, 1> HDF5##Ucname##DataSet1D;      \
  typedef vector<HDF5##Ucname##DataSet1D> HDF5##Ucname##DataSet1Ds; \
  typedef HDF5DataSetD<Ucname##Traits, 2> HDF5##Ucname##DataSet2D;      \
  typedef vector<HDF5##Ucname##DataSet2D> HDF5##Ucname##DataSet2Ds; \
  typedef HDF5DataSetD<Ucname##Traits, 3> HDF5##Ucname##DataSet3D;      \
  typedef vector<HDF5##Ucname##DataSet3D> HDF5##Ucname##DataSet3Ds

  /** \name Basic data set types
       \ingroup hdf5
       @{
  */
  IMP_RMF_FOREACH_TYPE(IMP_RMF_DECLARE_DATA_SET);
  /** @} */
#endif
} /* namespace RMF */

#endif /* IMPLIBRMF_HDF_5DATA_SET_D_H */
