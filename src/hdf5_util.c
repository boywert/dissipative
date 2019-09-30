/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/hdf5_util.c
 * \date        09/2014
 * \author      Federico Marinacci
 * \brief Contains the wrapper functions to the HDF5 library functions. 
 * \details The wrapper functions explicitly check for error conditions 
 *  and terminate the run if such conditions occur.  The HDF5 error handler 
 *  is disabled in case of termination not to repeat the error message of the 
 *  handler again at the program exit.
 * 
 * \par Major modifications and contributions:
 * 
 * - 11.09.2014 first file version
 */



#include "allvars.h"
#include "proto.h"

#ifdef HAVE_HDF5
#ifndef HDF5UTIL_H
#define HDF5UTIL_H

#include <hdf5.h>

/*!
 * This routine wraps creating a file to give a nice error message
 */
hid_t my_H5Fcreate(const char *fname, unsigned int flags, hid_t fcpl_id, hid_t fapl_id)
{
  hid_t file_id = H5Fcreate(fname, flags, fcpl_id, fapl_id);

#ifndef TOLERATE_WRITE_ERROR
  if(file_id < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to create file %s\n", ThisTask, fname);
    }
#endif

  return file_id;
}

/*!
 * This routine wraps creating a group to give a nice error message
 */
hid_t my_H5Gcreate(hid_t loc_id, const char *groupname, size_t size_hint)
{
  hid_t group_id = H5Gcreate(loc_id, groupname, size_hint);

#ifndef TOLERATE_WRITE_ERROR
  if(group_id < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to create group %s\n", ThisTask, groupname);
    }
#endif

  return group_id;
}

/*!
 * This routine wraps creating a dataset to give a nice error message
 */
hid_t my_H5Dcreate(hid_t loc_id, const char *datasetname, hid_t type_id, hid_t space_id, hid_t dcpl_id)
{
  hid_t dataset_id = H5Dcreate(loc_id, datasetname, type_id, space_id, dcpl_id);

#ifndef TOLERATE_WRITE_ERROR
  if(dataset_id < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, Error detected in HDF5: unable to create dataset %s\n", ThisTask, datasetname);
    }
#endif

  return dataset_id;
}

/*!
 * This routine wraps writing a dataset to give a nice error message
 */
herr_t my_H5Dwrite(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void *buf, const char *datasetname)
{
#ifdef TOLERATE_WRITE_ERROR
  if(WriteErrorFlag)
    return 0;
#endif

  herr_t status = H5Dwrite(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to write dataset %s\n", ThisTask, datasetname);
    }
#endif

  return status;
}

/*!
 * This routine wraps creating an attribute to give a nice error message
 */
hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id)
{
  hid_t attribute_id = H5Acreate(loc_id, attr_name, type_id, space_id, acpl_id);

#ifndef TOLERATE_WRITE_ERROR
  if(attribute_id < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to create attribute %s\n", ThisTask, attr_name);
    }
#endif

  return attribute_id;
}

/*!
 * This routine wraps writing an attribute to give a nice error message
 */
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name)
{
#ifdef TOLERATE_WRITE_ERROR
  if(WriteErrorFlag)
    return 0;
#endif

  herr_t status = H5Awrite(attr_id, mem_type_id, buf);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to write attribute %s\n", ThisTask, attr_name);
    }
#endif

  return status;
}

/*!
 * This routine wraps creating a dataspace to give a nice error message
 */
hid_t my_H5Screate(H5S_class_t type)
{
  hid_t dataspace_id = H5Screate(type);

#ifndef TOLERATE_WRITE_ERROR
  if(dataspace_id < 0)
    {
      H5Eset_auto(NULL, NULL);
      switch (type)
        {
        case H5S_SCALAR:
          terminate("On Task %d, error detected in HDF5: unable to create a scalar dataspace\n", ThisTask);
          break;
        case H5S_SIMPLE:
          terminate("On Task %d, error detected in HDF5: unable to create a simple dataspace\n", ThisTask);
          break;
        default:
          terminate("On Task %d, error detected in HDF5: unknown dataspace type\n", ThisTask);
          break;
        }
    }
#endif

  return dataspace_id;
}

/*!
 * This routine wraps creating a simple dataspace to give a nice error message
 */
hid_t my_H5Screate_simple(int rank, const hsize_t * current_dims, const hsize_t * maximum_dims)
{
  hid_t dataspace_id = H5Screate_simple(rank, current_dims, maximum_dims);

#ifndef TOLERATE_WRITE_ERROR
  if(dataspace_id < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to create a simple dataspace\n", ThisTask);
    }
#endif

  return dataspace_id;
}

/*!
 * This routine wraps opening a file to give a nice error message
 */
hid_t my_H5Fopen(const char *fname, unsigned int flags, hid_t fapl_id)
{
  hid_t file_id = H5Fopen(fname, flags, fapl_id);

  if(file_id < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to open file %s\n", ThisTask, fname);
    }

  return file_id;
}

/*!
 * This routine wraps opening a group to give a nice error message
 */
hid_t my_H5Gopen(hid_t loc_id, const char *groupname)
{
  hid_t group = H5Gopen(loc_id, groupname);

#ifndef TOLERATE_WRITE_ERROR
  if(group < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to open group %s\n", ThisTask, groupname);
    }
#endif

  return group;
}

/*!
 * This routine wraps opening a dataset to give a nice error message
 */
hid_t my_H5Dopen(hid_t file_id, const char *datasetname)
{
  hid_t dataset = H5Dopen(file_id, datasetname);

#ifndef TOLERATE_WRITE_ERROR
  if(dataset < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to open dataset %s\n", ThisTask, datasetname);
    }
#endif

  return dataset;
}

/*!
 * This routine wraps opening a dataset. However, in contrast to my_H5Dpoen(), if the dataset
 * does not exist it does not terminate the run. This is useful while reading an ICs file 
 * because in that case a non-exisitng dataset is put to zero (see also read_ic.c)  
 */
hid_t my_H5Dopen_if_existing(hid_t file_id, const char *datasetname)
{
  /* save error handler and disable it */
  H5E_auto_t errfunc;
  void *client_data;
  H5Eget_auto(&errfunc, &client_data);
  H5Eset_auto(NULL, NULL);

  hid_t dataset = H5Dopen(file_id, datasetname);

  /* reset error handler */
  H5Eset_auto(errfunc, client_data);

  return dataset;
}

/*!
 * This routine wraps opening an attribute to give a nice error message
 */
hid_t my_H5Aopen_name(hid_t loc_id, const char *attr_name)
{
  hid_t attribute_id = H5Aopen_name(loc_id, attr_name);

#ifndef TOLERATE_WRITE_ERROR
  if(attribute_id < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to open attribute %s\n", ThisTask, attr_name);
    }
#endif

  return attribute_id;
}

/*!
 * This routine wraps reading a dataset to give a nice error message
 */
herr_t my_H5Dread(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, void *buf, const char *datasetname)
{
  herr_t status = H5Dread(dataset_id, mem_type_id, mem_space_id, file_space_id, xfer_plist_id, buf);
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to read dataset %s\n", ThisTask, datasetname);
    }
  return status;
}

/*!
 * This routine wraps reading an attribute to give a nice error message
 */
hid_t my_H5Dget_space(hid_t dataset_id, const char *datasetname)
{
  hid_t status = H5Dget_space(dataset_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to determine space for dataset %s\n", ThisTask, datasetname);
    }
#endif

  return status;
}


/*!
 * This routine wraps reading an attribute to give a nice error message
 */
herr_t my_H5Aread(hid_t attr_id, hid_t mem_type_id, void *buf, const char *attr_name, hssize_t size)
{
  hid_t hdf5_space = H5Aget_space(attr_id);
  hssize_t attr_size = H5Sget_simple_extent_npoints(hdf5_space);
  H5Sclose(hdf5_space);

  if(attr_size != size)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: mismatch in size for attribute %s, expected size = %lld, actual attribute size = %lld\n", ThisTask, attr_name, size, attr_size);
    }

  herr_t status = H5Aread(attr_id, mem_type_id, buf);
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to read attribute %s\n", ThisTask, attr_name);
    }
  return status;
}

/*!
 * This routine wraps closing an attribute to give a nice error message
 */
herr_t my_H5Sset_extent_simple(hid_t space_id, int rank, const hsize_t * current_size, const hsize_t * maximum_size, const char *attr_name)
{
  herr_t status = H5Sset_extent_simple(space_id, rank, current_size, maximum_size);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to set extent for attribute %s\n", ThisTask, attr_name);
    }
#endif

  return status;
}

/*!
 * This routine wraps closing an attribute to give a nice error message
 */
herr_t my_H5Aclose(hid_t attr_id, const char *attr_name)
{
  herr_t status = H5Aclose(attr_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to close attribute %s\n", ThisTask, attr_name);
    }
#endif

  return status;
}

/*!
 * This routine wraps closing a dataset to give a nice error message
 */
herr_t my_H5Dclose(hid_t dataset_id, const char *datasetname)
{
  herr_t status = H5Dclose(dataset_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to close dataset %s\n", ThisTask, datasetname);
    }
#endif

  return status;
}

/*!
 * This routine wraps closing a group to give a nice error message
 */
herr_t my_H5Gclose(hid_t group_id, const char *groupname)
{
  herr_t status = H5Gclose(group_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to close group %s\n", ThisTask, groupname);
    }
#endif

  return status;
}

/*!
 * This routine wraps closing a file to give a nice error message
 */
herr_t my_H5Fclose(hid_t file_id, const char *fname)
{
  herr_t status = H5Fclose(file_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to close file %s\n", ThisTask, fname);
    }
#endif
  return status;
}

/*!
 * This routine wraps closing a dataspace to give a nice error message
 */
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type)
{
  herr_t status = H5Sclose(dataspace_id);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      switch (type)
        {
        case H5S_SCALAR:
          terminate("On Task %d, error detected in HDF5: unable to close a scalar dataspace\n", ThisTask);
          break;
        case H5S_SIMPLE:
          terminate("On Task %d, error detected in HDF5: unable to close a simple dataspace\n", ThisTask);
          break;
        default:
          terminate("On Task %d, error detected in HDF5: unknown dataspace type\n", ThisTask);
          break;
        }
    }
#endif

  return status;
}

/*!
 * This routine wraps copying an existing datatype to give a nice error message
 */
hid_t my_H5Tcopy(hid_t type_id)
{
  hid_t datatype_id = H5Tcopy(type_id);
#ifndef TOLERATE_WRITE_ERROR
  if(datatype_id < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not properly copy datatype\n", ThisTask);
    }
#endif
  return datatype_id;
}

/*!
 * This routine wraps closing a datatype to give a nice error message
 */
herr_t my_H5Tclose(hid_t type_id)
{
  herr_t status = H5Tclose(type_id);
#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not properly close datatype\n", ThisTask);
    }
#endif
  return status;
}

/*!
 * This routine wraps selecting a hyperslab to give a nice error message
 */
herr_t my_H5Sselect_hyperslab(hid_t space_id, H5S_seloper_t op, const hsize_t * start, const hsize_t * stride, const hsize_t * count, const hsize_t * block)
{
  herr_t status = H5Sselect_hyperslab(space_id, op, start, stride, count, block);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not properly select the chosen hyperslab\n", ThisTask);
    }
#endif
  return status;
}

/*!
 * This routine wraps returning the size in bytes of a given datatype to give a nice error message
 */
size_t my_H5Tget_size(hid_t datatype_id)
{
  size_t size = H5Tget_size(datatype_id);

#ifndef TOLERATE_WRITE_ERROR
  if(size == 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: unable to determine the size of the given datatype\n", ThisTask);
    }
#endif
  return size;
}

/*!
 * This routine wraps setting the size in bytes of a given datatype to give a nice error message
 */
herr_t my_H5Tset_size(hid_t datatype_id, size_t size)
{
  herr_t status = H5Tset_size(datatype_id, size);

#ifndef TOLERATE_WRITE_ERROR
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not properly set the size of the given datatype\n", ThisTask);
    }
#endif

  return status;
}

#ifdef HDF5_FILTERS
/*!
 * This routine wraps checking if all hdf5 filters selected for plist_id are available
 * to give a nice error message
 */
htri_t my_H5Pall_filters_avail(hid_t plist_id)
{
  htri_t status = H5Pall_filters_avail(plist_id);
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not properly verify the availability of all filters\n", ThisTask);
    }
  return status;
}

/*!
 * This routine wraps creating the property list of the given property class identified by class_id
 * to give a nice error message
 */
hid_t my_H5Pcreate(hid_t class_id)
{
  hid_t plist_id = H5Pcreate(class_id);
  if(plist_id < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not create the property list associated to the given property class\n", ThisTask);
    }
  return plist_id;
}

/*!
 * This routine wraps closing a property list to give a nice error message
 */
herr_t my_H5Pclose(hid_t plist)
{
  herr_t status = H5Pclose(plist);
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not close the input property list\n", ThisTask);
    }
  return status;
}

/*!
 * This routine wraps setting the size of the chunks of a chunked dataset
 * to give a nice error message
 */
herr_t my_H5Pset_chunk(hid_t plist, int ndims, const hsize_t * dim)
{
  herr_t status = H5Pset_chunk(plist, ndims, dim);
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not set chunk size for the dataset\n", ThisTask);
    }
  return status;
}

/*!
 * This routine wraps setting the use of the shuffle filter to give a nice 
 * error message
 */
herr_t my_H5Pset_shuffle(hid_t plist_id)
{
  herr_t status = H5Pset_shuffle(plist_id);
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not set the shuffle filter in the properties list\n", ThisTask);
    }
  return status;
}

/*!
 * This routine wraps setting the use of the deflate compression (gzip) to 
 * give a nice error message
 */
herr_t my_H5Pset_deflate(hid_t plist_id, uint level)
{
  herr_t status = H5Pset_deflate(plist_id, level);
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not set the deflate compression in the properties list\n", ThisTask);
    }
  return status;
}

/*!
 * This routine wraps setting the use of the Fletcher32 checksum 
 * to give a nice error message
 */
herr_t my_H5Pset_fletcher32(hid_t plist_id)
{
  herr_t status = H5Pset_fletcher32(plist_id);
  if(status < 0)
    {
      H5Eset_auto(NULL, NULL);
      terminate("On Task %d, error detected in HDF5: could not set the Fletcher32 checksum in the properties list\n", ThisTask);
    }
  return status;
}
#endif

#endif
#endif
