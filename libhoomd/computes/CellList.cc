/*
Highly Optimized Object-oriented Many-particle Dynamics -- Blue Edition
(HOOMD-blue) Open Source Software License Copyright 2008, 2009 Ames Laboratory
Iowa State University and The Regents of the University of Michigan All rights
reserved.

HOOMD-blue may contain modifications ("Contributions") provided, and to which
copyright is held, by various Contributors who have granted The Regents of the
University of Michigan the right to modify and/or distribute such Contributions.

Redistribution and use of HOOMD-blue, in source and binary forms, with or
without modification, are permitted, provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions, and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of HOOMD-blue's
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR
ANY WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// $Id$
// $URL$
// Maintainer: joaander

#include <limits.h>
#include <boost/bind.hpp>
#include <boost/python.hpp>
#include <algorithm>

#include "CellList.h"

using namespace boost;
using namespace boost::python;
using namespace std;

/*! \param sysdef system to compute the cell list of
*/
CellList::CellList(boost::shared_ptr<SystemDefinition> sysdef)
    : Compute(sysdef),  m_nominal_width(Scalar(1.0f)), m_radius(1), m_max_cells(UINT_MAX), m_compute_tdb(false),
      m_flag_charge(false)
    {
    // allocation is deferred until the first compute() call - initialize values to dummy variables
    m_width = make_scalar3(0.0, 0.0, 0.0);
    m_dim = make_uint3(0,0,0);
    m_Nmax = 0;
    m_params_changed = true;
    m_particles_sorted = false;
    m_box_changed = false;
    GPUArray<unsigned int> conditions(3, exec_conf);
    m_conditions.swap(conditions);
    
    m_sort_connection = m_pdata->connectParticleSort(bind(&CellList::slotParticlesSorted, this));
    m_boxchange_connection = m_pdata->connectBoxChange(bind(&CellList::slotBoxChanged, this));
    }

CellList::~CellList()
    {
    m_sort_connection.disconnect();
    m_boxchange_connection.disconnect();
    }

/*! \returns Cell dimensions that match with the current width, box dimension, and max_cells setting
*/
uint3 CellList::computeDimensions()
    {
    uint3 dim;
    
    // calculate the bin dimensions
    const BoxDim& box = m_pdata->getBox();
    assert(box.xhi > box.xlo && box.yhi > box.ylo && box.zhi > box.zlo);
    dim.x = (unsigned int)((box.xhi - box.xlo) / (m_nominal_width));
    dim.y = (unsigned int)((box.yhi - box.ylo) / (m_nominal_width));
    
    if (m_sysdef->getNDimensions() == 2)
        {
        dim.z = 3;
    
        // decrease the number of bins if it exceeds the max
        if (dim.x * dim.y * dim.z > m_max_cells)
            {
            float scale_factor = powf(float(m_max_cells) / float(dim.x*dim.y*dim.z), 1.0f/2.0f);
            dim.x = int(float(dim.x)*scale_factor);
            dim.y = int(float(dim.y)*scale_factor);
            }
        }
    else
        {
        dim.z = (unsigned int)((box.zhi - box.zlo) / (m_nominal_width));

        // decrease the number of bins if it exceeds the max
        if (dim.x * dim.y * dim.z > m_max_cells)
            {
            float scale_factor = powf(float(m_max_cells) / float(dim.x*dim.y*dim.z), 1.0f/3.0f);
            dim.x = int(float(dim.x)*scale_factor);
            dim.y = int(float(dim.y)*scale_factor);
            dim.z = int(float(dim.z)*scale_factor);
            }
        }
    
    return dim;
    }

void CellList::compute(unsigned int timestep)
    {
    bool force = false;
    
    if (m_prof)
        m_prof->push("Cell");
    
    if (m_params_changed)
        {
        // need to fully reinitialize on any parameter change
        initializeAll();
        m_params_changed = false;
        force = true;
        }
    
    if (m_box_changed)
        {
        uint3 new_dim = computeDimensions();
        if (new_dim.x == m_dim.x && new_dim.y == m_dim.y && new_dim.z == m_dim.z)
            {
            // number of bins has not changed, only need to update width
            initializeWidth();
            }
        else
            {
            // number of bins has changed, need to fully reinitialize memory
            initializeAll();
            }
        
        m_box_changed = false;
        force = true;
        }
    
    if (m_particles_sorted)
        {
        // sorted particles simply need a forced update to get the proper indices in the data structure
        m_particles_sorted = false;
        force = true;
        }
    
    // only update if we need to
    if (shouldCompute(timestep) || force)
        {
        bool overflowed = false;
        do
            {
            computeCellList();
            
            overflowed = checkConditions();
            // if we overflowed, need to reallocate memory and reset the conditions
            if (overflowed)
                {
                initializeAll();
                cout << "Notice: cell list overflow, allocating " << m_Nmax << " slots per cell" << endl;
                resetConditions();
                }
            } while (overflowed);
        }
    
    if (m_prof)
        m_prof->pop();
    }

/*! \param num_iters Number of iterations to average for the benchmark
    \returns Milliseconds of execution time per calculation

    Calls computeCellList repeatedly to benchmark the compute
*/
double CellList::benchmark(unsigned int num_iters)
    {
    ClockSource t;
    
    // ensure that any changed parameters have been propagaged and memory allocated
    compute(0);
    
    // warm up run
    computeCellList();
    
#ifdef ENABLE_CUDA
    if (exec_conf->isCUDAEnabled())
        {
        cudaThreadSynchronize();
        CHECK_CUDA_ERROR();
        }
#endif
    
    // benchmark
    uint64_t start_time = t.getTime();
    for (unsigned int i = 0; i < num_iters; i++)
        computeCellList();
        
#ifdef ENABLE_CUDA
    if (exec_conf->isCUDAEnabled())
        cudaThreadSynchronize();
#endif
    uint64_t total_time_ns = t.getTime() - start_time;
    
    // convert the run time to milliseconds
    return double(total_time_ns) / 1e6 / double(num_iters);
    }

void CellList::initializeAll()
    {
    initializeWidth();
    initializeMemory();
    }

void CellList::initializeWidth()
    {
    if (m_prof)
        m_prof->push("init");
    
    // initialize dimensions and width
    m_dim = computeDimensions();

    const BoxDim& box = m_pdata->getBox();
    m_width.x = (box.xhi - box.xlo) / Scalar(m_dim.x);
    m_width.y = (box.yhi - box.ylo) / Scalar(m_dim.y);
    m_width.z = (box.zhi - box.zlo) / Scalar(m_dim.z);

    if (m_prof)
        m_prof->pop();

    }

void CellList::initializeMemory()
    {
    if (m_prof)
        m_prof->push("init");
    
    // if it is still set at 0, estimate Nmax
    if (m_Nmax == 0)
        {
        unsigned int estim_Nmax = (unsigned int)(ceilf(float(m_pdata->getN()*1.0f / float(m_dim.x*m_dim.y*m_dim.z))));
        m_Nmax = estim_Nmax + 8 - (estim_Nmax & 7);
        }
    else
        {
        // otherwise, round up to the nearest multiple of 8 if we are not already on one
        if ((m_Nmax & 7) != 0)
            m_Nmax = m_Nmax + 8 - (m_Nmax & 7);
        }
    
    // initialize indexers
    m_cell_indexer = Index3D(m_dim.x, m_dim.y, m_dim.z);
    m_cell_list_indexer = Index2D(m_Nmax, m_cell_indexer.getNumElements());
    m_cell_adj_indexer = Index2D((m_radius*2+1)*(m_radius*2+1)*(m_radius*2+1), m_cell_indexer.getNumElements());
    
    // allocate memory
    GPUArray<unsigned int> cell_size(m_cell_indexer.getNumElements(), exec_conf);
    m_cell_size.swap(cell_size);

    GPUArray<unsigned int> cell_adj(m_cell_adj_indexer.getNumElements(), exec_conf);
    m_cell_adj.swap(cell_adj);
    
    GPUArray<Scalar4> xyzf(m_cell_list_indexer.getNumElements(), exec_conf);
    m_xyzf.swap(xyzf);
    
    if (m_compute_tdb)
        {
        GPUArray<Scalar4> tdb(m_cell_list_indexer.getNumElements(), exec_conf);
        m_tdb.swap(tdb);
        }
    else
        {
        // array is no longer needed, discard it
        GPUArray<Scalar4> tdb;
        m_tdb.swap(tdb);
        }

    if (m_prof)
        m_prof->pop();
    
    initializeCellAdj();
    }

void CellList::initializeCellAdj()
    {
    if (m_prof)
        m_prof->push("init");
    
    ArrayHandle<unsigned int> h_cell_adj(m_cell_adj, access_location::host, access_mode::overwrite);
    
    // loop over all cells
    for (int k = 0; k < int(m_dim.z); k++)
        for (int j = 0; j < int(m_dim.y); j++)
            for (int i = 0; i < int(m_dim.x); i++)
                {
                unsigned int cur_cell = m_cell_indexer(i,j,k);
                unsigned int offset = 0;
                
                // loop over neighboring cells
                // need signed integer values for performing index calculations with negative values
                int r = int(m_radius);
                int mx = int(m_dim.x);
                int my = int(m_dim.y);
                int mz = int(m_dim.z);
                for (int nk = k-r; nk <= k+r; nk++)
                    for (int nj = j-r; nj <= j+r; nj++)
                        for (int ni = i-r; ni <= i+r; ni++)
                            {
                            int wrapi = ni % mx;
                            if (wrapi < 0)
                                wrapi += mx;
                            
                            int wrapj = nj % my;
                            if (wrapj < 0)
                                wrapj += my;
                            
                            int wrapk = nk % mz;
                            if (wrapk < 0)
                                wrapk += mz;
                            
                            unsigned int neigh_cell = m_cell_indexer(wrapi, wrapj, wrapk);
                            h_cell_adj.data[m_cell_adj_indexer(offset, cur_cell)] = neigh_cell;
                            offset++;
                            }
                
                // sort the adj list for each cell
                sort(&h_cell_adj.data[m_cell_adj_indexer(0, cur_cell)],
                     &h_cell_adj.data[m_cell_adj_indexer(offset, cur_cell)]);
                }
    
    if (m_prof)
        m_prof->pop();
    }

void CellList::computeCellList()
    {
    if (m_prof)
        m_prof->push("compute");
    
    // precompute scale factor
    Scalar3 scale = make_scalar3(Scalar(1.0) / m_width.x,
                                 Scalar(1.0) / m_width.y,
                                 Scalar(1.0) / m_width.z);
    
    // acquire the particle data
    ParticleDataArraysConst arrays = m_pdata->acquireReadOnly();
    const BoxDim& box = m_pdata->getBox();
    
    // access the cell list data arrays
    ArrayHandle<unsigned int> h_cell_size(m_cell_size, access_location::host, access_mode::overwrite);
    ArrayHandle<Scalar4> h_xyzf(m_xyzf, access_location::host, access_mode::overwrite);
    ArrayHandle<Scalar4> h_tdb(m_tdb, access_location::host, access_mode::overwrite);
    ArrayHandle<unsigned int> h_conditions(m_conditions, access_location::host, access_mode::readwrite);
    
    // shorthand copies of the indexers
    Index3D ci = m_cell_indexer;
    Index2D cli = m_cell_list_indexer;
    
    // clear the bin sizes to 0
    memset(h_cell_size.data, 0, sizeof(unsigned int) * m_cell_indexer.getNumElements());
    
    // for each particle
    for (unsigned int n = 0; n < arrays.nparticles; n++)
        {
        if (isnan(arrays.x[n]) || isnan(arrays.y[n]) || isnan(arrays.z[n]))
            {
            h_conditions.data[1] = n+1;
            continue;
            }
            
        // find the bin each particle belongs in
        unsigned int ib = (unsigned int)((arrays.x[n]-box.xlo)*scale.x);
        unsigned int jb = (unsigned int)((arrays.y[n]-box.ylo)*scale.y);
        unsigned int kb = (unsigned int)((arrays.z[n]-box.zlo)*scale.z);
        
        // need to handle the case where the particle is exactly at the box hi
        if (ib == m_dim.x)
            ib = 0;
        if (jb == m_dim.y)
            jb = 0;
        if (kb == m_dim.z)
            kb = 0;
            
        // sanity check
        assert(ib < (unsigned int)(m_dim.x) && jb < (unsigned int)(m_dim.y) && kb < (unsigned int)(m_dim.z));
        
        // record its bin
        unsigned int bin = ci(ib, jb, kb);
        // check if the particle is inside the dimensions
        if (bin >= ci.getNumElements())
            {
            h_conditions.data[2] = n+1;
            continue;
            }

        // setup the flag value to store
        Scalar flag;
        if (m_flag_charge)
            flag = arrays.charge[n];
        else
            flag = __int_as_scalar(n);

        // store the bin entries
        unsigned int offset = h_cell_size.data[bin];
        
        if (offset < m_Nmax)
            {
            h_xyzf.data[cli(offset, bin)] = make_scalar4(arrays.x[n], arrays.y[n], arrays.z[n], flag);
            if (m_compute_tdb)
                {
                h_tdb.data[cli(offset, bin)] = make_scalar4(__int_as_scalar(arrays.type[n]),
                                                            arrays.diameter[n],
                                                            __int_as_scalar(arrays.body[n]),
                                                            Scalar(0.0));
                }
            }
        else
            {
            h_conditions.data[0] = max(h_conditions.data[0], offset+1);
            }
        
        // increment the cell occupancy counter
        h_cell_size.data[bin]++;
        }
        
    m_pdata->release();
    
    if (m_prof)
        m_prof->pop();
    }
                
bool CellList::checkConditions()            
    {
    bool result = false;

    ArrayHandle<unsigned int> h_conditions(m_conditions, access_location::host, access_mode::read);

    // up m_Nmax to the overflow value, reallocate memory and set the overflow condition
    if (h_conditions.data[0] > m_Nmax)
        {
        m_Nmax = h_conditions.data[0];
        result = true;
        }                 

    // detect nan position errors
    if (h_conditions.data[1])
        {
        unsigned int n = h_conditions.data[1] - 1;
        cerr << endl << "***Error! Particle " << n << " has NaN for its position." << endl << endl;
        throw runtime_error("Error computing cell list");
        }

    // detect particles leaving box errors
    if (h_conditions.data[2])
        {
        unsigned int n = h_conditions.data[2] - 1;
        cerr << endl << "***Error! Particle " << n << " is no longer in the simulation box." << endl
             << endl;
        throw runtime_error("Error computing cell list");
        }

    return result;
    }

void CellList::resetConditions()
    {
    ArrayHandle<unsigned int> h_conditions(m_conditions, access_location::host, access_mode::overwrite);
    h_conditions.data[0] = 0;
    h_conditions.data[1] = 0;
    h_conditions.data[2] = 0;
    }



void export_CellList()
    {
    class_<CellList, boost::shared_ptr<CellList>, bases<Compute>, boost::noncopyable >
        ("CellList", init< boost::shared_ptr<SystemDefinition> >())
        .def("setNominalWidth", &CellList::setNominalWidth)
        .def("setRadius", &CellList::setRadius)
        .def("setMaxCells", &CellList::setMaxCells)
        .def("setComputeTDB", &CellList::setComputeTDB)
        .def("setFlagCharge", &CellList::setFlagCharge)
        .def("setFlagIndex", &CellList::setFlagIndex)
        .def("getWidth", &CellList::getWidth, return_internal_reference<>())
        .def("getDim", &CellList::getDim, return_internal_reference<>())
        .def("getNmax", &CellList::getNmax)
        .def("benchmark", &CellList::benchmark)
        ;
    }
