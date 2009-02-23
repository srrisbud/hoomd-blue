/*
Highly Optimized Object-Oriented Molecular Dynamics (HOOMD) Open
Source Software License
Copyright (c) 2008 Ames Laboratory Iowa State University
All rights reserved.

Redistribution and use of HOOMD, in source and binary forms, with or
without modification, are permitted, provided that the following
conditions are met:

* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names HOOMD's
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND
CONTRIBUTORS ``AS IS''  AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS  BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

// $Id$
// $UR$

/*! \file GaussianForceGPU.cc
	\brief Defines the GaussianForceGPU class
*/

#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4103 4244 )
#endif

#include "GaussianForceGPU.h"
#include "cuda_runtime.h"

#include <stdexcept>

#include <boost/python.hpp>
using namespace boost::python;

#include <boost/bind.hpp>

using namespace boost;
using namespace std;

/*! \param pdata Particle Data to compute forces on
 	\param nlist Neighborlist to use for computing the forces
	\param r_cut Cuttoff radius beyond which the force is 0
	\post memory is allocated and all parameters epsilon and sigma are set to 0.0
*/
GaussianForceGPU::GaussianForceGPU(boost::shared_ptr<ParticleData> pdata, boost::shared_ptr<NeighborList> nlist, Scalar r_cut) 
	: GaussianForceCompute(pdata, nlist, r_cut)
	{
	// can't run on the GPU if there aren't any GPUs in the execution configuration
	if (exec_conf.gpu.size() == 0)
		{
		cerr << endl << "***Error! Creating a GaussianForceGPU with no GPU in the execution configuration" << endl << endl;
		throw std::runtime_error("Error initializing GaussianForceGPU");
		}
	
	if (m_ntypes > 44)
		{
		cerr << endl << "***Error! GaussianForceGPU cannot handle " << m_ntypes << " types" << endl << endl;
		throw runtime_error("Error initializing GaussianForceGPU");
		}
		
	// default block size is the highest performance in testing on different hardware
	// choose based on compute capability of the device
	cudaDeviceProp deviceProp;
	int dev;
	exec_conf.gpu[0]->call(bind(cudaGetDevice, &dev));	
	exec_conf.gpu[0]->call(bind(cudaGetDeviceProperties, &deviceProp, dev));
	if (deviceProp.major == 1 && deviceProp.minor == 0)
		m_block_size = 320;
	else if (deviceProp.major == 1 && deviceProp.minor == 1)
		m_block_size = 256;
	else if (deviceProp.major == 1 && deviceProp.minor < 4)
		m_block_size = 384;
	else
		{
		cout << "***Warning! Unknown compute " << deviceProp.major << "." << deviceProp.minor << " when tuning block size for GaussianForceGPU" << endl;
		m_block_size = 64;
		}
		
	// ulf workaround setup
	#ifndef DISABLE_ULF_WORKAROUND
	// the ULF workaround is needed on GTX280 and older GPUS
	// it is not needed on C1060, S1070, GTX285, GTX295, and (hopefully) newer ones
	m_ulf_workaround = true;
	
	if (deviceProp.major == 1 && deviceProp.minor >= 3)
		m_ulf_workaround = false;
	if (string(deviceProp.name) == "GTX 280")
		m_ulf_workaround = true;
		
	if (m_ulf_workaround)
		cout << "Notice: ULF bug workaround enabled for GaussianForceGPU" << endl;
	#else
	m_ulf_workaround = false;
	#endif

	// allocate the coeff data on the GPU
	int nbytes = sizeof(float2)*m_pdata->getNTypes()*m_pdata->getNTypes();
	
	d_coeffs.resize(exec_conf.gpu.size());
	for (unsigned int cur_gpu = 0; cur_gpu < exec_conf.gpu.size(); cur_gpu++)
		{
		exec_conf.gpu[cur_gpu]->setTag(__FILE__, __LINE__);
		exec_conf.gpu[cur_gpu]->call(bind(cudaMalloc, (void **)((void *)&d_coeffs[cur_gpu]), nbytes));
		assert(d_coeffs[cur_gpu]);
		exec_conf.gpu[cur_gpu]->call(bind(cudaMemset, (void *)d_coeffs[cur_gpu], 0, nbytes));
		}
	// allocate the coeff data on the CPU
	h_coeffs = new float2[m_pdata->getNTypes()*m_pdata->getNTypes()];
	}
	

GaussianForceGPU::~GaussianForceGPU()
	{
	// free the coefficients on the GPU
	exec_conf.tagAll(__FILE__, __LINE__);
	for (unsigned int cur_gpu = 0; cur_gpu < exec_conf.gpu.size(); cur_gpu++)
		{
		assert(d_coeffs[cur_gpu]);
		exec_conf.gpu[cur_gpu]->call(bind(cudaFree, (void *)d_coeffs[cur_gpu]));
		}
	delete[] h_coeffs;
	}
	
/*! \param block_size Size of the block to run on the device
	Performance of the code may be dependant on the block size run
	on the GPU. \a block_size should be set to be a multiple of 32.
*/
void GaussianForceGPU::setBlockSize(int block_size)
	{
	m_block_size = block_size;
	}

/*! \post The parameters \a epsilon and \a sigma are set for the pairs \a typ1, \a typ2 and \a typ2, \a typ1.
	
	Setting the parameters for typ1,typ2 automatically sets the same parameters for typ2,typ1: there
	is no need to call this funciton for symmetric pairs. Any pairs that this function is not called
	for will have lj1 and lj2 set to 0.0.
	
	\param typ1 Specifies one type of the pair
	\param typ2 Specifies the second type of the pair
	\param epsilon First parameter used to calcluate forces
	\param sigma Second parameter used to calculate forces
*/
void GaussianForceGPU::setParams(unsigned int typ1, unsigned int typ2, Scalar epsilon, Scalar sigma)
	{
	assert(h_coeffs);
	if (typ1 >= m_ntypes || typ2 >= m_ntypes)
		{
		cerr << endl << "***Error! Trying to set Gaussian params for a non existant type! " << typ1 << "," << typ2 << endl << endl;
		throw runtime_error("GaussianForceGPU::setParams argument error");
		}
	
	// set coeffs in both symmetric positions in the matrix
	h_coeffs[typ1*m_pdata->getNTypes() + typ2] = make_float2(epsilon, sigma);
	h_coeffs[typ2*m_pdata->getNTypes() + typ1] = make_float2(epsilon, sigma);
	
	int nbytes = sizeof(float2)*m_pdata->getNTypes()*m_pdata->getNTypes();
	exec_conf.tagAll(__FILE__, __LINE__);
	for (unsigned int cur_gpu = 0; cur_gpu < exec_conf.gpu.size(); cur_gpu++)
		exec_conf.gpu[cur_gpu]->call(bind(cudaMemcpy, d_coeffs[cur_gpu], h_coeffs, nbytes, cudaMemcpyHostToDevice));
	}
		
/*! \post The gaussian forces are computed for the given timestep on the GPU. 
	The neighborlist's compute method is called to ensure that it is up to date
	before forces are computed.
 	\param timestep Current time step of the simulation
 	
 	Calls gpu_compute_gauss_forces to do the dirty work.
*/
void GaussianForceGPU::computeForces(unsigned int timestep)
	{
	// start by updating the neighborlist
	m_nlist->compute(timestep);
	
	// start the profile
	if (m_prof) m_prof->push(exec_conf, "Gauss pair");
	
	// The GPU implementation CANNOT handle a half neighborlist, error out now
	bool third_law = m_nlist->getStorageMode() == NeighborList::half;
	if (third_law)
		{
		cerr << endl << "***Error! GaussianForceGPU cannot handle a half neighborlist" << endl << endl;
		throw runtime_error("Error computing forces in GaussianForceGPU");
		}
	
	// access the neighbor list, which just selects the neighborlist into the device's memory, copying
	// it there if needed
	vector<gpu_nlist_array>& nlist = m_nlist->getListGPU();

	// access the particle data
	vector<gpu_pdata_arrays>& pdata = m_pdata->acquireReadOnlyGPU();
	gpu_boxsize box = m_pdata->getBoxGPU();
	
	// run the kernel on all GPUs in parallel
	exec_conf.tagAll(__FILE__, __LINE__);
	gauss_options opt;
	opt.r_cutsq = m_r_cut * m_r_cut;
	opt.block_size = m_block_size;
	opt.ulf_workaround = m_ulf_workaround;
	opt.shift_mode = m_shift_mode;
	
	for (unsigned int cur_gpu = 0; cur_gpu < exec_conf.gpu.size(); cur_gpu++)
		exec_conf.gpu[cur_gpu]->callAsync(bind(gpu_compute_gauss_forces, m_gpu_forces[cur_gpu].d_data, pdata[cur_gpu], box, nlist[cur_gpu], d_coeffs[cur_gpu], m_pdata->getNTypes(), opt));
	exec_conf.syncAll();
	
	m_pdata->release();
	
	// the force data is now only up to date on the gpu
	m_data_location = gpu;

	Scalar avg_neigh = m_nlist->estimateNNeigh();
	int64_t n_calc = int64_t(avg_neigh * m_pdata->getN());
	int64_t mem_transfer = m_pdata->getN() * (4 + 16 + 20) + n_calc * (4 + 16);
	int64_t flops = n_calc * (3+12+5+6+1+3+7);
	if (m_shift_mode == shift)
		flops += n_calc * 4;
	if (m_prof) m_prof->pop(exec_conf, flops, mem_transfer);
	}

void export_GaussianForceGPU()
	{
	class_<GaussianForceGPU, boost::shared_ptr<GaussianForceGPU>, bases<GaussianForceCompute>, boost::noncopyable >
		("GaussianForceGPU", init< boost::shared_ptr<ParticleData>, boost::shared_ptr<NeighborList>, Scalar >())
		.def("setBlockSize", &GaussianForceGPU::setBlockSize)
		;
	}

#ifdef WIN32
#pragma warning( pop )
#endif
