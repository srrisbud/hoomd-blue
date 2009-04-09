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
// $URL$

/*! \file NVTUpdater.h
	\brief Declares the NVTUpdater class
*/

#include "Updater.h"
#include "Integrator.h"
#include "Variant.h"
#include <vector>
#include <boost/shared_ptr.hpp>

#ifndef __NVTUPDATER_H__
#define __NVTUPDATER_H__

//! NVT Integration via the Nose-Hoover thermostat
/*! This updater performes constant N, constant volume, constant temperature (NVT) dynamics. Particle positions and velocities are 
	updated according to the Nose-Hoover algorithm. The forces that drive this motion are defined external to this class
	in ForceCompute. Any number of ForceComputes can be given, the resulting forces will be summed to produce a net force on 
	each particle.
	
	\ingroup updaters
*/
class NVTUpdater : public Integrator
	{
	public:
		//! Constructor
		NVTUpdater(boost::shared_ptr<SystemDefinition> sysdef, Scalar deltaT, Scalar tau, boost::shared_ptr<Variant> T);
		
		//! Take one timestep forward
		virtual void update(unsigned int timestep);
		
		//! Update the temperature
		/*! \param T New temperature to set
		*/
		virtual void setT(boost::shared_ptr<Variant> T) { m_T = T; }
				
		//! Update the tau value
		/*! \param tau New time constant to set
		*/		
		virtual void setTau(Scalar tau) { m_tau = tau; }
		
		//! Returns a list of log quantities this compute calculates
		virtual std::vector< std::string > getProvidedLogQuantities();
		
		//! Calculates the requested log value and returns it
		virtual Scalar getLogValue(const std::string& quantity, unsigned int timestep);		
	protected:
		Scalar m_tau;					//!< tau value for Nose-Hoover
		boost::shared_ptr<Variant> m_T;	//!< Temperature set point
		Scalar m_Xi;					//!< Friction coeff
		Scalar m_eta;					//!< Added degree of freedom
		bool m_accel_set;				//!< Flag to tell if we have set the accelleration yet
		Scalar m_curr_T;				//!< Current calculated temperature of the system
	};
	
//! Exports the NVTUpdater class to python
void export_NVTUpdater();

#endif
