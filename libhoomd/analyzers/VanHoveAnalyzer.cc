/*
Highly Optimized Object-oriented Many-particle Dynamics -- Blue Edition
(HOOMD-blue) Open Source Software License Copyright 2009-2014 The Regents of
the University of Michigan All rights reserved.

HOOMD-blue may contain modifications ("Contributions") provided, and to which
copyright is held, by various Contributors who have granted The Regents of the
University of Michigan the right to modify and/or distribute such Contributions.

You may redistribute, use, and create derivate works of HOOMD-blue, in source
and binary forms, provided you abide by the following conditions:

* Redistributions of source code must retain the above copyright notice, this
list of conditions, and the following disclaimer both in the code and
prominently in any materials provided with the distribution.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions, and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* All publications and presentations based on HOOMD-blue, including any reports
or published results obtained, in whole or in part, with HOOMD-blue, will
acknowledge its use according to the terms posted at the time of submission on:
http://codeblue.umich.edu/hoomd-blue/citations.html

* Any electronic documents citing HOOMD-Blue will link to the HOOMD-Blue website:
http://codeblue.umich.edu/hoomd-blue/

* Apart from the above required attributions, neither the name of the copyright
holder nor the names of HOOMD-blue's contributors may be used to endorse or
promote products derived from this software without specific prior written
permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR ANY
WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// Maintainer: joaander

/*! \file VanHoveAnalyzer.cc
    \brief Defines the VanHoveAnalyzer class
*/

#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4244 )
#endif

#include "VanHoveAnalyzer.h"
#include "HOOMDInitializer.h"

#ifdef ENABLE_MPI
#include "Communicator.h"
#endif

#include <boost/python.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
using namespace boost::python;
using namespace boost::filesystem;

#include <iomanip>
using namespace std;

/*! \param sysdef SystemDefinition containing the Particle data to analyze
  \param fname File name to write output to
  \param header_prefix String to print before the file header
  \param overwrite Will overwite an exiting file if true (default is to append)

  On construction, the initial coordinates of all parrticles in the system are recoreded. The file is opened
  (and overwritten if told to). Nothing is initially written to the file, that will occur on the first call to
  analyze()
*/
VanHoveAnalyzer::VanHoveAnalyzer(boost::shared_ptr<SystemDefinition> sysdef,
                                 std::string fname,
                                 const int num_windows,
                                 const int num_bins,
                                 const Scalar r2max,
                                 const std::string& header_prefix,
                                 bool overwrite)
    : Analyzer(sysdef), m_delimiter("\t"), m_header_prefix(header_prefix), m_appending(false),
      m_columns_changed(false), m_num_bins(num_bins), m_num_windows(num_windows), m_r2max(r2max)
    {
    m_exec_conf->msg->notice(5) << "Constructing VanHoveAnalyzer: " << fname << " " << header_prefix << " " << overwrite << endl;

    SnapshotParticleData snapshot(m_pdata->getNGlobal());

    m_pdata->takeSnapshot(snapshot);

#ifdef ENABLE_MPI
    // if we are not the root processor, do not perform file I/O
    if (m_comm && !m_exec_conf->isRoot())
        {
        return;
        }
#endif

    // open the file
    if (exists(fname) && !overwrite)
        {
        m_exec_conf->msg->notice(3) << "analyze.van_hove: Appending van_hove to existing file \"" << fname << "\"" << endl;
        m_file.open(fname.c_str(), ios_base::in | ios_base::out | ios_base::ate);
        m_appending = true;
        }
    else
        {
        m_exec_conf->msg->notice(3) << "analyze.van_hove: Creating new van_hove in file \"" << fname << "\"" << endl;
        m_file.open(fname.c_str(), ios_base::out);
        }

    if (!m_file.good())
        {
        m_exec_conf->msg->error() << "analyze.van_hove: Unable to open file " << fname << endl;
        throw runtime_error("Error initializing analyze.van_hove");
        }

    // resize the histogram
    m_van_hove.resize(m_num_bins);

    // zero counters/pointers
    m_R0_offset = 0;
    m_num_samples = 0;

    // record the initial particle positions by tag
    m_initial_x.resize(m_pdata->getNGlobal() * m_num_windows);
    m_initial_y.resize(m_pdata->getNGlobal() * m_num_windows);
    m_initial_z.resize(m_pdata->getNGlobal() * m_num_windows);
    BoxDim box = m_pdata->getGlobalBox();

    // for each particle in the data
    for (unsigned int tag = 0; tag < snapshot.size; tag++)
        {
        // save its initial position
        Scalar3 pos = snapshot.pos[tag];
        Scalar3 unwrapped = box.shift(pos, snapshot.image[tag]);
        m_initial_x[tag] = unwrapped.x;
        m_initial_y[tag] = unwrapped.y;
        m_initial_z[tag] = unwrapped.z;
        }
    }

VanHoveAnalyzer::~VanHoveAnalyzer()
    {
    m_exec_conf->msg->notice(5) << "Destroying VanHoveAnalyzer" << endl;
    }

/*!\param timestep Current time step of the simulation

    analyze() will first write out the file header if the columns have changed.

    On every call, analyze() will write calculate the VanHove for each group and write out a row in the file.
*/
void VanHoveAnalyzer::analyze(unsigned int timestep)
    {
    if (m_prof)
        m_prof->push("Analyze VanHove");

    // take particle data snapshot
    SnapshotParticleData snapshot(m_pdata->getNGlobal());

    m_pdata->takeSnapshot(snapshot);

#ifdef ENABLE_MPI
    // if we are not the root processor, do not perform file I/O
    if (m_comm && !m_exec_conf->isRoot())
        {
        if (m_prof) m_prof->pop();
        return;
        }
#endif

    // error check
    if (m_columns.size() == 0)
        {
        m_exec_conf->msg->warning() << "analyze.van_hove: No columns specified in the VanHove analysis" << endl;
        return;
        }

    // ignore writing the header on the first call when appending the file
    if (m_columns_changed && m_appending)
        {
        m_appending = false;
        m_columns_changed = false;
        }

    // write out the header only once if the columns change
    if (m_columns_changed)
        {
        writeHeader();
        m_columns_changed = false;
        }

    // write out the row every time
    if (m_R0_offset >= m_num_windows)
        {
        writeRow(timestep, snapshot);
        m_num_samples += 1;
        }

    if (m_prof)
        m_prof->pop();

    rollR0(snapshot);
    }

/*! \param delimiter New delimiter to set

    The delimiter is printed between every element in the row of the output
*/
void VanHoveAnalyzer::setDelimiter(const std::string& delimiter)
    {
    m_delimiter = delimiter;
    }

/*! \param group Particle group to calculate the VanHove of
    \param name Name to print in the header of the file

    After a column is added with addColumn(), future calls to analyze() will calculate the VanHove of the particles defined
    in \a group and print out an entry under the \a name header in the file.
*/
void VanHoveAnalyzer::addColumn(boost::shared_ptr<ParticleGroup> group, const std::string& name)
    {
    m_columns.push_back(column(group, name));
    m_columns_changed = true;
    m_histograms.resize(m_columns.size());
    for (unsigned int i=0; i<m_columns.size(); i++)
        m_histograms[i].resize(m_num_bins);
    }

void VanHoveAnalyzer::rollR0(const SnapshotParticleData& snapshot)
    {
    // for each particle in the data
    unsigned int N_particles = snapshot.size;
    BoxDim box = m_pdata->getGlobalBox();
    unsigned int offset = (m_R0_offset % m_num_windows) * N_particles;
    for (unsigned int tag = 0; tag < N_particles; tag++)
        {
        // save its initial position
        Scalar3 pos = snapshot.pos[tag];
        Scalar3 unwrapped = box.shift(pos, snapshot.image[tag]);
        m_initial_x[offset + tag] = unwrapped.x;
        m_initial_y[offset + tag] = unwrapped.y;
        m_initial_z[offset + tag] = unwrapped.z;
        }
    m_R0_offset += 1;
    }
/*! The entire header row is written to the file. First, timestep is written as every file includes it and then the
    columns are looped through and their names printed, separated by the delimiter.
*/
void VanHoveAnalyzer::writeHeader()
    {
    // write out the header prefix
    m_file << m_header_prefix;

    // timestep is always output
    m_file << "timestep"<< m_delimiter
           << "num_bins" << m_delimiter
           << "r_cut";

    if (m_columns.size() == 0)
        {
        m_exec_conf->msg->warning() << "analyze.van_hove: No columns specified in the VanHove analysis" << endl;
        return;
        }

    // only print the delimiter after the timestep if there are more columns
    m_file  << m_delimiter;

    // write all but the last of the quantities separated by the delimiter
    for (unsigned int i = 0; i < m_columns.size()-1; i++)
        m_file << m_columns[i].m_name << m_delimiter;
    // write the last one with no delimiter after it
    m_file << m_columns[m_columns.size()-1].m_name << endl;
    m_file.flush();
    }

/*! \param group Particle group to calculate the VanHove of
    Loop through all particles in the given group and calculate the VanHove over them.
    \returns The calculated VanHove
*/
void VanHoveAnalyzer::calcVanHove(boost::shared_ptr<ParticleGroup const> group,
                                  const SnapshotParticleData& snapshot,
                                  const unsigned int group_index)
    {

    std::vector<Scalar>& g_histogram = m_histograms[group_index];
    //clear the existing histogram
    if (m_num_samples == 0)
        {
        memset(&g_histogram[0], 0, sizeof(Scalar) * m_num_bins);
        }

    BoxDim box = m_pdata->getGlobalBox();

    // handle the case where there are 0 members gracefully
    if (group->getNumMembersGlobal() == 0)
        {
        m_exec_conf->msg->warning() << "analyze.van_hove: Group has 0 members, reporting a calculated van_hove of 0.0" << endl;
        return;
        }


    // for each particle in the group
    // the oldest set of R0 coordinates is one passed the current position of the m_R0_offset pointer
    unsigned int offset = (m_R0_offset) % (m_num_windows) * snapshot.size;
    for (unsigned int group_idx = 0; group_idx < group->getNumMembersGlobal(); group_idx++)
        {
        // get the tag for the current group member from the group
        unsigned int tag = group->getMemberTag(group_idx);
        Scalar3 pos = snapshot.pos[tag];
        int3 image = snapshot.image[tag];
        Scalar3 unwrapped = box.shift(pos, image);
        Scalar dx = unwrapped.x - m_initial_x[offset + tag];
        Scalar dy = unwrapped.y - m_initial_y[offset + tag];
        Scalar dz = unwrapped.z - m_initial_z[offset + tag];
        Scalar dr2 = dx*dx + dy*dy + dz*dz;
        if (dr2 < m_r2max)
            {
            g_histogram[(int)(m_num_bins * sqrt(dr2/m_r2max))] += 1;
            }
        }
    // 4\pi/3 dr^3 * N_particles (in group) * N_samples (in rolling average)
    Scalar dr = sqrt(m_r2max)/m_num_bins;
    Scalar dV = 0;
    Scalar denom_precompute  =  M_PI * 4.0 * group->getNumMembersGlobal() * (m_num_samples + 1) * dr * dr * dr / 3.0;

    // divide to complete the average
    for (unsigned int i=0; i< m_num_bins; i++)
        {
        dV =(3 * i * (i + 1) + 1);
        m_van_hove[i] = g_histogram[i]/(denom_precompute * dV);
        }
    }

/*! \param timestep current time step of the simulation

    Performs all the steps needed in order to calculate the VanHoves for all the groups in the columns and writes out an
    entire row to the file.
*/
void VanHoveAnalyzer::writeRow(unsigned int timestep, const SnapshotParticleData& snapshot)
    {
    if (m_prof) m_prof->push("VanHove");

    // The timestep is always output
    m_file << setprecision(10) << timestep
           << m_delimiter << m_num_bins
           << m_delimiter << sqrt(m_r2max);

    // quit now if there is nothing to log
    if (m_columns.size() == 0)
        {
        return;
        }

    // only print the delimiter after the timestep if there are more columns

    // write all but the last of the columns separated by the delimiter
    for (unsigned int i = 0; i < m_columns.size(); i++)
        {
        calcVanHove(m_columns[i].m_group, snapshot, i);
        m_file << setprecision(10);
        for (unsigned int j = 0; j < m_num_bins; j++)
            {
            m_file << m_delimiter << m_van_hove[j];
            }
        }
    m_file << std::endl;
    m_file.flush();

    if (!m_file.good())
        {
        m_exec_conf->msg->error() << "analyze.van_hove: I/O error while writing file" << endl;
        throw runtime_error("Error writting van_hove file");
        }

    if (m_prof) m_prof->pop();
    }

void export_VanHoveAnalyzer()
    {
    class_<VanHoveAnalyzer, boost::shared_ptr<VanHoveAnalyzer>, bases<Analyzer>, boost::noncopyable>
    ("VanHoveAnalyzer", init< boost::shared_ptr<SystemDefinition>, const std::string&,
                              const unsigned int, const unsigned int, const Scalar, 
                              const std::string&, bool >())
    .def("setDelimiter", &VanHoveAnalyzer::setDelimiter)
    .def("addColumn", &VanHoveAnalyzer::addColumn)
    ;
    }

#ifdef WIN32
#pragma warning( pop )
#endif
