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


#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4103 4244 )
#endif

#include <math.h>
#include "HOOMDDumpWriter.h"
#include "HOOMDInitializer.h"
#include "VanHoveAnalyzer.h"
#include "ParticleGroup.h"

#include <iostream>
#include <sstream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
using namespace boost::filesystem;
#include <boost/shared_ptr.hpp>
using namespace boost;

#include <fstream>
using namespace std;

//! Name the unit test module
#define BOOST_TEST_MODULE VanHoveAnalyzerTest
#include "boost_utf_configure.h"

/*! \file test_vanHove.cc
    \brief Unit tests for VanHoveAnaylzer
    \ingroup unit_tests
*/

//! Test basic functionality of HOOMDInitializer
BOOST_AUTO_TEST_CASE( HOOMD_van_hove )
    {
    // create a test input file
    ofstream f("test_input.xml");
    f << "<?xml version =\"1.0\" encoding =\"UTF-8\" ?>\n\
<hoomd_xml version=\"1.3\">\n\
<configuration time_step=\"0\" dimensions=\"3\">\n\
<box lx=\"20\" ly= \"20\" lz=\"20\" xy=\"0\" xz=\"0\" yz=\"0\"/>\n\
<position >\n\
1 1 1\n\
1 -1 1\n\
1 1 -1\n\
-1 1 -1\n\
</position>\n\
<image>\n\
0 0 0\n\
0 0 0\n\
0 0 0\n\
0 0 0\n\
</image>\n\
<mass>\n\
1.0\n\
1.0\n\
1.0\n\
1.0\n\
</mass>\n\
<diameter>\n\
1.0\n\
1.0\n\
1.0\n\
1.0\n\
</diameter>\n\
<type>\n\
0\n\
0\n\
0\n\
0\n\
</type>\n\
</configuration>\n\
</hoomd_xml>" << endl;
    f.close();

    // now that we have created a test file, load it up into a pdata
    boost::shared_ptr<ExecutionConfiguration> exec_conf(new ExecutionConfiguration(ExecutionConfiguration::CPU));
    HOOMDInitializer init(exec_conf,"test_input.xml");
    boost::shared_ptr<SnapshotSystemData> snapshot;
    snapshot = init.getSnapshot();
    boost::shared_ptr<SystemDefinition> sysdef(new SystemDefinition(snapshot));
    boost::shared_ptr<ParticleData> pdata = sysdef->getParticleData();

    unsigned int N_parts = pdata->getN();
    std::vector<unsigned int> member_tags(N_parts);
    for (unsigned int i = 0; i < N_parts; i++)
        {
        member_tags[i] = i;
        }

    boost::shared_ptr<ParticleGroup> all(new ParticleGroup(sysdef, member_tags));
    boost::shared_ptr<VanHoveAnalyzer> van_hove(new VanHoveAnalyzer(
                                                                   sysdef,
                                                                   "test_van_hove.log",
                                                                   10,
                                                                   10,
                                                                   100.0,
                                                                   "",
                                                                   true
                                                                   ));

    van_hove->addColumn(all, "all");
    //Analyze first twn steps (frozen particles)
    for (unsigned int i=1; i<=10; i++)
        {
        van_hove->analyze(i);
        }
    {
    //Analyze first twn steps (frozen particles)
    ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::read);
    for (int i = 0; i < 4; i++)
        {
        h_pos.data[i].x += 1;
        }
    }

    // verify all parameters
    for (unsigned int i=11; i<=21; i++)
        {
        van_hove->analyze(i);
        }

    f.close();
    //check line by line output. Should have a blip at r =1 for the first 10 lines (steps 11-20), then a blip at zero for step 21
    ifstream f_in("test_van_hove.log");
    string line;
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "timestep\tnum_bins\tr_cut\tall");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "11\t10\t10\t0\t0.03410463066\t0\t0\t0\t0\t0\t0\t0\t0");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "12\t10\t10\t0\t0.03410463066\t0\t0\t0\t0\t0\t0\t0\t0");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "13\t10\t10\t0\t0.03410463066\t0\t0\t0\t0\t0\t0\t0\t0");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "14\t10\t10\t0\t0.03410463066\t0\t0\t0\t0\t0\t0\t0\t0");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "15\t10\t10\t0\t0.03410463066\t0\t0\t0\t0\t0\t0\t0\t0");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "16\t10\t10\t0\t0.03410463066\t0\t0\t0\t0\t0\t0\t0\t0");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "17\t10\t10\t0\t0.03410463066\t0\t0\t0\t0\t0\t0\t0\t0");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "18\t10\t10\t0\t0.03410463066\t0\t0\t0\t0\t0\t0\t0\t0");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "19\t10\t10\t0\t0.03410463066\t0\t0\t0\t0\t0\t0\t0\t0");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "20\t10\t10\t0\t0.03410463066\t0\t0\t0\t0\t0\t0\t0\t0");
    getline(f_in, line);
    BOOST_CHECK_EQUAL(line, "21\t10\t10\t0.2387324146\t0\t0\t0\t0\t0\t0\t0\t0\t0");
    BOOST_REQUIRE(!f_in.bad());
    f_in.close();
    // clean up after ourselves
    remove_all("test_input.xml");
    //remove_all("test_van_hove.log");
    }

#ifdef WIN32
#pragma warning( pop )
#endif
