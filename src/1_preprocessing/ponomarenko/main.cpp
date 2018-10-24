/*
 *
 * This file is part of the Image Processing Framework.
 *
 * Copyright(c) 2011 Miguel Colom.
 * miguel.colom@cmla.ens-cachan.fr
 *
 * This file may be licensed under the terms of of the
 * GNU General Public License Version 2 (the ``GPL'').
 *
 * Software distributed under the License is distributed
 * on an ``AS IS'' basis, WITHOUT WARRANTY OF ANY KIND, either
 * express or implied. See the GPL for the specific language
 * governing rights and limitations.
 *
 * You should have received a copy of the GPL along with this
 * program. If not, go to http://www.gnu.org/licenses/gpl.html
 * or write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>
#include "framework/CFramework.h"
#include "framework/CInspector.h"
#include "algo.h"

#include <time.h>

int main(int argc, char **argv) {
  CFramework *fw = CFramework::get_framework();  
  try {

    // Parse command line
    AlgoOptions opts;
    char input_file[512];
    parseCommandLine(argc, argv, opts, input_file);

    // Load input image
	 {
		 CImage input;
		 input.load((char *)input_file);

		 // Call algorithm
		 int num_channels = input.get_num_channels();
		 int num_bins;
		 float *means, *stds;
		 algorithm(opts, input, means, stds, num_bins);//, means, stds, num_bins);

       // Print results
       for (int bin = 0; bin < num_bins; bin++) {
         // Means
         for (int ch = 0; ch < num_channels; ch++)
             printf("%f  ", means[ch*num_bins+bin]);

         // Standard deviations
         for (int ch = 0; ch < num_channels; ch++)
             printf("%f  ", stds[ch*num_bins+bin]);
         //
         printf("\n");
       }

		 delete[] means;
		 delete[] stds;
	 }
    delete fw;
  }
  catch(...) {    
    printf("A exception raised while running the algorithm.\n");
    printf("Backtrace:\n");
    fw->print_backtrace(stdout);
    // [ToDo]: call user handler for the exception
    exit(-1);
  }
}
