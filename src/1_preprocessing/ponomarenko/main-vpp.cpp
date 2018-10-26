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
extern "C" {
#include "vpp.h"
}

#include <time.h>

int main(int argc, char **argv) {
	CFramework *fw = CFramework::get_framework();
	try {

		// Parse command line
		AlgoOptions opts;
		char input_file[512];
		parseCommandLine(argc, argv, opts, input_file);

		// Init input stream
		int nx, ny, nc;
		FILE* vid_in = vpp_init_input(input_file, &nx, &ny, &nc);
		if (!vid_in) {
			fprintf(stderr, "%s: cannot initialize input '%s'\n", argv[0], input_file);
			return EXIT_FAILURE;
		}

		{
			// Allocate buffer for image
			CImage input(nx, ny, 16, nc);

			// Allocate buffer for vpp
			float *frame = (float *)malloc(nx*ny*nc*sizeof(float));

			int frame_num = 0;
			while (vpp_read_frame(vid_in, frame, nx, ny, nc)) {

				// Copy frame in CImage
				float *channels[nc];
				for (int c = 0; c < nc; ++c)
					channels[c] = input.get_channel(c);

				for (int y = 0, i1 = 0, i2 = 0; y < ny; ++y)
				for (int x = 0                ; x < nx; ++x, ++i1)
				for (int c = 0                ; c < nc; ++c, ++i2)
					channels[c][i1] = frame[i2];

				// Call algorithm
				int num_bins;
				float *means, *stds;
				algorithm(opts, input, means, stds, num_bins);

				// Print results
				printf("frame %d\n", frame_num);
				for (int bin = 0; bin < num_bins; bin++) {
					for (int ch = 0; ch < nc; ch++)
						printf("%f  ", means[ch*num_bins+bin]);
					for (int ch = 0; ch < nc; ch++)
						printf("%f  ", stds[ch*num_bins+bin]);
					printf("\n");
				}

				delete[] means;
				delete[] stds;
				frame_num++;
			}

			free(frame);
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
