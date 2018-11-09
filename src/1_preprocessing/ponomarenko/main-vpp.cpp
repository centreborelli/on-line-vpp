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
		char output_file[512];
		parseCommandLine(argc, argv, opts, input_file, output_file);

		// Init input stream
		int nx, ny, nc, nb;
		FILE* vid_in = vpp_init_input(input_file, &nx, &ny, &nc);
		if (!vid_in) {
			fprintf(stderr, "%s: cannot initialize input '%s'\n", argv[0], input_file);
			return EXIT_FAILURE;
		}

		// Set number of bins
		if (opts.num_bins <= 0) opts.num_bins = nx * ny / 42000;
		if (opts.num_bins <= 0) opts.num_bins = 1; // Force at least one bin
		nb = opts.num_bins;

		// Init output stream
		FILE* sigma_out = vpp_init_output(output_file, 2*nc, nb, 1);
		if (!sigma_out)
		{
			fprintf(stderr, "%s: cannot initialize output '%s'\n", argv[0], output_file);
			return EXIT_FAILURE;
		}

		{
			// Allocate buffer for image
			CImage input(nx, ny, 16, nc);

			// Allocate buffer for vpp
			float *frame = (float *)malloc(nx*ny*nc*sizeof(float));
			float *sigma_table = (float *)malloc(nb*2*nc*sizeof(float));

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

				// Write results on output buffer
				for (int bin = 0, i = 0; bin < nb; bin++) {
					for (int ch = 0; ch < nc; ch++, i++)
						sigma_table[i] = means[ch*num_bins+bin];
					for (int ch = 0; ch < nc; ch++, i++)
						sigma_table[i] = stds[ch*num_bins+bin];
				}

				//save the optical flow
				vpp_write_frame(sigma_out, sigma_table, 2*nc, nb, 1);

				delete[] means;
				delete[] stds;
				frame_num++;
			}

			free(frame);
			free(sigma_table);
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
